#!/usr/bin/env python3
"""Export corrected thesis N2LOopt to NuHamil's unweighted radial binary.

Run `make` first. Python dependencies: numpy and scipy. No Coulomb is added.
The symmetric choice is a specified symmetric limit, not a refitted force.
"""
from pathlib import Path
import argparse
import ctypes as ct
import hashlib
import json
import struct
import time
import numpy as np
from projector_math import channel_weights

HC = 197.3269804  # MeV fm; explicitly recorded, not fitted to a reference.
ROOT = Path(__file__).resolve().parent


def read_binary(path):
    raw=Path(path).read_bytes()
    n,jmax,nch=struct.unpack_from('<3i',raw)
    p=np.frombuffer(raw,'<f8',n,12).copy()
    w=np.frombuffer(raw,'<f8',n,12+8*n).copy()
    offset=12+16*n; channels={}
    for _ in range(nch):
        j,par,s,z,dim=struct.unpack_from('<5i',raw,offset);offset+=20
        if dim not in (n,2*n): raise ValueError('Invalid block dimension')
        channels[j,par,s,z]=np.frombuffer(raw,'<f8',dim*dim,offset).reshape((dim,dim),order='F').copy()
        offset+=8*dim*dim
    if offset!=len(raw): raise ValueError('Unexpected binary length')
    return p,w,jmax,channels


def channels_for(jmax,n):
    result={}
    for j in range(jmax+1):
      for s in (0,1):
        groups=[[j]] if s==0 else (([[j]] if j>0 else [])+[[1] if j==0 else [j-1,j+1]])
        for ls in groups:
          t=(ls[0]+s+1)%2
          for z in (-1,0,1):
            if z and t!=1: continue
            result[j,(-1)**ls[0],s,z]=(ls,np.zeros((len(ls)*n,len(ls)*n)))
    return dict(sorted(result.items()))


def generate(variant,outdir,n=100,kmax=8.,jmax=8,nangle=128):
    lib=ct.CDLL(str(ROOT/'libn2loopt.so'))
    lib.source2_initialize.argtypes=[ct.c_int]
    lib.source2_kernel.argtypes=[ct.c_double,ct.c_double,ct.c_int,ct.c_void_p,ct.c_void_p]
    lib.source2_initialize(variant)
    xp,w=np.polynomial.legendre.leggauss(n)
    p=(xp+1)*kmax/2;w=w*kmax/2
    x,qw=np.polynomial.legendre.leggauss(nangle)
    weights=channel_weights(jmax,x,qw)
    channels=channels_for(jmax,n)
    slots=[]; tensors=[]
    # A single dense contraction projects all channels for each momentum pair.
    for key,(ls,mat) in channels.items():
      j,par,s,z=key;t=(ls[0]+s+1)%2
      iz={-1:0,1:3}.get(z,1 if t==1 else 2)
      for io,lo in enumerate(ls):
        for ii,li in enumerate(ls):
          ww=np.zeros((nangle,4,4,4),order='F')
          sl=slice(0,1) if s==0 else slice(1,4)
          ww[:,sl,sl,iz]=weights[j,s,lo,li]*(-1.)**((lo-li)//2)*HC**3/(2*np.pi)**3
          tensors.append(ww.ravel(order='F'))
          slots.append((mat,io*n,ii*n))
    tensor=np.array(tensors)
    kernel=np.empty((nangle,4,4,4),order='F')
    start=time.time()
    for i,po in enumerate(p):
      for k,pi in enumerate(p):
        lib.source2_kernel(po*HC,pi*HC,nangle,x.ctypes.data,kernel.ctypes.data)
        # einsum avoids repeated multithreaded BLAS startup for tiny contractions.
        vals=np.einsum('ij,j->i',tensor,kernel.ravel(order='F'),optimize=False)
        for val,(mat,ro,co) in zip(vals,slots):mat[ro+i,co+k]=val
      if (i+1)%20==0: print(f'variant={variant} rows={i+1}/{n}, seconds={time.time()-start:.1f}',flush=True)
    mats={key:value[1] for key,value in channels.items()}
    if not all(np.isfinite(m).all() for m in mats.values()): raise ValueError('Nonfinite potential')
    herm=max(float(np.max(abs(m-m.T))) for m in mats.values())
    scale=max(float(np.max(abs(m))) for m in mats.values())
    charge=0.; count=0
    for key,m in mats.items():
      j,par,s,z=key
      if z!=-1:continue
      for zz in (0,1):
        charge=max(charge,float(np.max(abs(m-mats[j,par,s,zz]))));count+=m.size
    if herm>1e-10*max(scale,1.): raise ValueError(f'Hermiticity failed: {herm}')
    if variant==1 and charge>1e-12*max(scale,1.): raise ValueError(f'Isospin equality failed: {charge}')
    label='N2LOopt_source2'+('_IS' if variant else '')
    output=outdir/f'{label}_kmax{kmax:g}_N{n}_Jmax{jmax}.bin'
    with output.open('wb') as f:
      f.write(struct.pack('<3i',n,jmax,len(mats)))
      f.write(p.astype('<f8').tobytes());f.write(w.astype('<f8').tobytes())
      for key,m in mats.items():
        f.write(struct.pack('<5i',*key,len(m)))
        f.write(np.asarray(m,dtype='<f8').tobytes(order='F'))
    rp,rw,rj,rm=read_binary(output)
    assert np.array_equal(rp,p) and np.array_equal(rw,w) and rj==jmax
    assert all(np.array_equal(m,rm[key]) for key,m in mats.items())
    report=dict(file=str(output),sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
        bytes=output.stat().st_size,variant=variant,nmesh=n,kmax_fm_inv=kmax,jmax=jmax,
        angular_points=nangle,channels=len(mats),hbarc_MeV_fm=HC,
        units='momenta and weights fm^-1; potential MeV fm^3; unweighted',
        radial_phase='i^(Lin-Lout)',max_abs_potential=scale,
        max_hermiticity_residual=herm,max_T1_charge_difference=charge,
        T1_charge_entry_comparisons=count,seconds=time.time()-start,
        coulomb_included=False,
        generator_revision='source2-channel-nucleon-masses-2026-11-18',
        two_pion_inverse_mass_convention=('common nucleon mass' if variant else 'charge-channel nucleon mass'),
        nucleon_masses_MeV=({'pp':938.91865,'pn':938.91865,'nn':938.91865} if variant else
                            {'pp':938.272,'pn':938.91865,'nn':939.5653}),
        library_sha256=hashlib.sha256((ROOT/'libn2loopt.so').read_bytes()).hexdigest(),
        source_commit='c2d6a2fcbaa84e06de4141b8f33d2fe2dff01d56',
        symmetric_prescription=(dict(nucleon_mass_MeV=938.91865,pion_mass_MeV=138.039,
          singlet_contact_1e4_GeV_minus2=-0.151757290446826,
          description='Arithmetic-average masses/contact; no refit') if variant else None))
    output.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outdir',type=Path,required=True)
    ap.add_argument('--variant',choices=['original','symmetric','both'],default='both')
    ap.add_argument('--nmesh',type=int,default=100)
    ap.add_argument('--kmax',type=float,default=8.)
    ap.add_argument('--jmax',type=int,default=8)
    ap.add_argument('--nangle',type=int,default=128)
    args=ap.parse_args()
    if not (2<=args.nmesh<=1000 and 0<=args.jmax<=20 and args.kmax>0 and args.nangle>=32):ap.error('Invalid grid')
    args.outdir.mkdir(parents=True,exist_ok=True)
    variants={'original':[0],'symmetric':[1],'both':[0,1]}[args.variant]
    for variant in variants:
      print(json.dumps(generate(variant,args.outdir,args.nmesh,args.kmax,args.jmax,args.nangle),indent=2),flush=True)


if __name__=='__main__':main()
