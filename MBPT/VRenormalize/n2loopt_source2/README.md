# N2LOopt from the CompPhysics thesis source

This package evaluates the two-nucleon N2LOopt interaction from the Sean Miller
CCD prototype in `CompPhysics/ThesisProjects`, projects it into LSJ partial
waves, and exports NuHamil momentum-space input. It also provides the CENS
`/cpot/` interface. No Coulomb or three-nucleon interaction is included.

Upstream source:
https://github.com/CompPhysics/ThesisProjects/tree/c2d6a2fcbaa84e06de4141b8f33d2fe2dff01d56/doc/MSc/msc_students/former/SeanMiller/Thesis/Codes/CCD_prototype/CCD/Systems/NLO2opt

The unmodified Fortran files are preserved in `upstream/`; `vendor/` contains
the compiled versions. `source2-corrections.patch` records every change.
This is a local integration, not an upstream release or endorsement.

## Corrections required before use

1. Remove bare `WT` from `vnlo`: one-pion exchange already appears once as
   `WT*(sigma1.q)*(sigma2.q)` in `vlo`. The extra coefficient had wrong units
   and generated forbidden spin transitions at zero momentum.
2. Return a defined zero when the exponential regulator underflows. Remove
   the artificial `1e-10` regulator cutoff, which previously returned an
   uninitialized function value.
3. Use analytic zero-transfer limits of the L and A loop functions rather
   than modifying the input transfer to `1e-8`.
4. Rename `ang_mom_functions` to `source2_ang_mom_functions` to avoid a CENS
   module-name collision. This has no physical effect.

The wrapper removes the original finite-box momentum/energy scales. The
angular projector uses Condon-Shortley spin spherical harmonics, then applies
the radial phase `i^(Lin-Lout)` and `(2*pi)^(-3)`. Output is multiplied by
`hbarc^3`, with `hbarc=197.3269804 MeV fm`. No fitted normalization is used.

## Build and generate momentum-space input

Requires GNU Fortran, make, Python 3, NumPy and SciPy (tested with GNU 11.4,
NumPy 1.26.4, SciPy 1.15.3 in WSL Ubuntu).

```sh
cd /mnt/d/Code/CENS/MBPT/VRenormalize/n2loopt_source2
make -j4
python3 export_nuhamil.py --variant both --nmesh 100 --kmax 8 --jmax 8 --outdir /path/to/results
```

`--variant original` uses charge-channel nucleon masses in the two-pion
1/m terms as well as the relativistic prefactor: pp = 938.272 MeV,
pn = 938.91865 MeV (the source2 arithmetic average), nn = 939.5653 MeV.
This convention was applied on 2026-11-18 after benchmarking against the
standard NuHamil N2LOopt input. Charge-dependent singlet contacts and the
source2 high-precision LECs are preserved, with the bug fixes above.
`--variant symmetric` uses the explicit symmetric-limit prescription below.
Each binary has a JSON manifest containing hashes, units and matrix checks.
The momentum generator is independent of `renorm.ini`; it outputs the bare
potential, with no Vlowk or SRG evolution.

The validated configuration is N=100, kmax=8 fm^-1, Jmax=8 and 128 angular
points. Larger J or different grids should be convergence-tested. The
Fortran CENS interface accepts J=0..20 and uses 128 angular points. Its
module state and the upstream routine are serial interfaces, not thread-safe.

## Specified isospin-symmetric limit

Variant 1 sets all nucleon masses in the generator to 938.91865 MeV, all
pion masses to 138.039 MeV, and all charge-channel singlet LO contacts to
-0.151757290446826 in units of 10^4 GeV^-2. These are arithmetic averages of
the source2 values. Derived mass arrays are recomputed. The physical variant's channel-mass
override is disabled in this variant, so the common nucleon mass is used
throughout the generator. The triplet contact,
other contacts, c_i, cutoff and regulator are unchanged.

This removes the generator's charge dependence. It is not a new fit to NN
data and is not a unique definition of an isospin-symmetric N2LOopt. T=0
matrix elements also change because the pion and nucleon masses change.
The output is named `N2LOopt_source2_IS` to distinguish this prescription.

## CENS options

The parent VRenormalize build recognizes:

```
type_of_pot = n2loopt
csb_choice = csb
cib_choice = cib
coulomb_included = no-coulomb
```

For the symmetric potential use `type_of_pot = n2loopt-isospin`. Its symmetry
comes from the generator, not replacement of charge channels. The original
option retains the existing CENS no-csb/no-cib channel-substitution behavior.

CENS's older oscillator/Vlowk machinery still uses its separate, unequal
external nucleon-mass constants. Therefore choosing a symmetric potential
does not by itself make the entire legacy CENS Hamiltonian symmetric. Use
the momentum exporter for the verified symmetric NuHamil input.

## NuHamil

Set `input_nn_file` to the generated binary, `nmesh2=100`, `pmax2=8.0`, and
`J2max_NNint=8`. Use a distinct `NNInt` tag to avoid reusing a cached older
interaction. For symmetric calculations explicitly set `coul=.false.` and
`pn_same_mass=.true.`. The provided example inputs use these settings.

Format: little-endian int32 `(NMesh,Jmax,nchan)`, float64 momentum and weight
arrays in fm^-1, then each channel's int32 `(J,parity,S,Tz,ndim)` and a
Fortran-column-major matrix in MeV fm^3. Tz=-1/0/+1 denotes pp/pn/nn. Coupled
states are ordered L=J-1 then J+1; J=0 triplet contains only L=1. Matrix
elements are unweighted. N100/J8 has 54 channels and 8,162,692 bytes.

## Validation and limits

Validation includes analytic angular-projection checks, independent Python
and Fortran projectors, rotation covariance, Hermiticity, absence of T=0/T=1
mixing, all exported T=1 charge-channel entries, angular convergence, CENS
build/dispatch smoke tests, and actual NuHamil read/write and oscillator
conversion tests. Details are in the job output validation reports.

The corrected kernel agrees with an independent legacy N2LOopt routine at
about 5.4e-13 relative to the largest tested amplitude when its constants and
mass prescription are explicitly matched. An earlier source2 build used
the average nucleon mass in all two-pion 1/m terms; its largest HO discrepancy
from the standard file was 6.697 keV. Applying charge-channel nucleon masses
reduced that maximum to 16.21 eV for hw=20 MeV, emax=6, e2max=12, with Coulomb.
The standard momentum binary itself reproduces the supplied HO file exactly.
Small residual differences remain, so the updated source2 output is not
bitwise identical to the standard NuHamil N2LO_opt binary. The symmetric
binary is unchanged by the channel-mass update.
