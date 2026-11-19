"""LSJ angular projection of a rotationally scalar two-nucleon kernel.

The incoming relative momentum is along +z; the outgoing vector lies in
the xz plane, with cos(theta)=x.  Spin axes are the same fixed Cartesian
axes on both sides. Kernel axes are (..., theta, m_out, m_in), with each
spin projection in ascending order -S,...,+S. Spherical harmonics and
Clebsch-Gordan coefficients use the Condon-Shortley convention.

This module supplies angular integrals only. In particular it does NOT
include (2*pi)^-3, hbar*c, antisymmetrization factors, or radial weights.
An angle-independent spin-independent kernel C projects to 4*pi*C in an
S wave. A central kernel F(x) projects to 2*pi*integral(P_l(x)*F(x), dx).

These are raw momentum angular states. The real radial-wave convention
used by the source's named partial-wave contact LECs and NuHamil adds
i**(Lin-Lout), equal to (-1)**((Lout-Lin)//2) in same-parity channels.
That phase is deliberately left to the adapter, so it is applied once.
"""

from functools import lru_cache
from math import factorial, pi, sqrt

import numpy as np
from scipy.special import eval_legendre, lpmv


@lru_cache(maxsize=None)
def cg(j1, m1, j2, m2, j, m):
    """Real Clebsch-Gordan coefficient; all arguments are integers here."""
    if (m1 + m2 != m or abs(m1) > j1 or abs(m2) > j2 or abs(m) > j
            or j < abs(j1-j2) or j > j1+j2 or min(j1, j2, j) < 0):
        return 0.0
    f = factorial
    prefactor = sqrt((2*j+1) * f(j+j1-j2) * f(j-j1+j2)
                     * f(j1+j2-j) / f(j1+j2+j+1))
    prefactor *= sqrt(f(j+m)*f(j-m)*f(j1-m1)*f(j1+m1)
                      * f(j2-m2)*f(j2+m2))
    low = max(0, j2-j-m1, j1+m2-j)
    high = min(j1+j2-j, j1-m1, j2+m2)
    series = sum((-1)**k / (f(k)*f(j1+j2-j-k)*f(j1-m1-k)
                             * f(j2+m2-k)*f(j-j2+m1+k)
                             * f(j-j1-m2+k))
                 for k in range(low, high+1))
    return prefactor * series


def real_y_at_phi_zero(l, m, x):
    """Y_lm(theta,phi=0), including the Condon-Shortley phase."""
    x = np.asarray(x, dtype=float)
    if abs(m) > l:
        return np.zeros_like(x)
    a = abs(m)
    result = sqrt((2*l+1)/(4*pi)*factorial(l-a)/factorial(l+a)) * lpmv(a, l, x)
    return result * ((-1)**a if m < 0 else 1)


def angular_weights(j, s, lout, lin, x, quadrature_weights):
    """Return weights[theta,m_out+S,m_in+S], including dx quadrature.

    Derived by averaging the LSJ matrix element over magnetic J and then
    fixing the incoming direction to +z. This yields the factor
    8*pi**2/(2*J+1)*sqrt((2*Lin+1)/(4*pi)). No factors depending on either
    radial momentum are present.
    """
    if s not in (0, 1):
        raise ValueError("Only nucleon-pair spins S=0,1 are supported")
    if any(l < 0 or j < abs(l-s) or j > l+s for l in (lin, lout)):
        raise ValueError("Invalid LSJ triangle")
    x = np.asarray(x, dtype=float)
    qw = np.asarray(quadrature_weights, dtype=float)
    if x.ndim != 1 or qw.shape != x.shape:
        raise ValueError("x and quadrature weights must be matching 1D arrays")
    result = np.zeros((len(x), 2*s+1, 2*s+1))
    prefactor = 8*pi*pi/(2*j+1) * sqrt((2*lin+1)/(4*pi))
    for mi in range(-s, s+1):
        # Incoming orbital magnetic projection is zero, hence M_J = mi.
        incoming = cg(lin, 0, s, mi, j, mi)
        if incoming == 0:
            continue
        for mo in range(-s, s+1):
            ml = mi-mo
            result[:, mo+s, mi+s] = (
                prefactor * incoming * cg(lout, ml, s, mo, j, mi)
                * real_y_at_phi_zero(lout, ml, x) * qw)
    return result


def channel_weights(jmax, x, quadrature_weights):
    """Dictionary {(J,S,Lout,Lin): angular_weights} for all usual NN waves.

    Includes both orientations of coupled J-1/J+1 waves. For J=0 the
    triplet channel is only 3P0. No isospin/Pauli selection is applied.
    """
    result = {}
    for j in range(jmax+1):
        keys = [(j, 0, j, j)]
        if j > 0:
            keys.append((j, 1, j, j))
        ls = [1] if j == 0 else [j-1, j+1]
        keys.extend((j, 1, lo, li) for lo in ls for li in ls)
        for key in keys:
            result[key] = angular_weights(*key, x, quadrature_weights)
    return result


def project(kernel, weights):
    """Contract angular/spin axes, retaining any leading batch dimensions."""
    return np.einsum("...xab,xab->...", kernel, weights, optimize=True)


def self_test():
    """Analytic central checks and tensor reciprocity independent of a potential."""
    x, qw = np.polynomial.legendre.leggauss(96)
    weights = channel_weights(8, x, qw)
    maximum_central_error = 0.0
    for (j, s, lo, li), weight in weights.items():
        eye = np.eye(2*s+1)
        for power in range(11):
            kernel = eval_legendre(power, x)[:, None, None] * eye
            expected = 4*pi/(2*li+1) if lo == li == power else 0.0
            maximum_central_error = max(maximum_central_error,
                                        abs(project(kernel, weight)-expected))
    assert maximum_central_error < 1e-11, maximum_central_error

    # Real tensor (sigma_1.q)(sigma_2.q), projected from product-spin
    # basis [upup,updown,downup,downdown] to triplet [-1,0,+1]. A smooth
    # q^2-dependent scalar makes off-diagonal angular integrals nontrivial.
    sx = np.array([[0., 1.], [1., 0.]])
    sz = np.diag([1., -1.])
    transform = np.array([[0., 0., 1.], [0., 1/sqrt(2), 0.],
                          [0., 1/sqrt(2), 0.], [1., 0., 0.]])

    def tensor_kernel(pout, pin):
        values = []
        for xx in x:
            qx = pout*sqrt(1-xx*xx)
            qz = pout*xx-pin
            sigq = qx*sx+qz*sz
            values.append(transform.T @ np.kron(sigq, sigq) @ transform
                          * np.exp(-0.13*(qx*qx+qz*qz)))
        return np.asarray(values)

    forward = tensor_kernel(0.7, 1.3)
    reverse = tensor_kernel(1.3, 0.7)
    maximum_reciprocity_error = 0.0
    for (j, s, lo, li), weight in weights.items():
        if s == 1:
            error = abs(project(forward, weight)
                        - project(reverse, weights[j, s, li, lo]))
            maximum_reciprocity_error = max(maximum_reciprocity_error, error)
    assert maximum_reciprocity_error < 1e-11, maximum_reciprocity_error
    # In raw angular states the S->D tensor integral has a POSITIVE sign.
    # The source2 inverse conversion of its named mixing contact has the
    # opposite sign; the caller's radial-wave phase converts between them.
    sd_tensor = project(tensor_kernel(1.0, 0.0), weights[1, 1, 2, 0])
    sd_expected = 4*pi*sqrt(8)/3*np.exp(-0.13)
    sd_tensor_error = abs(sd_tensor-sd_expected)
    assert sd_tensor_error < 1e-11, sd_tensor_error
    return {"maximum_central_error": float(maximum_central_error),
            "maximum_tensor_reciprocity_error": float(maximum_reciprocity_error),
            "sd_tensor_sign_and_normalization_error": float(sd_tensor_error)}


if __name__ == "__main__":
    import json
    print(json.dumps(self_test(), indent=2))
