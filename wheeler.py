
import numpy as np
from scipy.linalg import eigh_tridiagonal

def wheeler(m):
    """
    David Lignell
    Wheeler algorithm for computing weights and abscissas from moments.
    From Marchisio and Fox (2013) Computational Models for Polydisperse and Multiphase systems.
    
    input m array of moments (size = 2N)
    returns w, x (weights and abscissas)
    """

    N2 = len(m)
    N  = int(N2/2)

    sigma  = np.zeros((N+1, N2))
    a      = np.zeros(N)
    b      = np.zeros(N)
    j_diag = np.zeros(N)
    j_ldiag = np.zeros(N)

    sigma[1,:N2] = m

    a[0] = m[1]/m[0]

    for k in range(1, N):
        l = np.arange(k,N2-k)
        sigma[k+1,l] = sigma[k,l+1]-a[k-1]*sigma[k,l]-b[k-1]*sigma[k-1,l]
        a[k] = -sigma[k,k]/sigma[k,k-1]+ sigma[k+1,k+1]/sigma[k+1,k]
        b[k] = sigma[k+1,k]/sigma[k,k-1]

    j_diag  = a
    j_ldiag = -np.sqrt(np.abs(b[1:]))

    x, v = eigh_tridiagonal(j_diag, j_ldiag)

    w = v[0,:]**2 * m[0]

    return w, x

