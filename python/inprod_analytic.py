"""
    Inner products module
    ======================

    Inner products between the truncated set of basis functions for the \

    ocean and atmosphere streamfunction fields.


    .. note :: These are calculated using the analytical expressions from \
	De Cruz, L., Demaeyer, J. and Vannitsem, S.: A modular \
        arbitrary-order ocean-atmosphere model: MAOOAM v1.0, Geosci. \ 
        Model Dev. Discuss. \
	and from \
	Cehelsky, P., & Tung, K. K. : Theories of multiple equilibria and \
	weather regimes-A critical reexamination. Part II: Baroclinic two-layer \
	models. Journal of the atmospheric sciences, 44(21), 3282-3303, 1987.

    .. note :: The python code is available here :  \
        `inprod_analytic.py <../_modules/inprod_analytic.html>`_ .

    :Example:

    >>> from inprod_analytic import init_inprod
    >>> init_inprod()

    Global variables
    -------------------

    >>> awavenum = np.empty(natm, dtype=object)
    >>> owavenum = np.empty(noc, dtype=object)
    >>> atmos = atm_tensors(natm)
    >>> ocean = ocean_tensors(noc)

    Dependencies
    ------------------------

    it uses the modules :

    >>> import numpy as np
    >>> from scipy.sparse import csr_matrix

    >>> from params_maooam import nbatm
    >>> from params_maooam import nboc
    >>> from params_maooam import natm
    >>> from params_maooam import noc
    >>> from params_maooam import n
    >>> from params_maooam import oms
    >>> from params_maooam import ams
    >>> from params_maooam import pi

    Classes
    -------

    * atm_wavenum(typ,P,N,H,Nx,Ny)
    * ocean_wavenum(P,H,Nx,Ny)
    * atm_tensors(natm)
    * ocean_tensors(noc)
"""

import numpy as np
from scipy.sparse import csr_matrix

from params_maooam import nbatm
from params_maooam import nboc
from params_maooam import natm
from params_maooam import noc
from params_maooam import n
from params_maooam import oms
from params_maooam import ams
from params_maooam import pi

sqrtof2 = 1.4142135381698608398437500000


class atm_wavenum:
    """
    Class to define atmosphere wavenumbers.

    Attributes :

    * typ (char) = 'A','K' or 'L'.
    * M (int)
    * P (int)
    * H (int)
    * Nx (int)
    * Ny (int)
    """

    def __init__(self, typ, P, M, H, Nx, Ny):
        """creates the wavenumbers"""
        self.typ = typ
        self.P = P
        self.M = M
        self.H = H
        self.Nx = Nx
        self.Ny = Ny

    def __repr__(self):
            return "P = {}, M= {},H={}, Nx= {}, Ny={}".format(
                self.P, self.M, self.H, self.Nx, self.Ny)


class ocean_wavenum:
    """
    Class to define ocean wavenumbers

    Attributes :

    * P (int)
    * H (int)
    * Nx (int)
    * Ny (int)

    """

    def __init__(self, P, H, Nx, Ny):
        """creates the wavenumbers"""
        self.P = P
        self.H = H
        self.Nx = Nx
        self.Ny = Ny

    def __repr__(self):
            return "P = {}, H= {}, Nx= {}, Ny={}".format(
                self.P, self.H, self.Nx, self.Ny)


class atm_tensors:
    """
    Class which contains all the coefficients
    a,c,d,s,b,g needed for the tensor computation :

    Attributes :

    * :math:`a_{i, j}`
    * :math:`c_{i, j}`
    * :math:`d_{i, j}`
    * :math:`s_{i, j}`
    * :math:`b_{i, j, k}`
    * :math:`g_{i, j, k}`

    Return :

    * The object will be name `atmos`.
    """

    def __init__(self, natm):

        self.a = np.zeros((natm, natm), dtype=float)
        self.c = np.zeros((natm, natm), dtype=float)
        self.d = np.zeros((natm, noc), dtype=float)
        self.s = np.zeros((natm, noc), dtype=float)
        self.b = np.zeros((natm, natm, natm))
        self.g = np.zeros((natm, natm, natm))

    def __repr__(self):
            return "matrice d : \n{} \n".format(self.d)

    # !--------------------------------------------------------!
    # ! 1. Inner products in the equations for the atmosphere  !
    # !--------------------------------------------------------!

    def calculate_a(self):
        r"""
        .. math::
	    a_{i, j} = (F_i, {\nabla}^2 F_j)

	.. note:: Eigenvalues of the Laplacian (atmospheric)
        """
        if natm == 0:
            exit("*** Problem with calculate_a : natm==0!***")
        for i in range(0, natm):
            ti = awavenum[i]
            self.a[i, i] = - (n**2) * ti.Nx**2 - ti.Ny**2

    def calculate_b(self):
        r"""
        .. math::
            b_{i, j, k} = (F_i, J(F_j, \nabla^2 F_k))

        .. note:: Atmospheric g and a tensors must be computed before \
            calling this routine.
        """

        if natm == 0:
            exit("*** Problem with calculate_b : natm==0!***")

        for i in range(0, natm):
            for j in range(0, natm):
                for k in range(0, natm):
                    val = self.a[k, k]*self.g[i, j, k]
                    self.b[i, j, k] = val

    def calculate_c_atm(self):
        """
        .. math::
            c_{i,j} = (F_i, \partial_x F_j)
        .. note:: Beta term for the atmosphere
            Strict function !! Only accepts KL type.
            For any other combination, it will not calculate anything.
        """

        if natm == 0:
            exit("*** Problem with calculate_c_atm : natm==0!***")

        for i in range(0, natm):
            for j in range(0, natm):
                val = 0.
                Ti = awavenum[i]
                Tj = awavenum[j]

                if ((Ti.typ, Tj.typ) == ('K', 'L')):
                    val = delta(Ti.M - Tj.H) * delta(Ti.P - Tj.P)
                    val = n * Ti.M * val

                if (val != 0.):
                    self.c[i, j] = val
                    self.c[j, i] = -val

    def calculate_d(self, ocean):
        r"""
        .. math::
            d_{i,j} = (F_i, \nabla^2 \eta_j)
        .. note:: Forcing of the ocean on the atmosphere.
            Atmospheric s tensor and oceanic M tensor must be computed
            before calling this routine !
        """

        if natm == 0:
            exit("*** Problem with calculate_d : natm==0!***")

        for i in range(0, natm):
            for j in range(0, noc):
                self.d[i, j] = self.s[i, j] * ocean.M[j, j]

    def calculate_g(self):
        """
        .. math::
            g_{i,j,k} = (F_i, J(F_j, F_k))
        .. note:: This is a strict function: it only accepts AKL, KKL
            and LLL types.
            For any other combination, it will not calculate anything.
        """

        if natm == 0:
            exit("*** Problem with calculate_h : natm==0!***")

        for i in range(0, natm):
            for j in range(0, natm):
                for k in range(0, natm):
                    Ti = awavenum[i]
                    Tj = awavenum[j]
                    Tk = awavenum[k]
                    val = 0.

                    if (Ti.typ, Tj.typ, Tk.typ) == ('A', 'K', 'L'):
                        vb1 = B1(Ti.P, Tj.P, Tk.P)
                        vb2 = B2(Ti.P, Tj.P, Tk.P)
                        val = -2 * (sqrtof2 / pi) * Tj.M * delta(Tj.M - Tk.H) \
                            * flambda(Ti.P + Tj.P + Tk.P)
                        if (val != 0):
                            val = val * (((vb1**2) / (vb1**2 - 1)) - ((vb2**2)
                            / (vb2**2 - 1)))

                    if ((Ti.typ, Tj.typ, Tk.typ) == ('K', 'K', 'L')):
                        vs1 = S1(Tj.P, Tk.P, Tj.M, Tk.H)
                        vs2 = S2(Tj.P, Tk.P, Tj.M, Tk.H)
                        val = vs1 * (delta(Ti.M - Tk.H - Tj.M) \
                            * delta(Ti.P - Tk.P + Tj.P) \
                            - delta(Ti.M - Tk.H - Tj.M) \
                            * delta(Ti.P + Tk.P - Tj.P) \
                            + (delta(Tk.H - Tj.M + Ti.M) \
                            + delta(Tk.H - Tj.M - Ti.M)) \
                            * delta(Tk.P + Tj.P - Ti.P)) \
                            + vs2 * (delta(Ti.M - Tk.H - Tj.M) \
                            * delta(Ti.P - Tk.P - Tj.P) \
                            + (delta(Tk.H - Tj.M - Ti.M) \
                            + delta(Ti.M + Tk.H - Tj.M)) \
                            * (delta(Ti.P - Tk.P + Tj.P) \
                            - delta(Tk.P - Tj.P + Ti.P)))

                    val = val*n

                    if (val != 0.):
                        self.g[i, j, k] = val
                        self.g[j, k, i] = val
                        self.g[k, i, j] = val
                        self.g[i, k, j] = -val
                        self.g[j, i, k] = -val
                        self.g[k, j, i] = -val

        for i in range(0, natm):
            for j in range(i+1, natm):
                for k in range(j+1, natm):

                    Ti = awavenum[i]
                    Ti = awavenum[i]
                    Tj = awavenum[j]
                    Tk = awavenum[k]

                    val = 0.

                    if ((Ti.typ, Tj.typ, Tk.typ) == ('L', 'L', 'L')):
                        vs3 = S3(Tj.P, Tk.P, Tj.H, Tk.H)
                        vs4 = S4(Tj.P, Tk.P, Tj.H, Tk.H)
                        val = vs3 * ((delta(Tk.H - Tj.H - Ti.H) \
                        - delta(Tk.H - Tj.H + Ti.H)) \
                        * delta(Tk.P + Tj.P - Ti.P) \
                        + delta(Tk.H + Tj.H - Ti.H) \
                        * (delta(Tk.P - Tj.P + Ti.P) \
                        - delta(Tk.P - Tj.P - Ti.P))) + vs4 \
                        * ((delta(Tk.H + Tj.H - Ti.H) \
                        * delta(Tk.P - Tj.P - Ti.P)) \
                        + (delta(Tk.H - Tj.H + Ti.H)\
                        - delta(Tk.H - Tj.H - Ti.H)) \
                        * (delta(Tk.P - Tj.P - Ti.P) \
                        - delta(Tk.P - Tj.P + Ti.P)))

                    val = val*n

                    if (val != 0.):
                        atmos.g[i, j, k] = val
                        atmos.g[j, k, i] = val
                        atmos.g[k, i, j] = val
                        atmos.g[i, k, j] = -val
                        atmos.g[j, i, k] = -val
                        atmos.g[k, j, i] = -val

    def calculate_s(self):
        """
        .. math::
            s_{i,j} = (F_i, \eta_j)
        .. note:: Forcing (thermal) of the ocean on the atmosphere.
        """

        if natm == 0:
            exit("*** Problem with calculate_s : natm==0!***")

        if noc == 0:
            exit("*** Problem with calculate_s : noc==0!***")

        for i in range(0, natm):
            for j in range(0, noc):

                Ti = awavenum[i]
                Dj = owavenum[j]

                val = 0.

                if (Ti.typ == 'A'):
                    val = flambda(Dj.H) * flambda(Dj.P + Ti.P)
                    if (val != 0.):
                        val = val * 8 * sqrtof2 * Dj.P / \
                        (pi**2 * (Dj.P**2 - Ti.P**2) * Dj.H)

                if (Ti.typ == 'K'):
                    val = flambda(2 * Ti.M + Dj.H) * delta(Dj.P - Ti.P)

                    if (val != 0):
                        val = val*4*Dj.H/(pi * (-4 * Ti.M**2 + Dj.H**2))

                if (Ti.typ == 'L'):
                    val = delta(Dj.P - Ti.P) * delta(2 * Ti.H - Dj.H)

                if (val != 0.):
                    self.s[i, j] = val


class ocean_tensors:
    """
    Class which contains all the coefficients k,m,n,w,o,c \
    needed for the tensor computation :

    Attributes :

    * :math:`K_{i,j}`
    * :math:`M_{i,j}`
    * :math:`N_{i,j}`
    * :math:`W_{i,j}`
    * :math:`O_{i,j,k}`
    * :math:`C_{i,j,k}`

    Return :

    * The object will be name ocean
    """

    def __init__(self, noc):

        self.K = np.zeros((noc, natm), dtype=float)
        self.M = np.zeros((noc, noc), dtype=float)
        self.N = np.zeros((noc, noc), dtype=float)
        self.W = np.zeros((noc, natm), dtype=float)

        self.O = np.zeros((noc, noc, noc), dtype=float)
        self.C = np.zeros((noc, noc, noc), dtype=float)

    def __repr__(self):
            return "matrice C : \n{}".format(self.C)

    # !--------------------------------------------------------!
    # ! 2. Inner products in the equations for the ocean       !
    # !--------------------------------------------------------!

    def calculate_K(self, atmos):
        r"""
        Forcing of the atmosphere on the ocean.

        .. math::
            K_{i,j} = (\eta_i, \nabla^2 F_j)

        .. note::
            atmospheric a and s tensors must be computed before calling
            this function !
        """

        if noc == 0:
            exit("*** Problem with calculate_K : noc==0!***")

        for i in range(0, noc):
            for j in range(0, natm):
                self.K[i, j] = atmos.s[j, i] * atmos.a[j, j]

    def calculate_M(self):
        r"""
        Forcing of the ocean fields on the ocean.

        .. math::
            M_{i,j} = (\eta_i, \nabla^2 \eta_j)

        """

        if noc == 0:
            exit("*** Problem with calculate_M : noc==0!***")

        for i in range(0, noc):
            Di = owavenum[i]
            self.M[i, i] = - (n**2) * Di.Nx**2 - Di.Ny**2

    def calculate_N(self):
        """
        Beta term for the ocean

        .. math::
            N_{i,j} = (\eta_i, \partial_x \eta_j)

        """

        if noc == 0:
            exit("*** Problem with calculate_N : noc==0!***")

        val = 0.

        for i in range(0, noc):
            for j in range(0, noc):
                Di = owavenum[i]
                Dj = owavenum[j]
                val = delta(Di.P - Dj.P) * flambda(Di.H + Dj.H)

                if (val != 0.):
                    self.N[i, j] = val * (-2) * Dj.H * Di.H * n / \
                    ((Dj.H**2 - Di.H**2) * pi)

    def calculate_O(self):
        """
        Temperature advection term (passive scalar)

        .. math::
            O_{i,j,k} = (\eta_i, J(\eta_j, \eta_k))
        """

        if noc == 0:
            exit("*** Problem with calculate_O : noc==0!***")

        val = 0.

        for i in range(0, noc):
            for j in range(i, noc):
                for k in range(i, noc):

                    Di = owavenum[i]
                    Dj = owavenum[j]
                    Dk = owavenum[k]

                    vs3 = S3(Dj.P, Dk.P, Dj.H, Dk.H)

                    vs4 = S4(Dj.P, Dk.P, Dj.H, Dk.H)

                    val = vs3*((delta(Dk.H - Dj.H - Di.H) \
                    - delta(Dk.H - Dj.H + Di.H)) \
                    * delta(Dk.P + Dj.P - Di.P) \
                    + delta(Dk.H + Dj.H - Di.H) \
                    * (delta(Dk.P - Dj.P + Di.P) \
                    - delta(Dk.P - Dj.P - Di.P))) \
                    + vs4 * ((delta(Dk.H + Dj.H - Di.H) \
                    * delta(Dk.P - Dj.P - Di.P)) \
                    + (delta(Dk.H - Dj.H + Di.H) \
                    - delta(Dk.H - Dj.H - Di.H)) \
                    * (delta(Dk.P - Dj.P - Di.P) \
                    - delta(Dk.P - Dj.P + Di.P))) \

                    val = val * n / 2

                    if (val != 0.):
                        self.O[i, j, k] = val
                        self.O[j, k, i] = val
                        self.O[k, i, j] = val
                        self.O[i, k, j] = -val
                        self.O[j, i, k] = -val
                        self.O[k, j, i] = -val

    def calculate_C_oc(self):
        r"""
        .. math::
            C_{i,j,k} = (\eta_i, J(\eta_j,\nabla^2 \eta_k))

        .. note :: Requires :math:`O_{i,j,k}` \

        and :math:`M_{i,j}` to be calculated beforehand.

        """

        if noc == 0:
            exit("*** Problem with calculate_C : noc==0!***")

        val = 0.

        for i in range(0, noc):
            for j in range(0, noc):
                for k in range(0, noc):

                    val = self.M[k, k] * self.O[i, j, k]

                    if (val != 0.):
                        self.C[i, j, k] = val

    def calculate_W(self, atmos):
        """
        Short-wave radiative forcing of the ocean.

        .. math::
            W_{i,j} = (\eta_i, F_j)

        .. note ::
            atmospheric s tensor must be computed before calling
            this function !
        """

        if noc == 0:
            exit("*** Problem with calculate_W : noc==0!***")

        if natm == 0:
            exit("*** Problem with calculate_W : natm==0!***")

        for i in range(0, noc):
            for j in range(0, natm):
                self.W[i, j] = atmos.s[j, i]


#  !-----------------------------------------------------!
#  !                                                     !
#  ! Definition of the Helper functions from Cehelsky    !
#  ! \ Tung                                              !
#  !                                                     !
#  !-----------------------------------------------------!


def B1(Pi, Pj, Pk):
    return (Pk+Pj)/float(Pi)


def B2(Pi, Pj, Pk):
    return (Pk-Pj)/float(Pi)


def delta(r):

    if (r == 0):
        return 1.

    else:
        return 0.


def flambda(r):

    if (r % 2 == 0):
        return 0.

    else:
        return 1.


def S1(Pj, Pk, Mj, Hk):
    return -((Pk * Mj + Pj * Hk)) / 2.


def S2(Pj, Pk, Mj, Hk):
    return (Pk * Mj - Pj * Hk) / 2.


def S3(Pj, Pk, Hj, Hk):
    return (Pk * Hj + Pj * Hk) / 2.


def S4(Pj, Pk, Hj, Hk):
    return (Pk * Hj - Pj * Hk) / 2.


# initialization of the variables
awavenum = np.empty(natm, dtype=object)
owavenum = np.empty(noc, dtype=object)
atmos = atm_tensors(natm)
ocean = ocean_tensors(noc)


def init_inprod():
    """creates and computes the inner products."""

    j = -1

    # awavenum definition
    # here is the constructor : def __init__(self,typ,P,M,H,Nx,Ny):

    for i in range(0, nbatm):

        if ams[i, 0] == 1:

            awavenum[j+1] = atm_wavenum(
                'A', ams[i, 1], 0, 0, 0, ams[i, 1])

            awavenum[j+2] = atm_wavenum(
                'K', ams[i, 1], ams[i, 0], 0, ams[i, 0], ams[i, 1])

            awavenum[j+3] = atm_wavenum(
                'L', ams[i, 1], 0, ams[i, 0], ams[i, 0], ams[i, 1])

            j = j+3

        else:

            awavenum[j+1] = atm_wavenum(
                'K', ams[i, 1], ams[i, 0], 0, ams[i, 0], ams[i, 1])

            awavenum[j+2] = atm_wavenum(
                'L', ams[i, 1], 0, ams[i, 0], ams[i, 0], ams[i, 1])

            j = j+2

    # owavenum definition
    # here is the constructor : def __init__(self,P,H,Nx,Ny):

    for i in range(0, nboc):
            owavenum[i] = ocean_wavenum(
                oms[i, 1], oms[i, 0], oms[i, 0]/2., oms[i, 1])

    # Computation of the atmospheric inner products tensors

    atmos.calculate_a()
    atmos.calculate_g()
    atmos.calculate_s()
    atmos.calculate_b()
    atmos.calculate_c_atm()

    # Computation of the oceanic inner products tensors

    ocean.calculate_M()
    ocean.calculate_N()
    ocean.calculate_O()
    ocean.calculate_C_oc()
    ocean.calculate_W(atmos)
    ocean.calculate_K(atmos)

    # A last atmospheric one that needs ocean.M

    atmos.calculate_d(ocean)

