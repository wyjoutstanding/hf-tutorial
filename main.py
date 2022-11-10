import os, typing
import numpy
import scipy

# Use eigh to diagonalize matrices
from scipy.linalg import eigh
from scipy.linalg import fractional_matrix_power

import sol

def solve_rhf(nelecs, hcore: numpy.ndarray, ovlp: numpy.ndarray, eri: numpy.ndarray,
              ene_nuc :float = 0.0, max_iter :int = 100, tol: float = 1e-8) -> float:
    """
    Solve the Hartree-Fock with SCF iterations.
    Reference: Szabo and Ostlund, 3.4.6. (p. 146, start from step 2)

    The SCF procedure is:
        - Obtain a guess at the density matrix.
        - Calculate the exchange and coulomb matrices from the density matrix
          and the two-electron repulsion integrals.
        - Add exchange and coulomb matrices to the core-Hamiltonian to obtain the
          Fock matrix.
        - Diagonalize the Fock matrix.
        - Select the occupied orbitals and calculate the new density matrix.
        - Compute the energy
        - Compute the errors and check convergence
            - If converged, return the energy
            - If not converged, return to second step with new density matrix

    Inputs:
        nelecs : tuple
            the number of alpha and beta electrons
        hcore : numpy.ndarray
            The core Hamiltonian matrix.
        ovlp : numpy.ndarray
            The overlap matrix.
        eri : numpy.ndarray
            The two-electron repulsion integrals.
        ene_nuc : float
            The nuclear repulsion energy.
        max_iter : int
            The maximum number of SCF iterations.
        tol : float, optional
            The convergence tolerance, by default 1e-8.
    """

    nelec_alph, nelec_beta = nelecs
    assert nelec_alph == nelec_beta, "This code only supports closed-shell systems."
    n_electrons = nelec_alph+nelec_beta

    nao = hcore.shape[0]

    assert hcore.shape == (nao, nao)
    assert ovlp.shape  == (nao, nao)
    assert eri.shape   == (nao, nao, nao, nao)
    print("type: ", hcore.dtype, ovlp.dtype, eri.dtype, ene_nuc.dtype)
    print("shape: ", hcore.shape, ovlp.shape, eri.shape, ene_nuc.shape)

    print("Great! We are ready to solve the Hartree-Fock equations...")

    dm_init = None # Initialize the density matrix here.
    # dm_init = numpy.random.randn(nao, nao)

    iter_scf     = 0
    is_converged = False
    is_max_iter  = False

    ene_err = 1.0
    dm_err  = 1.0

    ene_rhf = None
    ene_old = None
    ene_cur = None

    assert dm_init is None
    dm_old    = dm_init

    # number of occupied orbitals
    n_occ = n_electrons // 2
    print("n_occ:", n_occ)

    dm_cur = numpy.zeros((nao, nao))
    F = numpy.zeros((nao, nao))
    
    numpy.random.seed(10213)
    C = numpy.random.randn(nao, n_occ)
    while not is_converged and not is_max_iter:
        # Fill in the code here to perform SCF iterations.
        # print("C1:", C.shape)
        for alph in range(nao):
            for beta in range(nao):
                _t = 0.0
                for i in range(n_occ):
                    _t += C[alph, i] * C[beta, i]
                dm_cur[alph, beta] = _t * 2.0

        # _tmpC = C @ C.T
        # assert numpy.allclose(_tmpC, dm_cur, 1e-8)
        # print("dm_cur", dm_cur)
        # print("dm_old", dm_old)

        for u in range(nao):
            for v in range(nao):
                _t = 0.0
                for alph in range(nao):
                    for beta in range(nao):
                        _t += dm_cur[alph,beta] * (eri[alph,beta,u,v] - eri[alph,v,u,beta]/2.0)
                F[u,v] = hcore[u,v] + _t

        # coul =   numpy.einsum("pqrs,rs->pq", eri, dm_cur)
        # exch = - numpy.einsum("prqs,rs->pq", eri, dm_cur) / 2.0
        # F = hcore + coul + exch

        # Diagonalize the Fock matrix
        
        # origin method
        # _ovlp = fractional_matrix_power(ovlp, -0.5)
        # _F = _ovlp @ F @ _ovlp
        # E, _C = eigh(_F)
        # C = (_ovlp @ _C)[:, 0:n_occ]
        
        # fusion method
        E, _C = eigh(F, ovlp)
        C = _C[:, 0:n_occ]

        ene_cur = 0.5 * numpy.einsum("pq,pq->", hcore + F, dm_cur)
        # ene_cur = numpy.einsum("pq,pq->", hcore + F, dm_cur)
        ene_rhf = ene_cur + ene_nuc

        print("ene_rhf:", ene_rhf, ene_old)
        # Compute the errors
        if ene_old is not None:
            dm_err  = numpy.linalg.norm(dm_cur - dm_old)
            ene_err = abs(ene_cur - ene_old)
            print(f"SCF iteration {iter_scf:3d}, energy = {ene_rhf: 12.8f}, error = {ene_err: 6.4e}, {dm_err: 6.4e}")

        # dm_old = dm_cur 
        dm_old = numpy.copy(dm_cur)
        ene_old = ene_cur

        iter_scf += 1
        is_max_iter  = iter_scf >= max_iter
        is_converged = ene_err < tol and dm_err < tol
        # is_converged = dm_err < tol

    # ene_rhf = ene_cur
    if is_converged:
        print(f"SCF converged in {iter_scf} iterations.")
    else:
        if ene_rhf is not None:
            print(f"SCF did not converge in {max_iter} iterations.")
        else:
            raise RuntimeError("SCF is not running, fill in the code in the main loop.")

    return ene_rhf

def solve_uhf(nelecs, hcore: numpy.ndarray, ovlp: numpy.ndarray, eri: numpy.ndarray,
              ene_nuc :float = 0.0, max_iter :int = 100, tol: float = 1e-8) -> float:
    raise NotImplementedError

def main(inp: str) -> None:
    """
    The main function of the program.

    Inputs:
        inp : str
            the input can be either h2o, which gives the integrals for water at equilibrium geometry,
            or some
    """

    inp = inp.split('-')
    mol = inp[0]

    tol = 1e-8

    if len(inp) == 2:
        int_dir = f"./data/{mol}"
        if not os.path.exists(int_dir):
            raise RuntimeError("Molecule not supported.")

        int_dir = f"./data/{mol}/{float(inp[1]):.4f}"
        if not os.path.exists(int_dir):
            raise RuntimeError("Bond length not supported.")

        nelecs  = numpy.load(f"{int_dir}/nelecs.npy")
        hcore   = numpy.load(f"{int_dir}/hcore.npy")
        ovlp    = numpy.load(f"{int_dir}/ovlp.npy")
        eri     = numpy.load(f"{int_dir}/eri.npy")
        ene_nuc = numpy.load(f"{int_dir}/ene_nuc.npy")

        ene_rhf_ref = numpy.load(f"{int_dir}/ene_rhf.npy")
        # ene_uhf_ref = numpy.load(f"{int_dir}/ene_uhf.npy")

        # Implement your restricted Hartree-Fock algorithm here.
        ene_rhf = solve_rhf(nelecs, hcore, ovlp, eri, tol=tol, max_iter=100, ene_nuc=ene_nuc)
        # This is the solution for RHF from Junjie, uncomment this to run it.
        # ene_rhf = sol.solve_rhf(nelecs, hcore, ovlp, eri, tol=tol, max_iter=200, ene_nuc=ene_nuc)

        # Implement your unrestricted Hartree-Fock algorithm here.
        # ene_uhf = solve_uhf(nelecs, hcore, ovlp, eri, tol=tol, max_iter=100, ene_nuc=ene_nuc)

        print(f"RHF energy: {ene_rhf: 12.8f}, Ref: {ene_rhf_ref: 12.8f}, Err: {abs(ene_rhf - ene_rhf_ref): 6.4e}")
        assert abs(ene_rhf - ene_rhf_ref) < tol

        return ene_rhf

    else:
        raise RuntimeError("Invalid input.")

if __name__ == "__main__":
    # mol can be either h2, heh+ or h2o.
    # r is the bond length in Angstrom.
    # mol = "heh+"
    # mol = "h2o"
    mol = "h2"
    r   = 1.0
    inp = f"{mol}-{r:.4f}"
    ene = main(inp)
