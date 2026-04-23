cap — Complex Absorbing Potential (CAP) module
===============================================

Overview
--------

This module implements the Complex Absorbing Potential (CAP) method for computing
resonance energies and widths of metastable electronic states in the framework of
Quantum Package (QP2).

The CAP Hamiltonian is:

.. math::

   H_{\text{CAP}} = H - i \eta W

where :math:`H` is the physical Hamiltonian, :math:`\eta` is the CAP strength
parameter, and :math:`W` is a box-type absorbing potential that is non-zero only
beyond the onset radii ``cap_onset_{x,y,z}``.

Physical background
-------------------

Adding a CAP to the Hamiltonian transforms the resonance energy into a complex
eigenvalue 
:math:`E_r - i\Gamma/2`, where :math:`E_r` is the resonance energy and
:math:`\Gamma` is the resonance width (related to the lifetime). The CAP strength 
:math:`\eta` is varied and the resonance is identified by the stationary point of
:math:`\eta \, dE/d\eta` (the :math:`\eta`-trajectory method).

Key providers and subroutines
------------------------------

``cap.irp.f``
    Main driver. Loops over ``n_steps_cap`` values of ``eta_cap`` and prints the
    complex energies for all states.

``diagonalize_ci.irp.f``
    Diagonalizes the CAP-CI Hamiltonian using either Davidson or LAPACK. Performs
    c-normalization (biorthogonal inner product ``(u|v) = sum_i u_i v_i``) and
    computes first-order corrections.

``dav_general_complex.irp.f``
    Implementation of the complex Davidson algorithm for iterative diagonalization
    of large sparse complex matrices.

``algebra.irp.f``
    Linear algebra utilities: complex QR decomposition (non-conjugate Gram-Schmidt),
    LAPACK wrappers (``zgeev``, ``zggev``), inner products, normalization.

``slater_rules.irp.f``
    CAP matrix element :math:`\langle i | W | j \rangle` via Slater-Condon rules.

``one_e_rdm.irp.f``
    One-body density matrix ``dm(p,q,istate,jstate)`` in the MO basis.

``density_matrix.irp.f``
    One-body density matrix in the MO basis for the CAP wavefunction.

``full_pt2.irp.f``
    Second-order perturbation correction to the CAP-CI energy.

``mo_ints.irp.f``
    CAP integrals in the MO basis: ``mo_wx_cap``, ``mo_wy_cap``, ``mo_wz_cap``.

``ao_ints.irp.f``
    CAP integrals in the AO basis.

``write_wf.irp.f``
    Write the CAP wavefunction to disk.

``overlap.irp.f``
    Biorthogonal overlap between CAP eigenstates.

``providers.irp.f``
    IRP providers for CAP-specific quantities (``psi_cap_coef``, ``psi_energy_cap``,
    ``eta_cap``, etc.).

EZFIO parameters
----------------

- ``eta_cap``          : CAP strength parameter :math:`\eta` (a.u.).
- ``n_steps_cap``      : Number of :math:`\eta` values to scan.
- ``eta_step_size``    : Step size between consecutive :math:`\eta` values.
- ``cap_onset_x/y/z``  : Onset radii of the box CAP along each axis (a.u.).
- ``do_cap``           : Whether to run the CAP calculation.

Scripts
-------

See ``scripts/README.rst`` for the post-processing scripts.
