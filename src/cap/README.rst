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
    Main driver. ``cap_scan(l_energy, do_print)`` loops over ``n_steps_cap``
    values of ``eta_cap``, returns the complex energies (E + V_nn) for all states
    and optionally prints them. ``cap()`` is a wrapper that prints the tables.

``cap_traj_out.irp.f``
    ``cap_print_traj``: unified trajectory output shared by both CAP-CIPSI
    algorithms (see below).

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
- ``cap_cipsi_per_eta``: Selects the CAP-CIPSI algorithm in ``cap_fci``
  (see below). Default: ``false``.

CAP-CIPSI algorithms (``cap_fci``)
----------------------------------

Two algorithms are available to follow the :math:`\eta`-trajectory with a
selected CI, mirroring the two Complex Scaling (``cs``) algorithms:

``fixed_ndet`` (``cap_cipsi_per_eta = false``, default)
    ``run_stochastic_cap_cipsi`` (``cipsi/stochastic_cap_cipsi.irp.f``).
    Determinants are selected with the CAP PT2 at the base :math:`\eta`
    (``eta_cap``). At each CIPSI iteration, the other :math:`\eta` values of the
    grid are only re-diagonalized on the same determinant set. The PT2 is
    therefore only available at the base :math:`\eta`.

``per_eta`` (``cap_cipsi_per_eta = true``)
    ``run_stochastic_cap_cipsi_per_eta`` (``cipsi/stochastic_cap_cipsi_per_eta.irp.f``).
    Outer loop over the :math:`\eta` grid: for each :math:`\eta`, a complete
    CIPSI is converged, starting each time from the same initial determinant
    set (cold start). The determinant space, the energy and the PT2 are thus
    specific to each :math:`\eta`.

Both algorithms end with the same summary printed by ``cap_print_traj``:

- a banner ``CAP SCAN SUMMARY  algorithm = fixed_ndet | per_eta``;
- block A: convergence vs N_det, one row per (state, :math:`\eta`, iteration);
- block B: final CAP energies per :math:`\eta`;
- block C: machine-readable records, one line per (state, :math:`\eta`,
  iteration), starting with ``TRAJ_CAP``:
  ``iter n_det eta state ReE ImE ReEPT2 ImEPT2 RePT2 ImPT2 variance``.

In ``fixed_ndet`` mode the per-iteration CAP energy tables are no longer
printed (block C replaces them); ``diagonalize_h_cap`` still prints them.

Example::

    qp set cap cap_onset_x 5.0
    qp set cap cap_onset_y 5.0
    qp set cap cap_onset_z 5.0
    qp set cap eta_cap 0.0
    qp set cap eta_step_size 1.e-4
    qp set cap n_steps_cap 20
    qp set cap cap_cipsi_per_eta true     # or false for fixed_ndet
    qp run cap_fci > cap_fci.out

Scripts
-------

See ``scripts/README.rst`` for the post-processing scripts.
