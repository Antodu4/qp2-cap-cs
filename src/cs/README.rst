cs — Complex Scaling (CS) module
=================================

Overview
--------

This module implements the Complex Scaling (CS) method for computing resonance
energies and widths of metastable electronic states in the framework of
Quantum Package (QP2).

In complex scaling, the electron coordinates are rotated into the complex plane by
the angle :math:`\theta`:

.. math::

   \mathbf{r} \to e^{i\theta} \mathbf{r}

This transforms the resonance energy into a complex eigenvalue
:math:`E_r - i\Gamma/2`, where :math:`E_r` is the resonance energy and
:math:`\Gamma` is the resonance width (lifetime).

The CS Hamiltonian is:

.. math::

   H_{\text{CS}} = e^{-i\theta} \left( H + (e^{-i\theta} - 1) \, T \right)

where :math:`H` is the physical Hamiltonian and :math:`T` is the kinetic energy
operator. The resonance is identified as a stationary point of the complex energy
with respect to :math:`\theta`.

Key providers and subroutines
------------------------------

``cs.irp.f``
    Main driver. Loops over ``n_steps_cs`` values of ``theta_cs`` and prints the
    complex CS energies for all states.

``diagonalize_ci_cs.irp.f``
    Diagonalizes the CS-CI Hamiltonian using either Davidson or LAPACK. Performs
    c-normalization (biorthogonal inner product ``(u|v) = sum_i u_i v_i``),
    state selection by S², and computes first-order corrections via the one-body
    density matrix.

``H_cs_matrix.irp.f``
    Builds the full CS Hamiltonian matrix in the determinant basis:
    ``H_matrix_all_dets_cs(i,j)``.

``diagonalization_hs2_complex_cs.irp.f``
    Applies :math:`H_{\text{CS}}` and :math:`S^2` to trial vectors for the complex
    Davidson algorithm.

``slater_rules_cs.irp.f``
    CS kinetic matrix element :math:`\langle i | T | j \rangle` via Slater-Condon rules.

``one_e_rdm_cs.irp.f``
    One-body density matrix in the MO basis for the CS wavefunction.

``density_matrix_cs.irp.f``
    One-body density matrix for the CS wavefunction.

``full_pt2_cs.irp.f``
    Second-order perturbation correction to the CS-CI energy.

``overlap_cs.irp.f``
    Biorthogonal overlap between CS eigenstates.

``providers_cs.irp.f``
    IRP providers for CS-specific quantities (``psi_cs_coef``, ``psi_energy_cs``,
    ``theta_cs``, etc.).

``write_cs_wf.irp.f``
    Write the CS wavefunction to disk.

EZFIO parameters
----------------

- ``theta_cs``       : Complex scaling angle :math:`\theta` (radians).
- ``n_steps_cs``     : Number of :math:`\theta` values to scan.
- ``theta_step_size``: Step size between consecutive :math:`\theta` values.
- ``do_cs``          : Whether to run the CS calculation.
- ``cs_write_wf``    : Whether to write the CS wavefunction to disk after each step.

Scripts
-------

See ``scripts/README.rst`` for the post-processing scripts.
