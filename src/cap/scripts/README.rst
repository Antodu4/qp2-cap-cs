cap/scripts — Post-processing scripts for CAP calculations
===========================================================

Overview
--------

This directory contains Python scripts for parsing and analysing the output of
CAP (Complex Absorbing Potential) calculations produced by the ``cap`` module.

Scripts
-------

``classes.py``
    Data model used by the other scripts.

    - ``Energy``    : Complex energy (Re and Im) plus first-order correction at
      one :math:`\eta` value for one state.
    - ``Energies``  : Collection of ``Energy`` objects over a range of :math:`\eta`
      values.
    - ``State``     : Data for one electronic state at one CIPSI iteration:
      ``n_det``, real CI energy ``E``, PT2 correction, and an ``Energies`` object.
    - ``Step``      : One CIPSI iteration; collects ``State`` objects for all states.
    - ``Cipsi``     : All CIPSI iterations for one calculation; a list of ``Step``
      objects.

``read_output.py``
    Parses one or more CAP calculation output files and collects the data into
    a ``pandas`` DataFrame.  The format of each file is detected automatically.

    Usage::

        python read_output.py output1.dat [output2.dat ...]

    *New format* (``cap_fci`` output containing ``TRAJ_CAP`` records, printed by
    ``cap_print_traj`` for both the ``fixed_ndet`` and ``per_eta`` algorithms):

    - full trajectory, CSV with columns ``iter``, ``N_det``, ``eta``, ``state``,
      ``Re(E_cap)``, ``Im(E_cap)``, ``Re(E_cap+PT2)``, ``Im(E_cap+PT2)``,
      ``Re(PT2)``, ``Im(PT2)``, ``variance``;
    - converged values per :math:`\eta` (last iteration of each
      (state, :math:`\eta`)): ``eta``, then for each state ``Re(E)``, ``Im(E)``,
      ``Re(E+PT2)``, ``Im(E+PT2)``.

    In ``fixed_ndet`` mode the PT2 is only computed at the base :math:`\eta`
    (first value of the grid); for the other values PT2 = 0 and E+PT2 = E.
    The first-order corrections are not part of this format.

    *Legacy format* (per-iteration "CAP energies" / "Energy corrections" tables):
    output columns ``state``, ``eta``, ``Re(E_cap)``, ``Im(E_cap)``,
    ``Re(Corr_E_cap)``, ``Im(Corr_E_cap)``, ``N_det``, ``E``, ``PT2``.

``deriv.py``
    Computes the :math:`\eta`-derivative of a complex CAP energy and prints
    :math:`\eta \, dE/d\eta` — the key quantity for locating resonances via the
    :math:`\eta`-trajectory method.

    Usage::

        python deriv.py energy_file.dat

    The input file must have columns: ``eta  Re(E)  Im(E)  ...``

Dependencies
------------

- Python 3
- NumPy
- Pandas
- Matplotlib
