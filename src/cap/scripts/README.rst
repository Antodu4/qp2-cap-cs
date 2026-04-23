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
    a ``pandas`` DataFrame.  Prints a CSV-style summary (one row per
    state/eta/N_det combination) suitable for further analysis or plotting.

    Usage::

        python read_output.py output1.dat [output2.dat ...]

    Output columns: ``state``, ``eta``, ``Re(E_cap)``, ``Im(E_cap)``,
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
