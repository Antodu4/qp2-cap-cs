cs/scripts — Post-processing scripts for CS calculations
=========================================================

Overview
--------

This directory contains Python scripts for parsing and analysing the output of
CS (Complex Scaling) calculations produced by the ``cs`` module.

Scripts
-------

``classes.py``
    Data model used by the other scripts.

    - ``Energy``    : Complex CS energy (Re and Im) plus first-order correction at
      one :math:`\theta` value for one state.
    - ``Energies``  : Collection of ``Energy`` objects over a range of
      :math:`\theta` values.
    - ``State``     : Data for one electronic state at one CIPSI iteration:
      ``n_det``, real CI energy ``E``, PT2 correction, and an ``Energies`` object.
    - ``Step``      : One CIPSI iteration; collects ``State`` objects for all states.
    - ``Cipsi``     : All CIPSI iterations for one calculation; a list of ``Step``
      objects.

``read_cs_output.py``
    Parses one or more CS calculation output files and collects the data into
    a ``pandas`` DataFrame.  Prints a CSV-style summary (one row per
    state/theta/N_det combination) suitable for further analysis or plotting.

    Usage::

        python read_cs_output.py output1.dat [output2.dat ...]

    Output columns: ``state``, ``theta``, ``Re(E_cs)``, ``Im(E_cs)``,
    ``Re(Corr_E_cs)``, ``Im(Corr_E_cs)``, ``N_det``, ``E``, ``PT2``.

``deriv.py``
    Computes the :math:`\theta`-derivative of a complex CS energy and prints
    :math:`\theta \, dE/d\theta` — useful for locating resonances and verifying
    the :math:`\theta`-stability of the result.

    Usage::

        python deriv.py energy_file.dat

    The input file must have columns: ``theta  Re(E)  Im(E)  ...``

Dependencies
------------

- Python 3
- NumPy
- Pandas
- Matplotlib
