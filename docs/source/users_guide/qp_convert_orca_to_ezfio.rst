.. _qp_convert_orca_to_ezfio:

qp_convert_orca_to_ezfio
========================

.. program:: qp_convert_orca_to_ezfio

This Python script uses the `trexio-tools`_ Python library to gather the
geometry, |AOs| and |MOs| from output files of `ORCA`_, and
puts this data in an |EZFIO| database. Some constraints are necessary
in the output file : the run needs to be a single point |HF|, |DFT| or
|CAS| |SCF|.

Usage
-----

.. code:: bash

    qp_convert_orca_to_ezfio [-o EZFIO_DIR] FILE

.. option:: -o, --output=EZFIO_DIR

    Renames the |EZFIO| directory. If this option is not present, the
    default name will be :file:`FILE.ezfio`

.. note::

   All the parameters of the wave function need to be present in the
   output file : complete description of the |AO| basis set, full set of
   molecular orbitals, etc.

   The following keywords are necessary for ORCA ::

      %output
        Print[P_Basis] 2
        Print[P_MOs] 1
      end


Example
-------

.. code:: bash

   qp_convert_orca_to_ezfio h2o.out -o h2o
