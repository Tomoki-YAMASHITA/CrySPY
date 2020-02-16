=================
Utility
=================

[github] https://github.com/Tomoki-YAMASHITA/CrySPY/tree/master/utility

.. contents:: Contents


.. index::
   single: cryspy_analyzer_RS.ipynb

cryspy_analyzer_RS.ipynb
==========================

[github] https://github.com/Tomoki-YAMASHITA/CrySPY/blob/master/utility/cryspy_analyzer_RS.ipynb

Jupyter Notebook to visualize results for Random Search.


.. index::
   single: cryspy_analyzer_BO.ipynb

cryspy_analyzer_BO.ipynb
=================================

[github] https://github.com/Tomoki-YAMASHITA/CrySPY/blob/master/utility/cryspy_analyzer_BO.ipynb

Jupyter Notebook to visualize results for Bayesian Optimization


.. index::
   single: cryspy_analyzer_LAQA.ipynb

cryspy_analyzer_LAQA.ipynb
==========================

[github] https://github.com/Tomoki-YAMASHITA/CrySPY/blob/master/utility/cryspy_analyzer_LAQA.ipynb

Jupyter Notebook to visualize results for LAQA.


.. index::
   single: cryspy_analyzer_EA.ipynb

cryspy_analyzer_EA.ipynb
=================================

[github] https://github.com/Tomoki-YAMASHITA/CrySPY/blob/master/utility/cryspy_analyzer_EA.ipynb

Jupyter Notebook to visualize results for Evolutionary algorithm


.. index::
   single: pkl_data.ipynb

pkl_data.ipynb
================

[github] https://github.com/Tomoki-YAMASHITA/CrySPY/blob/master/utility/pkl_data.ipynb

Jupyter Notebook to analyze pickled data in CrySPY.



.. index::
   single: only_structure_generation

only_structure_generation
==============================================

``only_structure_generation`` directory includes:

- ``structure_generation.ipynb``
- ``sample_data_Si16_for_EA`` directory


structure_generation.ipynb
-----------------------------

[github] https://github.com/Tomoki-YAMASHITA/CrySPY/blob/master/utility/only_structure_generation/structure_generation.ipynb

Jupyter Notebook for tutorial to generate structures using CrySPY.
You can generate crystal structures with RS and EA partly using CrySPY  
(Here no CSP, only structure generation) 


sample_data_Si16_for_EA
---------------------------

Sample data used in ``structure_generation.ipynb`` for parents in EA.


.. index::
   single: kpt_check.py

kpt_check.py
===============
``kpt_check.py`` can check a k-point mesh with a given ``kppvol``. In this script, ``POSCAR``, ``CONTCAR``, and ``init_struc_data.pkl`` are readable.

For example, check a k-point mesh with kppvol = 100 for a POSCAR file.

.. code-block:: bash

   $ python kpt_check.py POSCAR 100
   a = 10.689217
   b = 10.689217
   c = 10.730846
       Lattice vector
   10.689217 0.000000 0.000000
   0.000000 10.689217 0.000000
   0.000000 0.000000 10.730846

   kppvol:  100
   k-points:  [2, 2, 2]

You can write a ``KPOINTS`` file with ``-w`` or ``--write`` option for VASP.

.. code-block:: bash

   $ python kpt_check.py -w POSCAR 100
   $ cat KPOINTS
   pymatgen 4.7.6+ generated KPOINTS with grid density = 607 / atom
   0
   Monkhorst
   2 2 2

In checking k-point meshes for init_struc_data.pkl, first five structures in init_struc_data.pkl are automatically checked in the default setting. You can change the number of structures using ``-n`` or ``--nstruc`` option.

.. code-block:: bash

   $ python kpt_check.py -n 3 init_struc_data.pkl 100


   # ---------- 0th structure
   a = 8.0343076893
   b = 8.03430768936
   c = 9.1723323373
       Lattice vector
   8.034308 0.000000 0.000000
   -4.017154 6.957915 0.000000
   0.000000 0.000000 9.172332

   kppvol:  100
   k-points:  [3, 3, 3]


   # ---------- 1th structure
   a = 9.8451944096
   b = 9.84519440959
   c = 6.8764313585
       Lattice vector
   9.845194 0.000000 0.000000
   -4.922597 8.526188 0.000000
   0.000000 0.000000 6.876431

   kppvol:  100
   k-points:  [3, 3, 4]


   # ---------- 2th structure
   a = 7.5760383679
   b = 7.57603836797
   c = 6.6507478296
       Lattice vector
   7.576038 0.000000 0.000000
   -3.788019 6.561042 0.000000
   0.000000 0.000000 6.650748

   kppvol:  100
   k-points:  [4, 4, 4]



.. index::
   single: spg_check.py

spg_check.py
=================
``spg_check.py`` can check space group information of a specified file. Structure.from_file() in pymatgen is used in this code. Supported formats include CIF, POSCAR/CONTCAR, ... etc. XXX.vasp file (POSCAR format) is also supported in this code.

.. seealso::
   `pymatgen <http://pymatgen.org/>`_

.. code-block:: bash

   $ python spg_check.py Al2O3.vasp
   (u'R-3c', 167)

You can change a tolerance value for checking the space group with ``-t`` or ``--tolerance`` option (default value is 0.1).

.. code-block:: bash

   $ python spg_check.py -t 0.001 Al2O3.vasp
   (u'R-3c', 167)




.. index::
   single: struc2cif.py

struc2cif.py
===================
``struc2cif.py`` can convert a structure file to a cif file using pymatgen. Structure.from_file() in pymatgen is used in this code. Supported formats include CIF, POSCAR/CONTCAR, ... etc. XXX.vasp file (POSCAR format) is also supported in this code. (input file name + '.cif') file is generated.

.. seealso::
   `pymatgen <http://pymatgen.org/>`_

.. code-block:: bash

   $ python struc2cif.py POSCAR

You can change a tolerance value for checking the space group with ``-t`` or ``--tolerance`` option (default value is 0.1).

.. code-block:: bash

   $ python struc2cif.py -t 0.001 POSCAR



.. index::
   single: qe2vasp_cif.py

qe2vasp_cif.py
===================
``qe2vasp_cif.py`` can generate structure data in VASP and cif formats from an input and output of QE.


You can obtain input structure data (``in_struc.vasp`` and ``in_struc.cif``) from an input of QE (pwscf.in)

.. code-block:: bash

   $ python qe2vasp_cif.py pwscf.in

You can obtain optimized structure data (``out_struc.vasp`` and ``out_struc.cif``) from an output of QE (pwscf.out). The input of QE (pwscf.in) is also required as a first argument.

.. code-block:: bash

   $ python qe2vasp_cif.py pwscf.in pwscf.out
