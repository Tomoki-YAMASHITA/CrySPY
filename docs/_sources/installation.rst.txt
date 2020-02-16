=================
Installation
=================

.. contents:: Contents


System requirements
====================


.. index::
   single: COMBO

Python
--------

- Python 3.x.x
- `COMBO <https://github.com/tsudalab/combo3>`_
- numpy
- pandas
- `pymatgen <http://pymatgen.org>`_

Tested with Python 3.7.5 on Mac and Python 3.6.7 (miniconda3-4.3.30) on Linux (super computer).



Structure optimizer
--------------------
At least one optimizer is required.

- `VASP <https://www.vasp.at>`_
- `Quantum ESPRESSO <http://www.quantum-espresso.org>`_
- `soiap <https://github.com/nbsato/soiap>`_
- `LAMMPS <http://lammps.sandia.gov>`_



.. index::
   single: find_wy

Others
--------

- `find_wy <https://github.com/nim-hrkn/find_wy>`_



Installation of find_wy
========================

m_tspace
---------

First you need compile ``m_tspace`` for ``find_wy``. Download the source code of ``m_tspace`` in an arbitrary directory.
For example:

.. code-block:: bash

    $ mkdir -p ~/local/src
    $ cd ~/local/src/
    $ git clone https://github.com/nim-hrkn/m_tspace.git

.. seealso::
   | https://github.com/nim-hrkn/m_tspace
   | https://github.com/nim-hrkn/m_tspace/wiki


Additional two files are required for ``m_tspace``. Download the following files in ``~/local/src/m_tspace``:

- tsp98.f: http://phoenix.mp.es.osaka-u.ac.jp/~tspace/tspace_main/tsp07/tsp98.f
- prmtsp.f: http://phoenix.mp.es.osaka-u.ac.jp/~tspace/tspace_main/tsp07/prmtsp.f


.. code-block:: bash

    $ cd m_tspace
    $ wget http://phoenix.mp.es.osaka-u.ac.jp/~tspace/tspace_main/tsp07/tsp98.f
    $ wget http://phoenix.mp.es.osaka-u.ac.jp/~tspace/tspace_main/tsp07/prmtsp.f


Edit the ``makefile`` and run the ``make`` command. If you use ``ifort``, you had better delete ``-check all`` option and use ``-O2`` option.

.. code-block:: bash

    $ emacs makefile
    $ head -n 4 makefile
    #FC=gfortran
    #FFLAGS=-g -cpp -DUSE_GEN -ffixed-line-length-255
    FC=ifort
    FFLAGS=-O2 -g -traceback -cpp -DUSE_GEN -132
    $ make

If you used ``gfortran``, you might face the following problem:

.. code-block:: bash

    tsp98.f:9839:32:
    
           CALL SUBGRP(MG,JG,MGT,JGT,NTAB,IND)
                                    1
    Error: Actual argument contains too few elements for dummy argument 'ntab' (12/48) at (1)
    make: *** [tsp98.o] Error 1

Then change the source file of tsp98.f like this (line 9925):

Before:

.. code-block:: bash
    :emphasize-lines: 13

    9913: C SUBROUTINE SUBGRP ====*====3====*====4====*====5====*====6====*====7
    9914: C
    9915: C    IF (JG(I),I=1,MG) IS A SUBGROUP OF (JGT(J),J=1,MGT) THEN 
    9916: C          TABLE (NTAB(I),I=1,MG) IS MADE HERE AND IND=0
    9917: C    ELSE 
    9918: C          IND=-1
    9919: C
    9920: C                 1993/12/25
    9921: C                   BY  S.TANAKA AND A. YANASE
    9922: C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
    9923: C
    9924:       SUBROUTINE SUBGRP(MG,JG,MGT,JGT,NTAB,IND)
    9925:       DIMENSION NTAB(48),JG(48),JGT(48)


After:

.. code-block:: bash
    :emphasize-lines: 13

    9913: C SUBROUTINE SUBGRP ====*====3====*====4====*====5====*====6====*====7
    9914: C
    9915: C    IF (JG(I),I=1,MG) IS A SUBGROUP OF (JGT(J),J=1,MGT) THEN 
    9916: C          TABLE (NTAB(I),I=1,MG) IS MADE HERE AND IND=0
    9917: C    ELSE 
    9918: C          IND=-1
    9919: C
    9920: C                 1993/12/25
    9921: C                   BY  S.TANAKA AND A. YANASE
    9922: C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
    9923: C
    9924:       SUBROUTINE SUBGRP(MG,JG,MGT,JGT,NTAB,IND)
    9925:       DIMENSION NTAB(12),JG(48),JGT(48)


You are supposed to obtain ``m_tsp.a``.



find_wy
--------

Download the source code of ``find_wy`` in an arbitrary directory.
For example:

.. code-block:: bash

    $ mkdir -p ~/local/src
    $ cd ~/local/src/
    $ git clone https://github.com/nim-hrkn/find_wy.git

.. seealso::
   | https://github.com/nim-hrkn/find_wy
   | https://github.com/nim-hrkn/find_wy/wiki

Edit the ``make.inc`` and set the path to ``m_tsp.a`` that you just prepared.

.. code-block:: bash

    $ cd find_wy
    $ emacs make.inc
    $ head -n 4 make.inc
    TSPPATH=~/local/src/m_tspace
    #INCPATH = -I $(TSPPATH)
    TSP=$(TSPPATH)/m_tsp.a

You can delete ``-check all`` option and use ``-O2`` option. Then run the ``make`` command.

.. code-block:: bash

    $ make

When you got the executable file of ``find_wy``, run the following command for test:

.. code-block:: bash

    $ ./find_wy input_sample/input_si4o8.txt

If there is no problem, ``POS_WY_SKEL_ALL.json`` file is generated.


Download
==========

You can put the source code of CrySPY in an arbitrary directory.
Here, let us put the source code in ``~/CrySPY_root/CrySPY-x.x.x`` (x.x.x means the version).

Git
-----

.. code-block:: bash

    $ mkdir ~/CrySPY_root
    $ cd ~/CrySPY_root
    $ git clone https://github.com/Tomoki-YAMASHITA/CrySPY.git CrySPY-x.x.x


Zip or tar.gz file
-------------------

You can also download the source as a zip or tar.gz file from github: https://github.com/Tomoki-YAMASHITA/CrySPY/releases


Setup
=========
Directory tree in ``~/CrySPY_root/CrySPY-x.x.x/``::

    CrySPY-x.x.x
    ├── CHANGELOG.md
    ├── CrySPY/
    │   ├── BO/
    │   ├── EA/
    │   ├── IO/
    │   ├── LAQA/
    │   ├── RS/
    │   ├── __init__.py
    │   ├── calc_dscrpt/
    │   ├── f-fingerprint/
    │   ├── find_wy/
    │   ├── gen_struc/
    │   ├── interface/
    │   ├── job/
    │   └── start/
    │   └── utility.py
    ├── LICENSE
    ├── README.md
    ├── cryspy.py
    ├── docs/
    ├── example/
    └── utility/


find_wy
----------
Put the executable file of ``find_wy`` in ``~/CrySPY_root/CrySPY-x.x.x/CrySPY/find_wy/``, so that the executable file path is ``~/CrySPY_root/CrySPY-x.x.x/CrySPY/find_wy/find_wy``.

.. code-block:: bash

    $ cd ~/CrySPY_root/CrySPY-x.x.x/CrySPY/find_wy
    $ cp ~/local/src/find_wy/find_wy .


.. index::
   single: f-fingerprint

f-fingerprint
---------------
If you use Bayesian optimization, compile ``cal_fingerpirnt`` program.

.. code-block:: bash

    $ cd ~/CrySPY_root/CrySPY-x.x.x/CrySPY/f-fingerprint
    $ emacs Makefile
    $ make

Make sure that the executable file of ``cal_fingerprint`` exists in ``~/CrySPY_root/CrySPY-x.x.x/CrySPY/f-fingerprint/``.

.. note::
   Check ``/your_cryspy_path/CrySPY/f-fingerprint/cal_fingerprint``
