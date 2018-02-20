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

- Python 2.7.x
- `COMBO <https://github.com/tsudalab/combo>`_
- numpy
- pandas
- `pymatgen <http://pymatgen.org>`_


Structure optimizer
--------------------
- `VASP <https://www.vasp.at/>`_ or `Quantum ESPRESSO <http://www.quantum-espresso.org/>`_ or `soiap <https://github.com/nbsato/soiap>`_



.. index::
   single: find_wy

Others
--------

- `find_wy <https://github.com/nim-hrkn/find_wy>`_


Download
==========

Git
-----

.. code-block:: bash

    $ git clone https://github.com/Tomoki-YAMASHITA/CrySPY.git CrySPY_root

You can use an arbitarary directory name. Here just use ``CrySPY_root``.

Zip file
---------
You can download the source as a zip file from github: https://github.com/Tomoki-YAMASHITA/CrySPY


Setup
=========
You can put the source code in an arbitrary directory. Here, let's put the source code in ``~/CrySPY_root/``.

Directory tree in ``~/CrySPY_root/``::

    CrySPY_root
    ├── CrySPY
    │   ├── BO
    │   ├── IO
    │   ├── __init__.py
    │   ├── calc_dscrpt
    │   ├── f-fingerprint
    │   ├── find_wy
    │   ├── gen_struc
    │   ├── interface
    │   ├── job
    │   └── start
    ├── README.md
    ├── cryspy.py
    ├── docs
    ├── example
    ├── tutorial
    └── utility


find_wy
----------
Put the executable file of **find_wy** in ``~/CrySPY_root/CrySPY/find_wy/``, so that the executable file path is ``~/CrySPY_root/CrySPY/find_wy/find_wy``.

.. seealso::
   You need to install `find_wy <https://github.com/nim-hrkn/find_wy>`_ in advance.

.. note::
   Check ``/your_cryspy_path/CrySPY/find_wy/find_wy``



.. index::
   single: f-fingerprint

f-fingerprint
---------------
If you want to use Bayesian optimization, compile **cal_fingerpirnt** program.

.. code-block:: bash

    $ cd ~/CrySPY/CrySPY/f-fingerprint
    $ (edit Makefile)
    $ make

Make sure that the executable file of **cal_fingerprint** exist in ``~/CrySPY_root/CrySPY/f-fingerprint/``.

.. note::
   Check ``/your_cryspy_path/CrySPY/f-fingerprint/cal_fingerprint``
