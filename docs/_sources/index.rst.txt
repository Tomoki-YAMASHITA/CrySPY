.. image:: ./figs/logo/cryspy_fix-03.png
  :width: 20cm
  :align: center

|

===========
CrySPY
===========

| CrySPY is a crystal structure prediction tool written in Python.
| You can automatically do the following things with CrySPY.

* Structure generation
* Submitting jobs for structure optimization
* Collecting data for structure optimization
* Selecting candidates using machine learning

You can also use CrySPY partly for the purpose of only struture generation and descriptor calculation for crystal structures.

| Source code: [github] https://github.com/Tomoki-YAMASHITA/CrySPY
| Download: [github] https://github.com/Tomoki-YAMASHITA/CrySPY/releases

News
===========

* [2020 February 16] version 0.8.0 released

    - Migrate to Python 3
    - CrySPY logo created
    - Change several variable names and data formats
    - Change style of output for energy: eV/cell --> eV/atom
    - No. of working directories corresponds to structure ID
    - New features

        + recalculation
        + manual select in BO

* [2018 December 5] version 0.7.0 released

    - New features

        + Evolutionary algorithm

* [2018 August 20] version 0.6.4 released
* [2018 July 2] version 0.6.3 released
* [2018 June 26] Version 0.6.2 released
* [2018 March 1] Version 0.6.1 released
* [2018 January 9]

    - Our paper on CrySPY has been published in Physical Review Materials, https://link.aps.org/doi/10.1103/PhysRevMaterials.2.013803



Reference
===========

Bayesian optimization
----------------------

- T. Yamashita, N. Sato, H. Kino, T. Miyake, K. Tsuda, and T. Oguchi, Phys. Rev. Mater. **2**, 013803 (2018).

    +  https://link.aps.org/doi/10.1103/PhysRevMaterials.2.013803


LAQA
-----

- K.Terayama, T. Yamashita, T. Oguchi, and K. Tsuda, npj Comput. Mater. **4**, 32 (2018).

    +  https://www.nature.com/articles/s41524-018-0090-y



Contents
============

.. toctree::
   :maxdepth: 1

   about
   installation
   tutorial
   input_file
   data_format
   utility
   tips
   faq




Indices
==================

* :ref:`genindex`

.. * :ref:`modindex`
.. * :ref:`search`
