.. index::
   single: cryspy.in

===========================
Input file: ``cryspy.in``
===========================

.. contents:: Contents




File format
=================

CrySPY uses the `configparser <https://docs.python.org/3/library/configparser.html>`_ module to read input file, ``cryspy.in`` .
``cryspy.in``  consists of sections, led by a ``[section]`` header and followed by ``name = value`` or ``name : value`` entries.
Section names and values are case sensitive, but names are not.
Lines beginning with ``#`` or ``;`` are ignored and may be used to provide comments.
Accepted bool values are ``1``, ``yes``, ``true``, and ``on``, which cause this method to return ``True``, and ``0``, ``no``, ``false``, and ``off``, which cause it to return False. These string values for bool are checked in a case-insensitive manner.
Some values are given in a space-separated manner.

.. seealso:: `configparser <https://docs.python.org/3/library/configparser.html>`_
.. attention::
   | section name: case sensitive
   | name: case insensitive
   | value: case sensitive except for bool



Example
=================

Random Search (RS)
----------------------------

``cryspy.in`` in CrySPY-x.x.x/example/soiap_bash_RS_Si16/
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

soiap, bash, Si\ `16`:sub:::

   [basic]
   algo = RS
   calc_code = soiap
   tot_struc = 6
   natot = 16
   atype = Si
   nat =  16
   nstage = 1
   njob = 1
   jobcmd = bash
   jobfile = job_cryspy

   [lattice]
   minlen = 5
   maxlen = 10
   dangle = 20
   mindist_1 = 1.8

   [soiap]
   soiap_infile = soiap.in
   soiap_outfile = soiap.out
   soiap_cif = initial.cif


``cryspy.in`` in CrySPY-x.x.x/example/VASP_qsub_RS_Na8Cl8/
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

VASP, qsub, Na\ `8`:sub: Cl\ `8`:sub:::

   [basic]
   algo = RS
   calc_code = VASP
   tot_struc = 10
   natot = 16
   atype = Na Cl
   nat = 8 8
   nstage = 2
   njob = 10
   jobcmd = qsub
   jobfile = job_cryspy

   [lattice]
   minlen = 5
   maxlen = 12
   dangle = 20
   mindist_1 = 2.2  1.8
   mindist_2 = 1.8  2.2

   [VASP]
   kppvol = 40 80


``cryspy.in`` in CrySPY-x.x.x/example/QE_qsub_RS_Si16/
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

QE, qsub, Si\ `16`:sub:::

   [basic]
   algo = RS
   calc_code = QE
   tot_struc = 5
   natot = 16
   atype = Si
   nat = 16
   nstage = 2
   njob = 5
   jobcmd = qsub
   jobfile = job_cryspy

   [lattice]
   minlen = 5
   maxlen = 10
   dangle = 20
   mindist_1 = 1.8

   [QE]
   qe_infile = pwscf.in
   qe_outfile = pwscf.out
   kppvol = 40 80



Bayesian Optimization (BO)
-----------------------------------

``cryspy.in`` in CrySPY-x.x.x/example/VASP_qsub_BO_Na8Cl8/
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

VASP, qsub, Na\ `8`:sub: Cl\ `8`:sub:::

   [basic]
   algo = BO
   calc_code = VASP
   tot_struc = 10
   natot = 16
   atype = Na Cl
   nat = 8 8
   nstage = 2
   njob = 2
   jobcmd = qsub
   jobfile = job_cryspy

   [lattice]
   minlen = 5
   maxlen = 12
   dangle = 20
   mindist_1 = 2.2  1.8
   mindist_2 = 1.8  2.2

   [VASP]
   kppvol = 40 80

   [BO]
   nselect_bo = 2
   dscrpt = FP
   score = TS
   fp_rmin = 0.5
   fp_rmax = 5.0
   fp_npoints = 10
   fp_sigma = 1.0


Look Ahead based on Quadratic Approximation (LAQA)
-----------------------------------------------------------------------




Evolutionary Algorithm (EA)
-------------------------------------

``cryspy.in`` in CrySPY-x.x.x/example/VASP_qsub_EA_Na8Cl8/
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

VASP, qsub, Na\ `8`:sub: Cl\ `8`:sub:::

   [basic]
   algo = EA
   calc_code = VASP
   tot_struc = 10
   natot = 16
   atype = Na Cl
   nat = 8 8
   nstage = 2
   njob = 10
   jobcmd = qsub
   jobfile = job_cryspy

   [lattice]
   minlen = 5
   maxlen = 12
   dangle = 20
   mindist_1 = 2.2  1.8
   mindist_2 = 1.8  2.2

   [VASP]
   kppvol = 40 80

   [EA]
   n_pop = 10
   n_crsov = 5
   n_perm = 2
   n_strain = 2
   n_rand = 1
   n_elite = 0
   n_fittest = 5
   slct_func = TNM
   t_size = 2
   maxgen_ea = 2


.. index::
   single: [basic]

[basic] section
==================

.. csv-table::
   :header: Name, Value, Default value, Description
   :widths: auto

   ``algo``, "``RS`` , ``BO``, ``LAQA``, ``EA``",  ,  Algorithm
   ``calc_code``, "``VASP``, ``QE``, ``soiap``, ``LAMMPS``",  , Caluculation code for structure optimization
   ``tot_struc``, int,  , Total number of structures
   ``natot``, int,  , Total number of atoms in a unit cell
   ``atype``, "atomic symbol [atomic symbol ...]",  , Atom type
   ``nat``, "int [int ...]",  , "Number of atoms in atom type1 [type2 ...]"
   ``nstage``, int,  , Number of calculation stages
   ``njob``, int,  , Number simultaneously submitted jobs
   ``jobcmd``, str ,  , "Specify a command to submit jobs, such as qsub"
   ``jobfile``, str,  , "Specify a jobfile to submit jobs for VASP, QE, and so on"


.. index::
   single: algo

``algo``
----------

Available algorithms for crystal structure prediction are:

- ``RS``: **R**\ andom **S**\ earch
- ``BO``: **B**\ ayesian **O**\ ptimization
- ``LAQA``: **L**\ ook **A**\ head based on **Q**\ uadratic **A**\ pproximation
- ``EA``: **E**\ volutionary **A**\ lgorithm

In using LAQA, automatically ``fs_step_flag`` = ``True`` in [option] section.





.. index::
   single: calc_code

``calc_code``
---------------

CrySPY is interfaced with:

- ``VASP``: **VASP** (https://www.vasp.at)
- ``QE``: **Q**\ uantum **E**\ spresso (http://www.quantum-espresso.org)
- ``soiap``: **soiap** (https://github.com/nbsato/soiap)
- ``LAMMPS``: **LAMMPS** (http://lammps.sandia.gov)



.. index::
   single: [lattice]
   single: minlen
   single: maxlen

[lattice] section
==================

.. csv-table::
   :header: Name, Value, Default value, Description
   :widths: auto

   ``minlen``, float,  ,  Minimum length of lattce vector [Å]
   ``maxlen``, float,  ,  Maximum length of lattce vector [Å]
   ``dangle``, float,  ,  "Delta angle for alpha, beta, and gamma in degree unit"
   ``mindist_?``, float [float ...], ,  Constraint on minimum interatomic distance [Å]


.. index::
   single: dangle

``dangle``
------------

``dangle``, :math:`\theta`, places constranits on the lattice parameters :math:`\alpha, \beta`, and :math:`\gamma` as follows:


Triclinic system
^^^^^^^^^^^^^^^^^^

.. math::
   \mathrm{(Type\; 1)} \;\;\; 90^\circ - \theta \leq \alpha, \beta, \gamma < 90^\circ \\
   \mathrm{(Type\; 2)} \;\;\; 90^\circ \leq \alpha, \beta, \gamma \leq 90^\circ + \theta



Monoclinic system
^^^^^^^^^^^^^^^^^^^

.. math::
   90^\circ \leq \beta \leq 90^\circ + \theta


Rhombohedral system
^^^^^^^^^^^^^^^^^^^^^

.. math::
   90^\circ - \theta \leq \alpha \leq 90^\circ + \theta


.. index::
   single: mindist

``mindist``
------------

A mindist matrix consists on ``mindist_1``, ``mindist_2`` ... . For example, in the case of YCo\ :sub:`5` \ (atype = ['Y', 'Co']),
suppose that ``mindist_1`` is  [2.0, 1,8] and ``mindist_2`` is [1.8, 1.5].
The mindist matrix is

.. math::
   \begin{pmatrix}
   2.0 & 1.8 \\
   1.8 & 1.5
   \end{pmatrix}

This means that minimum interatomic distances of Y-Y, Y-Co, and Co-Co are 2.0, 1.8, and 1.5, respectively.
A mindist matrix should be a symmetric matrix.

.. attention::
   mindist matrix: symmetric matrix


.. index::
   single: [VASP]

[VASP] section
==================

.. csv-table::
   :header: Name, Value, Default value, Description
   :widths: auto

   ``kppvol``, int [int ...],  ,  Grid density per Å\ `-3`:sup: of  reciprocal cell in each stage
   ``force_gamma``, bool, ``False`` ,  "If True, force gammma-centered mesh"




.. index::
   single: [QE]

[QE] section
==================

.. csv-table::
   :header: Name, Value, Default value, Description
   :widths: auto

   ``kppvol``, int [int ...],  ,  Grid density per Å\ `-3`:sup: of  reciprocal cell in each stage
   ``qe_infile``, str,  ,  Specify your QE input file name
   ``qe_outfile``, str,  ,  Specify your QE output file name




.. index::
   single: [soiap]

[soiap] section
==================

.. csv-table::
   :header: Name, Value, Default value, Description
   :widths: auto

   ``soiap_infile``, str,  ,  Specify your soiap input file name
   ``soiap_outfile``, str,  ,  Specify your soiap output file name
   ``soiap_cif``,  str,  ,  Specify your CIF-formatted soiap initial structure file name




.. index::
   single: [LAMMPS]

[LAMMPS] section
==================

.. csv-table::
   :header: Name, Value, Default value, Description
   :widths: auto

   ``lammps_infile``, str,  ,  Specify your LAMMPS input file name
   ``lammps_potential``,  "str [str ...], ``None``", ``None`` ,  "Specify your LAMMPS potential, if any"
   ``lammps_outfile``,  str,  ,  Specify your LAMMPS output file name
   ``lammps_data``,  str,  ,  Specify your LAMMPS data file name




.. index::
   single: [BO]

[BO] section
=================

.. csv-table::
   :header: Name, Value, Default value, Description
   :widths: auto

   ``nselect_bo``, int,  ,  Number of structures to be selected at once
   ``score``, "``TS``, ``EI``, ``PI``",  , Acquisition function
   ``num_rand_basis``, int, 0, "If 0: Gaussian process, else: number of basis function"
   ``cdev``, float, 0.001, Cutoff of deviation for standardization
   ``dscrpt``, ``FP`` ,  , Descriptor for structure
   ``fp_rmin``, float, 0.5, Minimum cutoff of *r* in *fingerprint*
   ``fp_rmax``, float, 5.0, Maximum cutoff of *r* in *fingerprint*
   ``fp_npoints``, int, 50, Number of discretized *r* points for each pair in *fingerprint*
   ``fp_sigma``, float, 0.2, Sigma parameter in Gaussian smearing function in Angstrom unit
   ``max_select_bo``, int, 0, Maximum generation
   ``manual_select_bo``, int [int ...], , structure IDs to be selected manually


.. index::
   single: score

``score``
-----------------

- ``TS``: Thompson Sampling
- ``EI``: Expectation Improvement
- ``PI``: Probability of Improvement


.. index::
   single: dscrpt

``dscrpt``
-----------------

- ``FP``: *F*-fingerprint of Oganov and Valle (J. Chem. Phys. 130, 104504 (2009))





.. index::
   single: [LAQA]

[LAQA] section
=================

.. csv-table::
   :header: Name, Value, Default value, Description
   :widths: auto

   ``nselect_laqa``, int,  ,  Number of structures to select at once
   ``weight_laqa``, float, 1.0 ,  weight of bias


.. index::
   single: weight_laqa

``weight_laqa``
-----------------
In LAQA, the score is evaluated by the following equation:

.. math::
   \mathrm{score} = -E + c\frac{F^2}{2\Delta F},

where :math:`c` is ``weight_laqa``, weight of bias.




.. index::
   single: [EA]

[EA] section
=================

.. csv-table::
   :header: Name, Value, Default value, Description
   :widths: auto

   ``n_pop``, int,  ,  Population after second generation
   ``n_crsov``, int, , Number of structure generated by crossover
   ``n_perm``, int, , Number of structure generated by permutation
   ``n_strain``, int, , Number of structure generated by strain
   ``n_rand``, int, , Number of structure generated randomly
   ``n_elite``, int, , Number of elite
   ``fit_reverse``, bool, ``False``, "If False, search minimum value"
   ``n_fittest``, int, ``None`` ,  Number of structure which can survive
   ``slct_func``, "``TNM``, ``RLT``", , Select function
   ``t_size``, int, 3, [Only if slct_func == TNM] Size in tournament selection
   ``a_rlt``, float, 2.0, [Only if slct_func == RLT] Parameter for linear scaling
   ``b_rlt``, float, 1.0, [Only if slct_func == RLT] Parameter for linear scaling
   ``crs_func``, "``OP``, ``TP``", ``OP``, **O**\ne **P**\oint crossover or **T**\wo **P**\oint crossover
   ``crs_lat``, "``equal``, ``random``", ``equal``, How to mix lattice vectors
   ``nat_diff_tole``, int, 4, Tolerance for difference in number of atoms in crossover
   ``ntimes``, int, 1, ntimes permutation
   ``sigma_st``, float, 0.5, Standard deviation for strain
   ``maxcnt_ea``, int, 100, Maximum number of trials in EA
   ``maxgen_ea``, int, 0, Maximum generation
..   ``restart_gen``, int, 0, Restart from specified generation


.. index::
   single: n_pop
   single: n_crsov
   single: n_perm
   single: n_strain
   single: n_rand

``n_pop``, ``n_crsov``, ``n_perm``, ``n_strain``, and ``n_rand``
-----------------------------------------------------------------

Population in first generation is decided by ``tot_struc``. After second generation, population in :math:`n`\ th generation corresponds to ``n_pop``.
For example, ``tot_struc`` = 40, ``n_pop`` = 20 --> population = [40, 20, 20, ...].
You have to set ``n_pop`` to satisfy the following equation:

.. math::
   n_\mathrm{pop} = n_\mathrm{crsov} + n_\mathrm{perm} + n_\mathrm{strain} + n_\mathrm{rand}


.. index::
   single: slct_func
   single: TNM
   single: RLT

``slct_func``
---------------
Available selection functions are tournament (``TNM``) and roulette (``RLT``).


Tournament (``TNM``)
^^^^^^^^^^^^^^^^^^^^^


Roulette (``RLT``)
^^^^^^^^^^^^^^^^^^^^^
If fit_reverse is False, fitness = -fitness.
Linear scaling of fitness:

.. math::
   f_i^{\prime} = \frac{a-b}{f_\mathrm{max} - f_\mathrm{min}}f_i + \frac{b f_\mathrm{max} - a f_\mathrm{min}}{f_\mathrm{max} - f_\mathrm{min}},

where :math:`a, b` are parameters of ``a_rlt`` and ``b_rlt``. :math:`f_i, f_\mathrm{max}, f_\mathrm{min}` are :math:`i` th, maximum, and minimum values of the fitness.

Probability of selecting individual :math:`i` is expressed as:

.. math::
   p_i = \frac{f_i^{\prime}}{\sum_{k}f_k^{\prime}}.



.. index::
   single: simga_EA

``sigma_st``
--------------

In strain operation, lattice vectors :math:`\bm{a}`  are transformed to :math:`\bm{a^\prime}` by applying a strain matrix:

.. math::
   \bm{a^\prime} =    \begin{pmatrix}
                           1+\eta_1 & \frac{1}{2}\eta_6 & \frac{1}{2}\eta_5  \\
                           \frac{1}{2}\eta_6 & 1+\eta_2 & \frac{1}{2}\eta_4  \\
                           \frac{1}{2}\eta_5 & \frac{1}{2}\eta_4 & 1+\eta_3
                      \end{pmatrix} \bm{a},

where :math:`\eta_i` are given by normal distribution with a mean of zero and a standard deviation of ``sigma_st``, :math:`N(0, \sigma^2_\mathrm{st})`.





.. index::
   single: [option]

[option] section
===================

.. csv-table::
   :header: Name, Value, Default value, Description
   :widths: auto

   ``maxcnt``, int,  50,  Maximum number of trials to determine atom positions
   ``stop_chkpt``, int , 0,  Program stops at a specified check point
   ``symprec``, float , 0.001 , Precision for symmetry finding
   ``spgnum``, "``all``, space group number, 0", ``all`` , "Constraint on space group. If all, 1--230. If 0, without space group information "
   ``load_struc_flag``, bool, ``False``, "If True, load initial structures from ``./data/pkl_data/init_struc_data.pkl``"
   ``stop_next_struc``, bool, ``False``, "If True, not submit next structures, but submit next stage and collect results"
   ``recalc``, int [int ...], , "Specify structure IDs if you recalculate or continue optimization" 
   ``append_struc_ea``, bool, ``False``, "If True, append structures by EA"
   ``energy_step_flag``, bool, ``False``, "If True, save energy_step_data in ``./data/pkl_data/energy_step_data.pkl``"
   ``struc_step_flag``, bool, ``False``, "If True, save struc_step_data in ``./data/pkl_data/struc_step_data.pkl``"
   ``fs_step_flag``, bool, ``False``, "If True, save fs_step_data (force and stress) in ``./data/pkl_data/fs_step_data.pkl``"

