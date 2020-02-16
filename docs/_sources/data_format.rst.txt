===========================
Data format
===========================

.. contents:: Contents



Common data
==============
All data on CrySPY are saved in ``./data/pkl_data/`` as pickled files.
You can easily access the data using Python.
If you use Jupyter Notebook, ``CrySPY-x.x.x/utility/pkl_data.ipynb`` would be useful to analyze the data.

.. seealso::
   :doc:`./utility`


.. index::
   single: init_struc_data.pkl

``init_struc_data.pkl`` 
---------------------------
``init_struc_data.pkl`` includes initial structure data, ``init_struc_data``.

.. code-block:: Python

   import pickle
   with open('init_struc_data.pkl', 'rb') as f
       init_struc_data = pickle.load(f)


.. index::
   single: init_struc_data


``init_struc_data`` 
^^^^^^^^^^^^^^^^^^^^^^
- Type: dict

    + The keys are structre IDs
    + The values are structure data in pymatgen format

- String form: {0: struc1, 1: struc1, ...}

.. code-block:: Python

   # ---------- e.g., struc ID 7
   print(init_struc_data[7])


.. index::
   single: opt_struc_data.pkl

``opt_struc_data.pkl`` 
---------------------------
``opt_struc_data.pkl`` includes optimized structure data, ``opt_struc_data``.

.. code-block:: Python

   import pickle
   with open('opt_struc_data.pkl', 'rb') as f
       opt_struc_data = pickle.load(f)


.. index::
   single: opt_struc_data


``opt_struc_data`` 
^^^^^^^^^^^^^^^^^^^^^^
- Type: dict

    + The keys are structre IDs
    + The values are structure data in pymatgen format

- String form: {0: struc1, 1: struc1, ...}

.. code-block:: Python

   # ---------- e.g., struc ID 7
   print(opt_struc_data[7])




.. index::
   single: rslt_data.pkl

``rslt_data.pkl``
-------------------
``rslt_data.pkl`` includes result data, ``rslt_data``.

.. code-block:: Python

   import pickle
   with open('rslt_data.pkl', 'rb') as f
       rslt_data = pickle.load(f)


.. index::
   single: rslt_data

``rslt_data``
^^^^^^^^^^^^^^^^

- Type: DataFrame

    + Row labels are structure IDs

- String form:

.. code-block::

      Spg_num  Spg_sym  Spg_num_opt Spg_sym_opt  E_eV_atom  Magmom   Opt
   0      168       P6          191      P6/mmm  -3.826916     NaN  done
   1       95   P4_322           95      P4_322  -3.978478     NaN  done
   2      191   P6/mmm          191      P6/mmm  -2.289173     NaN  done
   3      113  P-42_1m          113     P-42_1m  -3.773191     NaN  done
   4      128   P4/mnc          123      P4/mmm  -3.296904     NaN  done

.. code-block:: Python

   # ---------- e.g., first 10 data
   print(rslt_data.head(10))



.. index::
   single: kpt_data.pkl

``kpt_data.pkl``
------------------
``kpt_data.pkl`` includes k-point data, ``kpt_data``.

.. code-block:: Python

   import pickle
   with open('kpt_data.pkl', 'rb') as f
       kpt_data = pickle.load(f)


.. index::
   single: kpt_data

``kpt_data``
^^^^^^^^^^^^^^^

- Type: dict

    + The keys are structure IDs
    + The values are list of k-mesh in each stage

- String form: {0: [[2, ,2 ,2], [4, 4, 4,], [6, 6, 6,], ...], 1: ...}


.. code-block:: Python

   # ---------- e.g., k-mesh of ID 7, stage 2
   # kpt_data[ID][stage]
   # kpt_data[ID][0] <-- stage 1
   # kpt_data[ID][1] <-- stage 2
   is = 2
   print(kpt_data[7][is-1])




Random Search
==============

.. index::
   single: RS_id_data.pkl

``RS_id_data.pkl``
------------------
``RS_id_data.pkl`` includes

- ``id_queueing``: queueing IDs
- ``id_running``: running IDs

.. code-block:: Python

   import pickle
   with open('RS_id_data.pkl', 'rb') as f
       id_queueing, id_running = pickle.load(f)



.. index::
   single: id_queueing (RS)

``id_queueing``
^^^^^^^^^^^^^^^

- Type: list
- String form: [5, 6, 7, 8, 9]



.. index::
   single: id_running (RS)

``id_running``
^^^^^^^^^^^^^^^

- Type: list
- String form: [0, 1, 2, 3, 4]





Bayesian Optimization
======================

.. index::
   single: BO_id_data.pkl

``BO_id_data.pkl``
------------------
``BO_id_data.pkl`` includes

- ``n_selection``: number of selection
- ``id_queueing``: queueing IDs
- ``id_running``: running IDs
- ``id_select_hist``: history of ID selection


.. code-block:: Python

   import pickle
   with open('BO_id_data.pkl', 'rb') as f
       n_selection, id_queueing, id_running, id_select_hist = pickle.load(f)


.. index::
   single: n_selection

``n_selection``
^^^^^^^^^^^^^^^

- Type: int
- String form: 1


.. index::
   single: id_queueing (BO)

``id_queueing``
^^^^^^^^^^^^^^^

- Type: list
- String form: [5, 6, 7, 8, 9]



.. index::
   single: id_running (BO)

``id_running``
^^^^^^^^^^^^^^^

- Type: list
- String form: [0, 1, 2, 3, 4]


.. index::
   single: id_select_hist (BO)

``id_select_hist``
^^^^^^^^^^^^^^^^^^^

- Type: list
- String form: [[5, 3, 9, 0, 7], ..., [8, 6, 4, 1, 2]]

    + [[list of first selection], [list of second selection], ...]




.. index::
   single: BO_data.pkl

``BO_data.pkl``
------------------
``BO_data.pkl`` includes

- ``init_dscrpt_data``: descriptor data for initial structures
- ``opt_dscrpt_data``: descriptor data for optimized structures
- ``bo_mean``: predictive mean in BO
- ``bo_var``: predictive variance in BO 
- ``bo_score``: score in BO


.. code-block:: Python

   import pickle
   with open('BO_data.pkl', 'rb') as f
       init_dscrpt_data, opt_dscrpt_data, bo_mean, bo_var, bo_score = pickle.load(f)




.. index::
   single: init_dscrpt_data

``init_dscrpt_data``
^^^^^^^^^^^^^^^^^^^^^^

- Type: dict

  + The keys are structre IDs
  + The values are descriptor data for initial structures in numpy.ndarray

- String form: {0: array([descriptor array of ID 0]), 1: array([descriptor array of ID 1]), ...}


.. code-block:: Python

   # ---------- init_dscrpt_data of ID 3
   init_dscrpt_data[3]



.. index::
   single: opt_dscrpt_data

``opt_dscrpt_data``
^^^^^^^^^^^^^^^^^^^^^^

- Type: dict

  + The keys are structre IDs
  + The values are descriptor data for optimized structures in numpy.ndarray

- String form: {0: array([descriptor array of ID 0]), 1: array([descriptor array of ID 1]), ...}


.. code-block:: Python

   # ---------- opt_dscrpt_data of ID 3
   opt_dscrpt_data[3]



.. index::
   single: bo_mean

``bo_mean``
^^^^^^^^^^^^^^^^^^^^^^

- Type: dict

  + The keys are selection No.
  + The values are dict of predictive mean

- String form: {2: {0: 3.93, 1: 3.92, 2: 3.94, ...}, 3: {...}, ...}

.. code-block:: Python

   # ---------- predictive mean data for each candidate at 2nd selection
   #            (1st selection is random)
   bo_mean[2]


.. index::
   single: bo_var

``bo_var``
^^^^^^^^^^^^^^^^^^^^^^

- Type: dict

  + The keys are selection No.
  + The values are dict of predictive variance

- String form: {2: {0: 0.014, 1: 0.013, 2: 0.018, ...}, 3: {...}, ...}

.. code-block:: Python

   # ---------- predictive variance data for each candidate at 2nd selection
   #            (1st selection is random)
   bo_var[2]



.. index::
   single: bo_score

``bo_score``
^^^^^^^^^^^^^^^^^^^^^^

- Type: dict

  + The keys are selection No.
  + The values are dict of score

- String form: {2: {0: 4.076, 1: 3.995, 2: 4.010, ...}, 3: {...}, ...}

.. code-block:: Python

   # ---------- score data for each candidate at 2nd selection
   #            (1st selection is random)
   bo_score[2]




LAQA
======================

.. index::
   single: LAQA_id_data.pkl

``LAQA_id_data.pkl``
-----------------------
``LAQA_id_data.pkl`` includes

- ``id_queueing``: queueing IDs
- ``id_running``: running IDs
- ``id_select_hist``: history of ID selection

.. code-block:: Python

   import pickle
   with open('LAQA_id_data.pkl', 'rb') as f
       id_queueing, id_running, id_select_hist = pickle.load(f)


.. index::
   single: id_queueing (LAQA)

``id_queueing``
^^^^^^^^^^^^^^^

- Type: list
- String form: [5, 6, 7, 8, 9]



.. index::
   single: id_running (LAQA)

``id_running``
^^^^^^^^^^^^^^^

- Type: list
- String form: [0, 1, 2, 3, 4]


.. index::
   single: id_select_hist (LAQA)

``id_select_hist``
^^^^^^^^^^^^^^^^^^^

- Type: list
- String form: [[5, 3, 9, 0, 7], ..., [8, 6, 4, 1, 2]]

    + [[list of first selection], [list of second selection], ...]


.. note::
   ``id_select_hist`` does not include 0th selection (all initial structures), start with 1st selection.




.. index::
   single: LAQA_data.pkl

``LAQA_data.pkl``
------------------
``LAQA_data.pkl`` includes

- ``tot_step_select``: total number of optimization steps in each selection
- ``laqa_step``: number of optimization steps in each ID
- ``laqa_struc``: list of structure data in each ID
- ``laqa_energy``: list of energy data in each ID
- ``laqa_bias``: list of bias data in each ID
- ``laqa_score``: list of score data in each ID

.. code-block:: Python

   import pickle
   with open('LAQA_data.pkl', 'rb') as f
       tot_step_select, laqa_step, laqa_struc, laqa_energy, laqa_bias, laqa_score = pickle.load(f)


.. index::
   single: tot_step_select

``tot_step_select``
^^^^^^^^^^^^^^^^^^^^^

- Type: list

    + len(``tot_step_select``) = len(``id_select_hist``) + 1
    + ``tot_step_select`` includes 0th selection

- String form: [2000, 200, 200, ...]

   + [0th, 1st, 2nd, ...]



.. code-block:: Python

   # ---------- total number of optimization steps (all steps)
   print('Total steps: {}'.format(sum(tot_step_select)))
   tot_step_select   # ---------- total number of optimization steps (all steps)

.. note::
   ``tot_step_select`` includes 0th selection (all initial structures) unlike ``id_select_hist``



.. index::
   single: laqa_step

``laqa_step``
^^^^^^^^^^^^^^^^^^^^^

- Type: dict

    + The keys are structure ID
    + The values are list of number of optimization steps

- String form: {0: [20, 7], 1:[20, 20, 20, 5], ...}


.. code-block:: Python

   print(laqa_step[7])
   # ---------- total steps in ID 7
   print('Total number of optimization steps in ID 7: {}'.format(sum(LAQA_step[7])))



.. index::
   single: laqa_struc

``laqa_struc``
^^^^^^^^^^^^^^^^^^^^^

- Type: dict

    +  The keys are structure ID
    +  The values are list of structure data in pymatgen format

- String form: {0: [list of structures], 1:[list of structures], ...}

    + len(laqa_struc[7]) == len(laqa_step[7])

| Latest structure data in each job are save in laqa_struc.
| If the optimization finished, laqa_struc[7][-1] is equal to opt_struc_data[7]
|
| For example,
| number of iteration for optimization = 5 (NSW = 5 in VASP input)  
|     5 opt. step --> save latest struc. --> 5 opt. step --> save latest struc. --> ...
|
| laqa_step[ID] = [5, 5, 5, ...]
| laqa_struc[ID] = [a struc_data, a struc_data, ...]
|
| So, 4 structure data are discarded in each job.
| If you want to save full structure data step by step,
| use ``struc_step_flag = True`` in cryspy.in.


.. code-block:: Python

   # ---------- latest structure of ID 7
   print(laqa_struc[7][-1])



.. index::
   single: laqa_energy

``laqa_energy``
^^^^^^^^^^^^^^^^^^^^^

- Type: dict

    + The keys are structure ID
    + The values are list of energy data

- String form: {0: [-3.287, -3.330], 1:[-3.105, -3.194, -3.233, -3.347], ...}

    + len(laqa_energy[7]) == len(laqa_step[7])

| Latest energy data in each job are save in laqa_energy.
|
| For example,
| number of iteration for optimization = 5 (NSW = 5 in VASP input)  
|     5 opt. step --> save latest energy --> 5 opt. step --> save latest energy --> ...
|
| laqa_step[ID] = [5, 5, 5, ...]
| laqa_energy[ID] = [an energy_data, an energy_data, ...]
|
| So, 4 energy data are discarded in each job.
| If you want to save full energy data step by step,
| use ``energy_step_flag = True`` in cryspy.in.



.. code-block:: Python

   # ---------- energy list of ID 7
   print(laqa_energy[7])




.. index::
   single: laqa_bias

``laqa_bias``
^^^^^^^^^^^^^^^^^^^^^

- Type: dict

    + The keys are structure ID
    + The values are list of bias data

- String form: {0: [0.059, 0.003], 1:[0.501, 0.210, 0.984, 0.758], ...}

    + len(laqa_bias[7]) == len(laqa_step[7])


.. code-block:: Python

   # ---------- bias list of ID 7
   print(laqa_bias[7])
   # ---------- latest bias of ID 7
   print(laqa_bias[7][-1])


.. index::
   single: laqa_score

``laqa_score``
^^^^^^^^^^^^^^^^^^^^^

- Type: dict

    + The keys are structure ID
    + The values are list of score data

- String form: {0: [inf, 3.346, -inf], 1:[3.606, 3.404, 4.217, -inf], ...}

    + len(laqa_score[7]) == len(laqa_step[7]) + 1

|     ``laqa_score`` includes 0th score (= plus infinity)
|     If the optimization finished, -inf is appended to the score list
|

.. code-block:: Python

   # ---------- score list of ID 7
   print(laqa_score[7])
   # ---------- latest score of ID 7
   print(laqa_score[7][-1])


.. note::
   ``laqa_score`` includes 0th score (= inf) unlike ``laqa_energy`` and ``LAQA_bias``, so len(LAQA_score[7]) is not equal to len(laqa_energy[7]).




Evolutionary algorithm
======================

.. index::
   single: EA_id_data.pkl

``EA_id_data.pkl``
------------------
``EA_id_data.pkl`` includes

- ``gen``: current generation
- ``id_queueing``: queueing IDs
- ``id_running``: running IDs


.. code-block:: Python

   import pickle
   with open('EA_id_data.pkl', 'rb') as f
       gen, id_queueing, id_running = pickle.load(f)


.. index::
   single: gen

``gen``
^^^^^^^^^^^^^^^

- Type: int
- String form: 1



.. index::
   single: id_queueing (EA)

``id_queueing``
^^^^^^^^^^^^^^^

- Type: list
- String form: [5, 6, 7, 8, 9]



.. index::
   single: id_running (EA)

``id_running``
^^^^^^^^^^^^^^^

- Type: list
- String form: [0, 1, 2, 3, 4]




.. index::
   single: EA_data.pkl

``EA_data.pkl``
------------------
``EA_data.pkl`` includes

- ``elite_struc``: elite structure data
- ``elite_fitness``: fitness of elite structures
- ``ea_info``: information on generational changes
- ``ea_origin``: information on origins (parents)

.. code-block:: Python

   import pickle
   with open('EA_data.pkl', 'rb') as f
       elite_struc, elite_fitness, ea_info, ea_origin = pickle.load(f)



.. index::
   single: elite_struc

``elite_struc``
^^^^^^^^^^^^^^^

- Type: dict

    + The keys are elite structre IDs
    + The values are elite structure data in pymatgen format

- String form: {0: struc0, 4: struc4, ...}



.. index::
   single: elite_fitness

``elite_fitness``
^^^^^^^^^^^^^^^^^^^

- Type: dict

    + The keys are elite structre IDs
    + The values are fitness of elite structures

- String form: {4: -4.101055417556523, 0: -4.061872594010355}



.. index::
   single: ea_info

``ea_info``
^^^^^^^^^^^^^^^

- Type: DataFrame
- String form:

.. code-block::

    Gen  Population  Crossover  Permutation  Strain  Random  Elite crs_func crs_lat slct_func
      1          10          0            0       0      10      0       OP   equal       TNM
      2          10          5            0       3       2      2       OP   equal       TNM



.. index::
   single: ea_origin

``ea_origin``
^^^^^^^^^^^^^^^

- Type: DataFrame
- String form:

.. code-block::

    Gen  Struc_ID  Operation  Parent
      1         0     random    None
      1         1     random    None
      1         2     random    None
      1         3     random    None
      1         4     random    None
      1         5     random    None
      1         6     random    None
      1         7     random    None
      1         8     random    None
      1         9     random    None
      2        10  crossover  (9, 5)
      2        11  crossover  (9, 4)
      2        12  crossover  (7, 4)
      2        13  crossover  (4, 5)
      2        14  crossover  (9, 7)
      2        15     strain    (0,)
      2        16     strain    (4,)
      2        17     strain    (9,)
      2        18     random    None
      2        19     random    None
      2         4      elite   elite
      2         0      elite   elite





Optional data
======================

.. index::
   single: energy_step_data.pkl

``energy_step_data.pkl``
---------------------------
``energy_step_data.pkl`` includes energy_step_data

.. code-block:: Python

   import pickle
   with open('energy_step_data.pkl', 'rb') as f
       energy_step_data = pickle.load(f)


.. index::
   single: energy_step_data

``energy_step_data``
^^^^^^^^^^^^^^^^^^^^^

- Type: dict

    + The keys are structure ID
    + The values are energy-step numpy.ndarray

- String form: {0:  [ [array(stage1, step1), array(stage1, step2), ...], [array(stage2, step1), array(stage2, step2), ...], ... ]}


.. code-block:: Python

   # energy_step_data[ID][stage][step]
   # energy_step_data[ID][0] <-- stage 1
   # energy_step_data[ID][1] <-- stage 2
   #
   # in LAQA
   # energy_step_data[ID][0] <-- 1st selection
   # energy_step_data[ID][1] <-- 2nd selection
   # ---------- energy-step data of ID 7, stage 2, step 5
   energy_step_data[7][2-1][5-1]


.. note::
   stage and step start from 1 unlike ID




.. index::
   single: struc_step_data.pkl

``struc_step_data.pkl``
---------------------------
``struc_step_data.pkl`` includes struc_step_data

.. code-block:: Python

   import pickle
   with open('struc_step_data.pkl', 'rb') as f
       struc_step_data = pickle.load(f)


.. index::
   single: struc_step_data

``struc_step_data``
^^^^^^^^^^^^^^^^^^^^^

- Type: dict

    + The keys are structure ID
    + The values are structure-step list

- String form: {0:  [ [ (stage1, step1), (stage1, step2), ...], [(stage2, step1), (stage2, step2), ...], ...]}


.. code-block:: Python

   # struc_step_data[ID][stage][step]
   # struc_step_data[ID][0] <-- stage 1
   # struc_step_data[ID][1] <-- stage 2
   #
   #
   # in LAQA
   # struc_step_data[ID][0] <-- 1st selection
   # struc_step_data[ID][1] <-- 2nd selection
   # ---------- structure-step data of ID 7, stage 2, step 5
   sturc_step_data[7][2-1][5-1]


.. note::
   The stage and step start from 1 unlike ID





.. index::
   single: fs_step_data.pkl

``fs_step_data.pkl``
---------------------------
``fs_step_data.pkl`` includes

- ``force_step_data``: force-step data
- ``stress_step_data``: stress-step data

.. code-block:: Python

   import pickle
   with open('fs_step_data.pkl', 'rb') as f
       force_step_data, stress_step_data = pickle.load(f)


.. index::
   single: force_step_data

``force_step_data``
^^^^^^^^^^^^^^^^^^^^^

- Type: dict

    + The keys are structure ID
    + The values are force-step numpy.ndarray

- String form: {0:  [ [array(stage1, step1), array(stage1, step2), ...], [array(stage2, step1), array(stage2, step2), ...], ... ]}

.. code-block:: Python

   # force_step_data[ID][stage][step]
   # force_step_data[ID][0] <-- stage 1
   # force_step_data[ID][1] <-- stage 2
   #
   # in LAQA
   # force_step_data[ID][0] <-- 1st selection
   # force_step_data[ID][1] <-- 2nd selection
   # ---------- force-step data of ID 7, stage 2, step 5
   force_step_data[7][2-1][5-1]


.. note::
   The stage and step start from 1 unlike ID



.. index::
   single: stress_step_data

``stress_step_data``
^^^^^^^^^^^^^^^^^^^^^

- Type: dict

    + The keys are structure ID
    + The values are stress-step numpy.ndarray

- String form: {0:  [ [array(stage1, step1), array(stage1, step2), ...], [array(stage2, step1), array(stage2, step2), ...], ... ]}

.. code-block:: Python

   # stress_step_data[ID][stage][step]
   # stress_step_data[ID][0] <-- stage 1
   # stress_step_data[ID][1] <-- stage 2
   #
   # in LAQA
   # stress_step_data[ID][0] <-- 1st selection
   # stress_step_data[ID][1] <-- 2nd selection
   # ---------- stress-step data of ID 7, stage 2, step 5
   stress_step_data[7][2-1][5-1]


.. note::
   The stage and step start from 1 unlike ID
