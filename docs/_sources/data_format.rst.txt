===========================
Data format
===========================

.. contents:: Contents



Common data
==============
All data on CrySPY are saved in ``./data/pkl_data/`` as pickled files.
You can easily access the data using Python.
If you use Jupyter Notebook, ``/your_cryspy_path/utility/pkl_data.ipynb`` will be useful to analyze the data.

.. seealso::
   :doc:`./utility`




.. index::
   single: init_struc_data.pkl
   single: opt_struc_data.pkl

``init_struc_data.pkl`` / ``opt_struc_data.pkl``
--------------------------------------------------
``init_struc_data.pkl`` (``opt_struc_data.pkl``) includes initial (optimized) structure data, ``init_struc_data`` (``opt_struc_data``).

.. code-block:: Python

   from __future__ import print_function
   import pickle
   with open('init_struc_data.pkl', 'rb') as f
       init_struc_data = pickle.load(f)



.. index::
   single: init_struc_data
   single: opt_struc_data


``init_struc_data`` / ``opt_struc_data``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| Type: dict
| String form: {0: struc1, 1: struc1, ...}
|     The keys are structre ID
|     The values are structure data in pymatgen format
|

.. code-block:: Python

   # ---------- e.g., struc ID 7
   print(init_struc_data[7])



.. index::
   single: rslt_data.pkl

``rslt_data.pkl``
-------------------
``rslt_data.pkl`` includes result data, ``rslt_data``.

.. code-block:: Python

   from __future__ import print_function
   import pickle
   with open('rslt_data.pkl', 'rb') as f
       rslt_data = pickle.load(f)



.. index::
   single: rslt_data

``rslt_data``
^^^^^^^^^^^^^^^^
| Type: DataFrame
|

.. code-block:: Python

   # ---------- e.g., first 10 data
   print(rslt_data.head(10))



.. index::
   single: kpt_data.pkl

``kpt_data.pkl``
------------------
``kpt_data.pkl`` includes k-point data, ``kpt_data``.

.. code-block:: Python

   from __future__ import print_function
   import pickle
   with open('kpt_data.pkl', 'rb') as f
       kpt_data = pickle.load(f)



.. index::
   single: kpt_data

``kpt_data``
^^^^^^^^^^^^^^^
| Type: dict
| String form: {0: [[2, ,2 ,2], [4, 4, 4,], [6, 6, 6,], ...], 1: ...}
|     The keys are structure ID
|     The values indicate k-mesh in each stage as list
|

.. code-block:: Python

   # ---------- e.g., k-mesh of ID 7, stage 2
   # kpt_data[ID][stage]
   # kpt_data[ID][0] <-- stage 1
   # kpt_data[ID][1] <-- stage 2
   print(kpt_data[7][1])




Random Search
==============

.. index::
   single: RS_id_data.pkl

``RS_id_data.pkl``
------------------
``RS_id_data.pkl`` includes

- ``next_id``: next structure ID to calculate
- ``id_done``: finished structure ID

.. code-block:: Python

   import pickle
   with open('RS_id_data.pkl', 'rb') as f
       next_id, id_done = pickle.load(f)



.. index::
   single: next_id

``next_id``
^^^^^^^^^^^^^^^
| Type: int
| String form: 5
|


.. index::
   single: id_done(RS)

``id_done``
^^^^^^^^^^^^^^^
| Type: 1d array
| String form: [0 1 2 3 4]
|




Bayesian Optimization
======================

.. index::
   single: BO_id_data.pkl

``BO_id_data.pkl``
------------------
``BO_id_data.pkl`` includes

- ``gen``: generation
- ``non_error_id``: non-error structure ID
- ``id_to_calc``: structure ID to calculate in the current generation
- ``id_done``: finished structure ID

.. code-block:: Python

   import pickle
   with open('BO_id_data.pkl', 'rb') as f
       gen, non_error_id, id_to_calc, id_done = pickle.load(f)


.. index::
   single: gen

``gen``
^^^^^^^^^^^^^^^
| Type: int
| String form: 1
|


.. index::
   single: non_error_id

``non_error_id``
^^^^^^^^^^^^^^^^^^
| Type: 1d array
| String form: [0 1 2 3 4 5 6 7 8 9]
|


.. index::
   single: id_to_calc(BO)

``id_to_calc``
^^^^^^^^^^^^^^^
| Type: 1d array
| String form: [8 6 4 1 2]
|


.. index::
   single: id_done(BO)

``id_done``
^^^^^^^^^^^^^^^
| Type: 1d array
| String form: [0 9 3 5 7]
|



.. index::
   single: BO_data.pkl

``BO_data.pkl``
------------------
``BO_data.pkl`` includes

- ``descriptors``: descriptor data
- ``targets``: target(=energy) data

.. code-block:: Python

   import pickle
   with open('BO_data.pkl', 'rb') as f
       descriptors, targets = pickle.load(f)




.. index::
   single: descriptors

``descriptors``
^^^^^^^^^^^^^^^
| Type: 2d array
| String form: [[descriptor array of ID 0], [descriptor array of ID 1], [descriptor array of ID 3], ....]
|     len(``descriptors``) = len(``non_error_id``)
|     If your calculation for ID 2 failed, the descriptor data of ID 2 is deleted like this example.
|

.. code-block:: Python

   # ---------- how to access descriptor of ID 3
   # descriptors[3] does not always correspond to the data of ID 3!
   descriptors[np.where(non_error_id == 3)[0][0]]




.. index::
   single: targets

``targets``
^^^^^^^^^^^^^^^
| Type: 1d array
| String form: [-10.45, -8.789, ....]
|     len(``targets``) = len(``id_done``)
|     The order of ID in ``targets`` follows the order of ``id_done``
|
| e.g.,
| id_done = [7, 3, 0, 1, ...]
| targets = [energy of ID 7, energy of ID 3, energy of ID 0, energy of ID 1, ...]
|

.. code-block:: Python

   # ---------- how to access target of ID 3
   # targets[3] does not correspond to the data of ID 3!
   targets[np.where(id_done == 3)[0][0]]





LAQA
======================

.. index::
   single: LAQA_id_data.pkl

``LAQA_id_data.pkl``
-----------------------
``LAQA_id_data.pkl`` includes

- ``id_to_calc``: structure ID to calculate in the current generation
- ``id_select_hist``: history of ID selection
- ``id_done``: finished structure ID

.. code-block:: Python

   import pickle
   with open('LAQA_id_data.pkl', 'rb') as f
       id_to_calc, id_select_hist, id_done = pickle.load(f)


.. index::
   single: id_to_calc(LAQA)

``id_to_calc``
^^^^^^^^^^^^^^^
| Type: list
| String form: [8, 6, 4, 1, 2]
|


.. index::
   single: id_select_hist

``id_select_hist``
^^^^^^^^^^^^^^^^^^^
| Type: list
| String form: [[5, 3, 9, 0, 7], ..., [8, 6, 4, 1, 2]]
|    [[list of first selection], [list of second selection], ...]
|

.. note::
   ``id_select_hist`` does not include 0th selection (all initial structures), start with 1st selection.




.. index::
   single: id_done(LAQA)

``id_done``
^^^^^^^^^^^^^^^
| Type: list
| String form: [5, 3, 9, 0, 7]
|



.. index::
   single: LAQA_data.pkl

``LAQA_data.pkl``
------------------
``LAQA_data.pkl`` includes

- ``tot_step_select``: total number of optimization steps in each selection
- ``LAQA_step``: number of optimization steps in each ID
- ``LAQA_struc``: list of structure data in each ID
- ``LAQA_energy``: list of energy data in each ID
- ``LAQA_bias``:
- ``LAQA_score``:

.. code-block:: Python

   import pickle
   with open('LAQA_data.pkl', 'rb') as f
       tot_step_select, LAQA_step, LAQA_struc, LAQA_energy, LAQA_bias, LAQA_score = pickle.load(f)


.. index::
   single: tot_step_select

``tot_step_select``
^^^^^^^^^^^^^^^^^^^^^
| Type: list
| String form: [2000, 200, 200, ...]
|     len(``tot_step_select``) = len(``id_select_hist``) + 1
|     ``tot_step_select`` includes 0th selection
|     [0th, 1st, 2nd, ...]
|


.. code-block:: Python

   # ---------- total number of optimization steps (all steps)
   print('Total steps: {}'.format(sum(tot_step_select)))
   # ---------- up to 5 selection. Note that tot_step_select includes 0th selection
   print('Number of steps up to 5 selection: {}'.format(sum(tot_step_select[:5+1])))

.. note::
   ``tot_step_select`` includes 0th selection (all initial structures) unlike ``id_select_hist``



.. index::
   single: LAQA_step

``LAQA_step``
^^^^^^^^^^^^^^^^^^^^^
| Type: dict
| String form: {0: [20, 7], 1:[20, 20, 20, 5], ...}
|     The keys are structure ID
|     The values are list of number of optimization steps
|

.. code-block:: Python

   print(LAQA_step[7])
   # ---------- total steps in ID 7
   print('Total number of optimization steps in ID 7: {}'.format(sum(LAQA_step[7])))



.. index::
   single: LAQA_struc

``LAQA_struc``
^^^^^^^^^^^^^^^^^^^^^
| Type: dict
| String form: {0: [list of structures], 1:[list of structures], ...}
|     The keys are structure ID
|     The values are list of structure data in pymatgen format
|     len(LAQA_struc[7]) == len(LAQA_step[7])
|     If the optimization finished, LAQA_struc[7][-1] is equal to opt_struc_data[7]
|

.. code-block:: Python

   # ---------- latest structure of ID 7
   print(LAQA_struc[7][-1])



.. index::
   single: LAQA_energy

``LAQA_energy``
^^^^^^^^^^^^^^^^^^^^^
| Type: dict
| String form: {0: [-3.287, -3.330], 1:[-3.105, -3.194, -3.233, -3.347], ...}
|     The keys are structure ID
|     The values are list of energy data
|     len(LAQA_energy[7]) == len(LAQA_step[7])
|

.. code-block:: Python

   # ---------- energy list of ID 7
   print(LAQA_energy[7])
   # ---------- latest energy of ID 7
   print(LAQA_energy[7][-1])



.. index::
   single: LAQA_bias

``LAQA_bias``
^^^^^^^^^^^^^^^^^^^^^
| Type: dict
| String form: {0: [0.059, 0.003], 1:[0.501, 0.210, 0.984, 0.758], ...}
|     The keys are structure ID
|     The values are list of bias data
|     len(LAQA_bias[7]) == len(LAQA_step[7])
|

.. code-block:: Python

   # ---------- bias list of ID 7
   print(LAQA_bias[7])
   # ---------- latest bias of ID 7
   print(LAQA_bias[7][-1])


.. index::
   single: LAQA_score

``LAQA_score``
^^^^^^^^^^^^^^^^^^^^^
| Type: dict
| String form: {0: [inf, 3.346, -inf], 1:[3.606, 3.404, 4.217, -inf], ...}
|     The keys are structure ID
|     The values are list of score data
|     len(LAQA_score[7]) == len(LAQA_step[7]) + 1
|     ``LAQA_score`` includes 0th score (= plus infinity)
|     If the optimization finished, -inf is appended to the score list
|

.. code-block:: Python

   # ---------- score list of ID 7
   print(LAQA_score[7])
   # ---------- latest score of ID 7
   print(LAQA_score[7][-1])


.. note::
   ``LAQA_score`` includes 0th score (= inf) unlike ``LAQA_energy`` and ``LAQA_bias``, so len(LAQA_score[7]) is not equal to len(LAQA_energy[7]).







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
| Type: dict
| String form: {0:  [ [array(stage1, step1), array(stage1, step2), ...], [array(stage2, step1), array(stage2, step2), ...], ... ]}
|     The keys are structure ID
|     The values are energy-step array
|

.. code-block:: Python

   # energy_step_data[ID][stage][step]
   # energy_step_data[ID][0] <-- stage 1
   # energy_step_data[ID][1] <-- stage 2
   #
   # ---------- energy-step data of ID 7, stage 2, step 8
   energy_step_data[7][2-1][8-1]


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
| Type: dict
| String form: {0:  [ [ (stage1, step1), (stage1, step2), ...], [(stage2, step1), (stage2, step2), ...], ...]}
|     The keys are structure ID
|     The values are structure-step list
|

.. code-block:: Python

   # struc_step_data[ID][stage][step]
   # struc_step_data[ID][0] <-- stage 1
   # struc_step_data[ID][1] <-- stage 2
   #
   # ---------- structure-step data of ID 7, stage 2, step 8
   sturc_step_data[7][2-1][8-1]


.. note::
   stage and step start from 1 unlike ID





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
| Type: dict
| String form: {0:  [ [array(stage1, step1), array(stage1, step2), ...], [array(stage2, step1), array(stage2, step2), ...], ... ]}
|     The keys are structure ID
|     The values are force-step array
|

.. code-block:: Python

   # force_step_data[ID][stage][step]
   # force_step_data[ID][0] <-- stage 1
   # force_step_data[ID][1] <-- stage 2
   #
   # ---------- force-step data of ID 7, stage 2, step 8
   force_step_data[7][2-1][8-1]


.. note::
   stage and step start from 1 unlike ID



.. index::
   single: stress_step_data

``stress_step_data``
^^^^^^^^^^^^^^^^^^^^^
| Type: dict
| String form: {0:  [ [array(stage1, step1), array(stage1, step2), ...], [array(stage2, step1), array(stage2, step2), ...], ... ]}
|     The keys are structure ID
|     The values are stress-step array
|

.. code-block:: Python

   # stress_step_data[ID][stage][step]
   # stress_step_data[ID][0] <-- stage 1
   # stress_step_data[ID][1] <-- stage 2
   #
   # ---------- stress-step data of ID 7, stage 2, step 8
   stress_step_data[7][2-1][8-1]


.. note::
   stage and step start from 1 unlike ID
