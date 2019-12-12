#!/usr/bin/env python
# -*- coding: utf-8 -*-
# --------------------------------------------------
#
# This code partly includes program codes in COMBO
# (https://github.com/tsudalab/combo)
# which is distributed under the MIT License.
#
# --------------------------------------------------

from __future__ import print_function

import combo
import numpy as np

from ..IO import read_input as rin


# ---------- inheritance
class policy_cryspy(combo.search.discrete.policy):

    def specified_search(self, specified_actions, max_num_probes, num_search_each_probe=1,
                         simulator=None, is_disp=True):
        N = int(num_search_each_probe)
        if int(max_num_probes) * N > len(self.actions):
            raise ValueError('max_num_probes * num_search_each_probe must \
                be smaller than the length of candidates')
        for n in xrange(0, max_num_probes):
            action = self.get_specified_action(N, specified_actions)
            return action

    def get_specified_action(self, N, specified_actions):
        action = specified_actions
        self.actions = self.delete_actions(action)
        return action


def bayes_opt(s_act, descriptors, targets, nselect):
    # ---------- start COMBO part
    # ------ standardization
    # X = combo.misc.centering(descriptors)
    dev = np.std(descriptors, axis=0)
    nonzero_indx = np.where(dev > rin.cdev)[0]
    X = (descriptors[:, nonzero_indx] - np.mean(descriptors[:, nonzero_indx], axis=0)) / dev[nonzero_indx]

    # ------ Declaring the policy by
    policryspy = policy_cryspy(test_X=X)

    # ------ pick up data, already optimized
    actions = policryspy.specified_search(specified_actions=s_act, max_num_probes=1)

    # ------ write
    policryspy.write(actions, -targets)    # Combo trys to find max values --> use "-" (minus)

    # ------ for cryspy.out and log
    out_log(policryspy.history)

    # ------ Bayes_search
    actions = policryspy.bayes_search(max_num_probes=1,
                                      num_search_each_probe=nselect,
                                      score=rin.score, num_rand_basis=rin.num_rand_basis)
    return actions


# ---------- just for cryspy.out
def out_log(history):
    n = history.total_num_search
    index = np.argmax(history.fx[0:n])

    with open('cryspy.out', 'a') as fout:
        fout.write('current best E = {0}\n'.format(-history.fx[index]))
        fout.write('\n')

    print('current best E = {0}\n'.format(-history.fx[index]))
    print('\n')
