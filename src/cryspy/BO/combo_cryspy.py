'''
Modify Policy class in combo3 for CrySPY

 This code partly includes program codes in COMBO
 (https://github.com/tsudalab/combo)
 which is distributed under the MIT License.
'''

from logging import getLogger

import combo
import numpy as np


logger = getLogger('cryspy')

# ---------- inheritance
class Policy_cryspy(combo.search.discrete.policy):

    def specified_search(self, specified_actions, max_num_probes,
                         num_search_each_probe=1,
                         simulator=None, is_disp=True):
        N = int(num_search_each_probe)
        if int(max_num_probes) * N > len(self.actions):
            logger.error('max_num_probes * num_search_each_probe must'
                        ' be smaller than the length of candidates')
            raise SystemExit(1)
        for n in range(0, max_num_probes):
            action = self.get_specified_action(N, specified_actions)
            return action

    def get_specified_action(self, N, specified_actions):
        action = specified_actions
        self.actions = self.delete_actions(action)
        return action

    def bayes_search_cryspy(self, training=None, max_num_probes=None,
                            num_search_each_probe=1,
                            predictor=None, is_disp=True,
                            simulator=None, score='TS', interval=0,
                            num_rand_basis=0):
        if max_num_probes is None:
            max_num_probes = 1
#            simulator = None
        is_rand_expans = False if num_rand_basis == 0 else True
        self.training = self._set_training(training)
        if predictor is None:
            self.predictor = self._init_predictor(is_rand_expans)
        else:
            self.predictor = predictor
        N = int(num_search_each_probe)
        for n in range(max_num_probes):
            if combo.search.utility.is_learning(n, interval):
                self.predictor.fit(self.training, num_rand_basis)
                self.test.Z = self.predictor.get_basis(self.test.X)
                self.training.Z = self.predictor.get_basis(self.training.X)
                self.predictor.prepare(self.training)
            else:
                try:
                    self.predictor.update(self.training, self.new_data)
                except:
                    self.predictor.prepare(self.training)
            if num_search_each_probe != 1:
                combo.search.utility.show_start_message_multi_search(
                    self.history.num_runs, score)
            K = self.config.search.multi_probe_num_sampling
            alpha = self.config.search.alpha
            # ---------- start cryspy
            (action, cryspy_mean,
             cryspy_var, cryspy_score) = self.get_actions_cryspy(score, N, K, alpha)
            return action, cryspy_mean, cryspy_var, cryspy_score
            # ---------- end cryspy
#            if simulator is None:
#                return action
#            t, X = call_simulator(simulator, action)
#            self.write(action, t, X)
#            if is_disp:
#                combo.search.utility.show_search_results(self.history, N)
#        return copy.deepcopy(self.history)

    def get_actions_cryspy(self, mode, N, K, alpha):
#        f = self.get_score_cryspy(mode, self.predictor, self.training, alpha)
        # ---------- start cryspy
        f, cryspy_mean, cryspy_var = self.get_score_cryspy(
            mode, self.predictor, self.training, alpha)
        cryspy_score = f
        # ---------- end cryspy
        temp = np.argmax(f)
        action = self.actions[temp]
        self.actions = self.delete_actions(temp)
        chosed_actions = np.zeros(N, dtype=int)
        chosed_actions[0] = action
        for n in range(1, N):
            f = self.get_marginal_score(mode, chosed_actions[0:n], K, alpha)
            temp = np.argmax(np.mean(f, 0))
            chosed_actions[n] = self.actions[temp]
            self.actions = self.delete_actions(temp)
        return chosed_actions, cryspy_mean, cryspy_var, cryspy_score

    def get_score_cryspy(self, mode, predictor=None, training=None, alpha=1):
        self._set_training(training)
        self._set_predictor(predictor)
        actions = self.actions
        test = self.test.get_subset(actions)
        if mode == 'EI':
            f = combo.search.score.EI(predictor, training, test)
        elif mode == 'PI':
            f = combo.search.score.PI(predictor, training, test)
        elif mode == 'TS':
            f = combo.search.score.TS(predictor, training, test, alpha)
        else:
            logger.error('mode must be EI, PI or TS.')
            raise SystemExit(1)
        # ---------- start cryspy
        cryspy_mean = predictor.get_post_fmean(training, test)
        cryspy_var = predictor.get_post_fcov(training, test)
        return f, cryspy_mean, cryspy_var
        # ---------- end cryspy
#        return f
