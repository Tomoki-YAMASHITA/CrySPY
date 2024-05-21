import numpy as np


class SelectParents:
    '''
    select parents

    # ---------- args
    ranking (list): [ID, ...], ranking of IDs based on fitness without duplication

    # ---------- instance methods
    self.set_tournament(self, t_size)
    self.set_roulette(self, fitness, a_rlt, b_rlt, fit_reverse)
    self.get_parents(n_parent) <-- self._tournament or self._roulette
    '''

    def __init__(self, ranking):
        # ---------- self
        self.ranking = ranking

    def set_tournament(self, t_size):
        '''
        setting for tournament selection
        '''
        # ---------- self
        self.t_size = t_size
        # ---------- set self.get_parents() <-- self._tournament or self._roulette
        self.get_parents = self._tournament

    def _tournament(self, n_parent):
        '''
        tournament selection

        # ---------- args
        n_parent (int): number of parents

        # ---------- return
        parent_id (list): ID of selected parents
        '''
        parent_id = []
        while len(parent_id) < n_parent:
            t_indx = np.random.choice(len(self.ranking), self.t_size, replace=False)
            if parent_id:    # not allow the same parent in crossover
                if parent_id[0] == self.ranking[min(t_indx)]:
                    continue
            parent_id.append(self.ranking[min(t_indx)])
        return parent_id

    def set_roulette(self, fitness, a_rlt, b_rlt, fit_reverse):
        '''
        setting for roulette selection

        # ---------- args
        fitness (dict) : {ID: fitness, ...}
        a_rlt (float): a parameter in linear scaling
        b_rlt (float): b parameter in linear scaling
        fit_reverse (bool): if False, lower fitness is better

        # ---------- instance variables
        self.cum_fit (np.array): cumulative fitness
        '''
        # ----------self
        self.a_rlt = a_rlt
        self.b_rlt = b_rlt
        self.fit_reverse = fit_reverse
        # ---------- calculate cumulative fitness
        fit_deduped = np.array([fitness[i] for i in self.ranking])
        fit_deduped = self._linear_scaling(fit_deduped)
        self.cum_fit = np.cumsum(fit_deduped/fit_deduped.sum())
        # ---------- set self.get_parents <-- self._tournament or self._roulette
        self.get_parents = self._roulette

    def _roulette(self, n_parent):
        '''
        roulette selection

        # ---------- args
        n_parent (int): number of parents

        # ---------- return
        parent_id (list): ID of selected parents
        '''
        # ---------- select parents
        parent_id = []
        while len(parent_id) < n_parent:
            # indx_array: if all false, array is vacant
            indx_array = np.where(self.cum_fit < np.random.rand())[0]
            # select_indx: consider vacant or not
            select_indx = indx_array[-1] + 1 if indx_array.size != 0 else 0
            if parent_id:    # not allow same parent for crossover
                if parent_id[0] == self.ranking[select_indx]:
                    continue
            parent_id.append(self.ranking[select_indx])
        return parent_id

    def _linear_scaling(self, fit_deduped):
        # ------ in case of int
        a = float(self.a_rlt)
        b = float(self.b_rlt)
        # ---------- for fit_reverse
        if not self.fit_reverse:
            fit_deduped = -fit_deduped
        # ---------- scaling
        fmax = float(fit_deduped.max())
        fmin = float(fit_deduped.min())
        # ------ in case the same values
        if fmax == fmin:
            return fit_deduped
        fit_deduped = (a - b)/(fmax - fmin)*fit_deduped + (
            (b*fmax - a*fmin)/(fmax - fmin))
        # ---------- return
        return fit_deduped
