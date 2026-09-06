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

    def __init__(self, ranking, rng=None):
        # ---------- self
        self.ranking = ranking
        self.rng = rng

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
        # ---------- initialize rng
        if self.rng is None:
            rng = np.random.default_rng()
        else:
            rng = self.rng
        # ---------- select parents
        parent_id = []
        for _ in range(n_parent):
            available_indices = [
                index
                for index, candidate_id in enumerate(self.ranking)
                if candidate_id not in parent_id
            ]
            t_size_eff = min(self.t_size, len(available_indices))
            t_indx = rng.choice(available_indices, t_size_eff, replace=False)
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
        self.weights (np.array): roulette weights
        '''
        # ----------self
        self.a_rlt = a_rlt
        self.b_rlt = b_rlt
        self.fit_reverse = fit_reverse
        # ---------- calculate weights
        fit_deduped = np.array([fitness[i] for i in self.ranking], dtype=float)
        self.weights = self._linear_scaling(fit_deduped)
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
        # ---------- initialize rng
        if self.rng is None:
            rng = np.random.default_rng()
        else:
            rng = self.rng
        # ---------- select parents
        parent_id = []
        for _ in range(n_parent):
            available_indices = [
                index
                for index, candidate_id in enumerate(self.ranking)
                if candidate_id not in parent_id
            ]
            available_weights = self.weights[available_indices]
            available_weights /= available_weights.max()
            select_indx = rng.choice(
                available_indices,
                p=available_weights / available_weights.sum(),
            )
            parent_id.append(self.ranking[select_indx])
        return parent_id

    def _linear_scaling(self, fit_deduped):
        # ------ in case of int
        a = float(self.a_rlt)
        b = float(self.b_rlt)
        # ---------- check input
        if fit_deduped.size == 0 or not np.all(np.isfinite(fit_deduped)):
            raise ValueError('Roulette requires non-empty, finite fitness')
        # ---------- for fit_reverse
        if not self.fit_reverse:
            fit_deduped = -fit_deduped
        # ---------- scaling
        fmax = float(fit_deduped.max())
        fmin = float(fit_deduped.min())
        # ------ in case the same values
        if fmax == fmin:
            return np.ones_like(fit_deduped, dtype=float)
        fit_deduped = b + (a - b) * (
            (fit_deduped - fmin) / (fmax - fmin)
        )
        if not np.all(np.isfinite(fit_deduped)) or np.any(fit_deduped <= 0):
            raise ValueError('Roulette weights must be finite and positive')
        # ---------- return
        return fit_deduped
