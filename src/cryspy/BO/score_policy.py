import physbo


class ScoreRecordingPolicy(physbo.search.discrete.Policy):
    def __init__(self, test_X, initial_data):
        super().__init__(test_X=test_X, initial_data=initial_data)
        self.selected_scores = []
        self.first_scores = None

    def _find_champion(self, f):
        champion, local_champion, local_index = super()._find_champion(f)

        # ---------- record all candidate scores for the first selection
        if self.first_scores is None:
            self.first_scores = (self.actions.copy(), f.copy())

        # ---------- record the score used for selection
        if champion != -1:
            self.selected_scores.append(
                (int(champion), float(f[local_index]))
            )

        return champion, local_champion, local_index
