from em_7_score import em_score
from numpy import max, sum


def em_run(self):
    # ML and scoring

    def em_max_sum(self):
        # These index the best results below
        self.max_likely_results["Max Final Score"] = max(
            self.max_likely_results["Final Scores"])
        self.max_likely_results["Max Sum Scores"] = sum(
            self.max_likely_results["Final Scores"])

    def em_best_result(self):
        # Last fxn within the em_core for loop
        self.max_likely_results[
            "Max Final Sequence"] = self.max_likely_results[
                "Final Scores"].index(
                    self.max_likely_results["Max Final Score"])
        self.max_likely_results[
            "Max Final Position"] = self.max_likely_results["Final Positions"][
                self.max_likely_results["Max Final Sequence"]]
        self.max_likely_results["Max Final Motif"] = self.max_likely_results[
            "Final Motifs"][self.max_likely_results["Max Final Sequence"]]

    em_score(self)
    em_max_sum(self)
    em_best_result(self)
