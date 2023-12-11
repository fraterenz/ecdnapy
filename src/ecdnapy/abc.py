from scipy import stats
from futils import snapshot
from typing import List
from .realisation import RealisationDistribution


def summary_statistic_wasserstein(
    sims: List[RealisationDistribution],
    target: snapshot.Histogram,
    target_name: str,
):
    """
    Return a list of records (list of dict) with the summary statistics and
    other quantities such as the parameters used.
    We compute the ecDNA distributions from the simulations `sims` and compare
    that against the the distributions from the data.

    """
    all_params = []

    for i, my_ecdna in enumerate(sims):
        # uniformise such that they have the same support which is required by
        # the wasserstein metric
        target_uniformised, sim_uniformised = snapshot.Uniformise.uniformise_histograms(
            [target, my_ecdna.distribution]
        ).make_histograms()

        assert len(target_uniformised) == len(sim_uniformised)

        v_values, v_weights = list(sim_uniformised.keys()), list(
            sim_uniformised.values()
        )
        u_values, u_weights = list(target_uniformised.keys()), list(
            target_uniformised.values()
        )

        params = my_ecdna.parameters.into_dict()
        # compute the summary statistic
        params["wasserstein"] = stats.wasserstein_distance(
            u_values, v_values, u_weights, v_weights
        )
        params["patient"] = target_name
        all_params.append(params)

    return all_params
