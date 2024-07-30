from scipy import stats
from futils import snapshot
from typing import Dict, List, Union
from .realisation import RealisationDistribution
from futils import abc


def summary_statistics(
    sims: List[RealisationDistribution],
    target: snapshot.Histogram,
    target_name: str,
    stats: abc.Stats,
    ) -> List[Dict[str, Union[str, float, int]]]:
    """
    Return a list of records (list of dict) with the summary statistics and
    other quantities such as the parameters used.
    We compute the ecDNA distributions from the simulations `sims` and compare
    that against the the distributions from the data.

    """
    wasserstein = abc.Wasserstein()
    all_params = []

    for i, my_ecdna in enumerate(sims):
        params = dict(**my_ecdna.parameters)
        params["patient"] = target_name
        # compute the summary statistics
        for s in stats:
            params[s.__class__.__name__] = s.distance(target, my_ecdna.distribution)
        all_params.append(params)

    return all_params
