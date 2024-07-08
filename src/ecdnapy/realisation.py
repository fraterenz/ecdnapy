"""Deal with data of the realisations of the birth-death process simulated in
rust with binary `ecdna-evo`.
"""
import numpy as np
import re
import json
import pandas as pd
from pathlib import Path
import sys
from scipy import stats
from typing import Any, Dict, List, Set, Union
from futils import snapshot, parsing


def load_histogram(path: Path) -> snapshot.Histogram:
    try:
        hist = snapshot.histogram_from_file(path)
    except json.JSONDecodeError as e:
        print(f"Error in opening {path} {e}")
        sys.exit(1)
    return hist


class RealisationDistribution:
    def __init__(self, distribution: snapshot.Histogram, params: parsing.Parameters) -> None:
        self.distribution = distribution
        self.parameters = params
        # cache
        self.distribution_array: Union[np.ndarray, None] = None

    def mean(self) -> float:
        if self.distribution_array is None:
            self.distribution_array = snapshot.array_from_hist(self.distribution)
        return self.distribution_array.mean()

    def var(self) -> float:
        if self.distribution_array is None:
            self.distribution_array = snapshot.array_from_hist(self.distribution)
        return self.distribution_array.var()

    def entropy(self) -> float:
        if self.distribution_array is None:
            self.distribution_array = snapshot.array_from_hist(self.distribution)
        return float(stats.entropy(self.distribution_array))


def realisation_distribution_from_path(path: Path) -> RealisationDistribution:
    assert path.is_file(), f"cannot find ecDNA distribution file {path}"
    return RealisationDistribution(load_histogram(path), parsing.parameters_from_path(path))


def load_ecdnas_from_folder(
    path2dir: Path, max2load: Union[int, None] = None
) -> List[RealisationDistribution]:
    assert path2dir.is_dir()
    realisations = []

    for path in path2dir.iterdir():
        i = 0
        for p in path.glob("*.json"):
            if max2load and i >= max2load:
                break
            realisations.append(realisation_distribution_from_path(p))
    print(f"loaded {len(realisations)} files from {path2dir}")
    return realisations
