"""Deal with data of the realisations of the birth-death process simulated in
rust with binary `ecdna-evo`.
"""

import re
import json
import pandas as pd
from pathlib import Path
import sys
from typing import Any, Dict, List, Set, Union
from futils import snapshot





def load_histogram(path: Path) -> snapshot.Histogram:
    try:
        hist = snapshot.histogram_from_file(path)
    except json.JSONDecodeError as e:
        print(f"Error in opening {path} {e}")
        sys.exit(1)
    return hist


class RealisationDistribution:
    def __init__(self, distribution: snapshot.Histogram, params: Parameters) -> None:
        self.distribution = distribution
        self.parameters = params


def realisation_distribution_from_path(path: Path) -> RealisationDistribution:
    assert path.is_file(), f"cannot find ecDNA distribution file {path}"
    return RealisationDistribution(load_histogram(path), parameters_from_path(path))


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
