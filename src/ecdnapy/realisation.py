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


class Parameters:
    def __init__(
        self,
        path: Path,
        sample: int,
        population: int,
        b0: float,
        b1: float,
        d0: float,
        d1: float,
        idx: int,
    ):
        self.sample = sample
        self.path = path
        self.cells = population
        self.b0 = b0
        self.b1 = b1
        self.d0 = d0
        self.d1 = d1
        self.idx = idx

    def into_dict(self) -> Dict[str, Any]:
        return self.__dict__

    def stringify(self, some_params: Set[str]) -> str:
        return ", ".join(
            [f"{k}={v}" for k, v in self.into_dict().items() if k in some_params]
        )


class ParametersFile:
    def __init__(
        self,
        b0: float,
        b1: float,
        d0: float,
        d1: float,
        idx: int,
    ):
        self.b0 = b0
        self.b1 = b1
        self.d0 = d0
        self.d1 = d1
        self.idx = int(idx)

    def into_dict(self) -> Dict[str, Any]:
        return self.__dict__


def parameters_from_path(path: Path) -> Parameters:
    """Assume something like
    100000samples100000population/ecdna/1dot1b0_1b1_0d0_0d1_0idx.json
    """
    parts = path.parts
    match_sample = re.compile(r"^(\d+)(samples)(\d+)(population)$", re.IGNORECASE)
    sample, cells = 0, 0
    for part in parts:
        matched = match_sample.search(part)
        if matched:
            sample = int(matched.group(1))
            cells = int(matched.group(3))

    assert sample > 0, f"cannot find a value for samples from {path}"
    assert cells > 0, f"cannot find a value for cells from {path}"

    params_file = parse_filename_into_parameters(path)
    return Parameters(path, sample, cells, **params_file.__dict__)


def parse_filename_into_parameters(filename: Path) -> ParametersFile:
    match_nb = re.compile(r"(\d+\.?\d*)([a-z]+[0-1]?)", re.IGNORECASE)
    filename_str = filename.stem
    filename_str = filename_str.replace("dot", ".").split("_")
    my_dict = dict()

    for ele in filename_str:
        matched = match_nb.search(ele)
        if matched:
            my_dict[matched.group(2)] = float(matched.group(1))
        else:
            raise ValueError(f"could not parse the filename into parameters {filename}")

    return ParametersFile(**my_dict)


def params_into_dataframe(params: List[Parameters]) -> pd.DataFrame:
    df = pd.DataFrame.from_records([param.into_dict() for param in params])
    df.idx = df.idx.astype(int)
    df.cells = df.cells.astype(int)
    df["samples"] = df["samples"].astype(int)
    return df


def params_files_into_dataframe(params: List[ParametersFile]) -> pd.DataFrame:
    df = pd.DataFrame.from_records([param.into_dict() for param in params])
    df.idx = df.idx.astype(int)
    df.cells = df.cells.astype(int)
    return df


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
        if max2load and i >= max2load:
            break
        for p in path.glob("*.json"):
            realisations.append(realisation_distribution_from_path(p))
    print(f"loaded {len(realisations)} files from {path2dir}")
    return realisations
