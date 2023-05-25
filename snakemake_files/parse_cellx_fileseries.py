from dataclasses import dataclass
import json
import os
from pathlib import Path
from typing import List, Optional

import xmltodict
import yaml


@dataclass
class CellX:
    si: int
    sj: int
    filepaths: List[Path]


def parse_cellx_fileseries(filepath: Path) -> List[CellX]:
    """
    Parse a CellX fileseries.xml file and returns a list of tuples for all series (si) and images (sj) found as well as
    their expected relative output paths as a list.

    Args:
        filepath (`str`)   Full path + file name to the fileseries file.

    Returns:
        out (`list`):    List for each series found (si) and each frame (sj) found as well as  their expected
                            relative output paths as a list.
    """
    out = []  # list which gets returned
    series = []  # series number (si)
    series_file_set = []  # all file sets (images) for a given series number (si)
    fluo_types = []  # series sub folder name (si)
    cellx_result_dir = ""  # CellX output directory defined in fileseries
    # definition of all output files produced by CellX for a given si, sj
    expected_output_files = [
        "cells_{j:05d}.mat",
        "cells_{j:05d}.txt",
        "control_{j:05d}.png",
        "mask_{j:05d}.mat",
        "seeding_{j:05d}.png",
    ]

    try:
        with open(filepath) as fd:
            doc = xmltodict.parse(fd.read())

        # Load the fileseries
        file_series_list = doc["CellXFiles"]["CellXTimeSeries"]
        if not isinstance(file_series_list, list):
            file_series_list = [file_series_list]
        for file_series in file_series_list:
            for key, value in file_series.items():
                if key == "@id":
                    series.append(yaml.safe_load(json.loads(json.dumps(value))))
                elif key == "@fluotypes":
                    fluo_types.append(value)
                elif key == "CellXResultDir":
                    cellx_result_dir = value
                elif key == "CellXFileSet":
                    series_file_set.append(value)

                    # Extract frame numbers of current series:
                    frames_in_file_set = []

                    # One image is organized slightly different in the dict
                    if isinstance(series_file_set[0], list):
                        file_sets = series_file_set[0]
                    else:
                        file_sets = series_file_set

                    for i in range(0, len(file_sets)):
                        for key2, value2 in file_sets[i].items():
                            if key2 == "@frame":
                                frames_in_file_set.append(
                                    yaml.safe_load(json.loads(json.dumps(value2)))
                                )

                        # Build expected relative paths for final CellX output:
                        file_set_paths = []
                        for k in range(0, len(expected_output_files)):
                            file_set_paths.append(
                                Path(
                                    cellx_result_dir,
                                    fluo_types[-1],
                                    expected_output_files[k].format(
                                        j=frames_in_file_set[-1]
                                    ),
                                )
                            )

                            out.append(
                                CellX(
                                    si=series[-1],
                                    sj=frames_in_file_set[-1],
                                    filepaths=file_set_paths,
                                )
                            )

        return out

    except Exception as e:
        print(e)


def find_si_sj(res: List[CellX], filepaths: List[Path]):
    """
    Find si and sj from knowing the output filepaths.
    
    """
    idx = ([el.filepaths for el in res]).index([Path(el) for el in filepaths])
    return {"si": res[idx].si, "sj": res[idx].sj}


def find_cellxpath(cellx_from_config: Optional[str]) -> Path:
    """
    Find the path of CellX.
    
    One can use os.name or sys.platform to find more about the running OS.
    """
    if cellx_from_config is None:
        if os.getenv("CELLXROOT") is not None:
            cellxpath = Path(os.environ["CELLXROOT"], "CellX.sh")
        else:
            raise ValueError(f"'--config cellx' is missing")
    else:
        cellxpath = Path(cellx_from_config)
    return cellxpath


def test_parse_cellx_fileseries():
    assert len(parse_cellx_fileseries(Path(".") / "test" / "fileseries.xml")) == 100


def test_parse_cellx_fileseries_0():
    assert len(parse_cellx_fileseries(Path(".") / "test" / "fileseries_0.xml")) == 5
