# SPDX-FileCopyrightText: 2026 Duncan McDougall <duncan.mcdougall@rfi.ac.uk>
#
# SPDX-License-Identifier: LicenseRef-RFI-Apache-2.0-Commons-clause
from pathlib import Path
import logging
import itertools

import pytest

from thermo_line_reader import (
    ProcessArgs,
    process,
    PixelMetric,
    TimeBounds,
    parse_args,
    data_convert,
)
import xarray as xr


def default_options(out_path) -> ProcessArgs:
    process_args = ProcessArgs()

    process_args.in_path = (
        Path(__file__).parent / "data" / "20260904_MouseTesticle__Y200um_20ums/"
    )
    process_args.out_path = out_path / "data.nxs"

    process_args.pixel_width = 1.3
    process_args.pixel_metric = PixelMetric.TIME

    process_args.time_bounds = TimeBounds.SUBSET

    process_args.micron_per_second = 20
    process_args.micron_per_line = 200

    process_args.write_unidec = True
    process_args.write_imzml = True
    process_args.write_tic = True

    process_args.show_progress = False

    return process_args


@pytest.mark.parametrize(
    ("unidec", "imzml", "tic"),
    [
        (unidec, imzml, tic)
        for unidec, imzml, tic in itertools.product(
            [True, False], [True, False], [True, False]
        )
    ],
)
def test_boolean_options(tmp_path, unidec, imzml, tic):

    logging.getLogger(data_convert.__name__).setLevel(logging.INFO)

    process_args = default_options(tmp_path)

    process_args.write_unidec = unidec
    process_args.write_imzml = imzml
    process_args.write_tic = tic

    process_args = parse_args(process_args)

    process(process_args, {})

    assert (tmp_path / "data.nxs").exists()
    assert (tmp_path / "data.imzML").exists() == imzml
    assert (tmp_path / "data.ibd").exists() == imzml
    assert (tmp_path / "data.total_ion_count.png").exists() == tic
    assert (tmp_path / "data.total_spectrum.txt").exists() == unidec
    assert (tmp_path / "data.unidec.hdf5").exists() == unidec

    array = xr.open_dataarray(
        tmp_path / "data.nxs", engine="nexus", entry_path="/entry/data/images"
    )
    assert array.shape == (207, 29, 4001)
    assert "x" in array.coords
    assert "y" in array.coords
    assert "time" in array.coords
    assert "mz" in array.coords
