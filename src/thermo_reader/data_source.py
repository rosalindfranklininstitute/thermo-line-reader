# SPDX-FileCopyrightText: 2026 RFI
#
# SPDX-License-Identifier: Apache-2.0
import datetime as dt
import logging
import math
import os
import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Callable, NamedTuple, override

import numpy as np
from xartus.lib import (
    AbstractDataSource,
    Axis,
    Chunk,
    DataShape,
    MultiCOO,
    MultiLinspace,
    Shape,
    Signal,
)
from xartus.lib.dtypes import Float1D32, Int1D32, Int3D32
from xartus.lib.writers import AbstractWriter, SumAndMax

from .load_thermo import (
    RawFileReaderAdapter,
)
from .msi_instrument import MSInstrumentData

logger = logging.getLogger(__name__)


class PixelMetric(Enum):
    TIME = "time"
    DISTANCE = "distance"


class TimeBounds(Enum):
    SUBSET = "subset"
    SUPERSET = "superset"


class RawLine(NamedTuple):
    file: Path
    number: int
    ignore: bool


@dataclass
class MinMax:
    min: float
    max: float


class DataLines(NamedTuple):
    lines: list[RawLine]
    scan_times: list[list[float]]
    scan_dates: list[dt.datetime]
    mass: MinMax
    times: MinMax
    mass_resolution: float

    @staticmethod
    def from_lines(lines: list[RawLine], time_bounds: TimeBounds) -> "DataLines":
        first_number = lines[0].number
        all_creation_times = []
        all_times = []
        mm_mass: None | MinMax = None
        mm_time: None | MinMax = None
        mass_resolution: None | float = None

        for ii, line in enumerate(lines):
            lines[ii] = RawLine(line.file, line.number - first_number, False)

            raw_file = RawFileReaderAdapter.FileFactory(str(lines[ii].file))
            if not raw_file.IsOpen or raw_file.IsError:
                raise RuntimeError(
                    "Unable to access the RAW file using the RawFileReader class!"
                )
            if raw_file.IsError:
                raise IOError(
                    "Error opening ({}) - {}".format(raw_file.FileError, lines[ii].file)
                )
            if raw_file.InAcquisition:
                raise IOError(
                    "RAW file still being acquired - {}".format(lines[ii].file)
                )

            instrument_data = MSInstrumentData(raw_file, 1)
            line_times = instrument_data.get_scan_times()
            if len(line_times) == 0:
                logger.warning(f"{lines[ii].file.name}: -- No data found, ignoring")
                lines[ii] = RawLine(lines[ii].file, lines[ii].number, True)
                continue

            tmp_min_mz, tmp_max_mz = instrument_data.get_mass_range()
            mass_resolution = instrument_data.stored_mass_resolution()
            tmp_min_s = np.min(line_times)
            tmp_max_s = np.max(line_times)
            logger.debug(
                f"{lines[ii].file.name}: ({tmp_min_s / 60: >.4f} - {tmp_max_s / 60: >.4f}min) {len(line_times)} scans"
            )

            if mm_mass is None or mm_time is None or mass_resolution is None:
                mm_mass = MinMax(tmp_min_mz, tmp_max_mz)
                mm_time = MinMax(np.min(line_times), np.max(line_times))
                mass_resolution = instrument_data.stored_mass_resolution()
            else:
                match time_bounds:
                    case TimeBounds.SUBSET:
                        mm_mass.min = max(mm_mass.min, tmp_min_mz)
                        mm_mass.max = min(mm_mass.max, tmp_max_mz)
                        mm_time.min = max(mm_time.min, tmp_min_s)
                        mm_time.max = min(mm_time.max, tmp_max_s)
                    case TimeBounds.SUPERSET:
                        mm_mass.min = min(mm_mass.min, tmp_min_mz)
                        mm_mass.max = max(mm_mass.max, tmp_max_mz)
                        mm_time.min = min(mm_time.min, tmp_min_s)
                        mm_time.max = max(mm_time.max, tmp_max_s)
                mass_resolution = min(
                    mass_resolution, instrument_data.stored_mass_resolution()
                )
            all_times.append(line_times)
            all_creation_times.append(instrument_data.creatoin_time)

        lines = [ll for ll in lines if not ll.ignore]

        if mm_mass is None or mm_time is None or mass_resolution is None:
            raise ValueError(f"No data found in {time_bounds.value} of lines.")

        return DataLines(
            lines=lines,
            scan_times=all_times,
            scan_dates=all_creation_times,
            mass=mm_mass,
            times=mm_time,
            mass_resolution=mass_resolution,
        )

    def get_line_data(
        self, line_index: int, instrument_index
    ) -> tuple[MSInstrumentData, list[float]]:

        raw_file = RawFileReaderAdapter.FileFactory(str(self.lines[line_index].file))
        return MSInstrumentData(raw_file, instrument_index), self.scan_times[line_index]


class ThermoDataSource(AbstractDataSource):
    def __init__(
        self,
        in_path: Path,
        filename_prefix: str | None,
        time_bounds: TimeBounds,
        pixel_metric: PixelMetric,
        pixel_width: float,
        micron_per_second: float,
        micron_per_line: float,
        sampling: MultiLinspace | None = None,
    ):
        sampling = sampling if sampling is not None else MultiLinspace()

        filename_prefix = (
            filename_prefix
            if filename_prefix is not None
            else self._find_prefix(in_path)
        )

        self.data_lines = DataLines.from_lines(
            self._read_lines_raw(in_path, filename_prefix), time_bounds
        )

        mass_count = math.ceil(
            (self.data_lines.mass.max - self.data_lines.mass.min)
            / self.data_lines.mass_resolution
        )
        self.mz_edges = sampling.get_edges(
            self.data_lines.mass.min, self.data_lines.mass.max, mass_count
        )
        match pixel_metric:
            case PixelMetric.TIME:
                delta_t = pixel_width
                delta_m = pixel_width * micron_per_second
            case PixelMetric.DISTANCE:
                delta_t = pixel_width / micron_per_second
                delta_m = pixel_width

        time_count = (
            math.ceil((self.data_lines.times.max - self.data_lines.times.min) / delta_t)
            + 1
        )
        self.time_values = np.array(
            [ii * delta_t + self.data_lines.times.min for ii in range(time_count)]
        )
        self.x_values = np.array([ii * delta_m for ii in range(time_count)])
        self.y_values = np.array(
            [ii * micron_per_line for ii in range(len(self.data_lines.lines))]
        )
        logger.info(" Lines:")
        logger.info(f"      : Count: {len(self.data_lines.lines)}")
        logger.info(f"      : Height: {micron_per_line}um")
        line_time = np.diff(np.percentile(self.time_values, [0, 100]))[0]
        logger.info(f"      : Time: {line_time / 60:.2f}min ({line_time:.1f}s)")
        line_scan_count = MinMax(
            min(len(scans) for scans in self.data_lines.scan_times),
            max(len(scans) for scans in self.data_lines.scan_times),
        )
        logger.info(
            f"      : Scans: {line_scan_count.min} - {line_scan_count.max} scans per line."
        )
        logger.info(" Pixels:")
        logger.info(f"       : Speed: {micron_per_second}um/s")
        logger.info(f"       : Width: {delta_t}s and {delta_m}um")
        logger.info("")

        logger.info(f" Pixel size: {delta_m:.1f}x{micron_per_line:.1f}um.")
        logger.info(
            f" Image size: {np.diff(np.percentile(self.x_values, [0, 100]))[0]}x{np.diff(np.percentile(self.y_values, [0, 100]))[0]}um."
        )

        self.total_shape = (
            len(self.x_values),
            len(self.y_values),
            len(self.mz_edges),
        )
        logger.info(f" Image shape: {self.total_shape}.")

        start = min(self.data_lines.scan_dates)
        end = max(self.data_lines.scan_dates) + dt.timedelta(
            seconds=self.data_lines.times.max
        )

        self.experiment_data = dict(
            pixel_width=delta_m,
            pixel_time=delta_t,
            pixel_height=micron_per_line,
            creation_date=str(start),
            end_date=str(end),
            duration=str(end - start),
            files=[str(line.file.name) for line in self.data_lines.lines],
            file_numbers=[line.number for line in self.data_lines.lines],
        )

    @override
    def __enter__(self):
        pass

    @override
    def __exit__(self, exc_type, exc_value, traceback):
        pass

    def _find_prefix(self, in_path: Path) -> str:
        paths = []
        for fle in in_path.iterdir():
            if not fle.is_file():
                continue
            if fle.suffix.lower() != ".raw":
                continue
            paths.append(fle.name)
        return os.path.commonprefix(paths)

    def _read_lines_raw(self, in_path: Path, filename_prefix: str) -> list[RawLine]:
        digits = re.compile("\\d+")
        error = False
        lines: list[RawLine] = []
        for fle in in_path.iterdir():
            if not fle.is_file():
                continue
            if fle.suffix.lower() != ".raw":
                continue

            parts = [
                int(num)
                for num in digits.findall(fle.name.removeprefix(filename_prefix))
            ]
            if len(parts) != 1:
                if len(parts) == 0:
                    logger.warning(
                        f"Could not parse filename {fle.name}: no line number included."
                    )
                else:
                    logger.warning(
                        f"Could not parse filename {fle.name}: Multiple numbers in name."
                    )
                error = True
            lines.append(RawLine(fle, parts[0], False))

        if error:
            raise RuntimeError("Could not parse all filenames")

        lines.sort(key=lambda x: x.number)
        array = np.array([line.number for line in lines])
        diffs = np.diff(array)
        for ii, diff in enumerate(diffs[1:]):
            if diff != diffs[0]:
                logger.warning(
                    f"The numbers {lines[ii].number} and {lines[ii + 1].number} are not consistently {diffs[0]} apart."
                )
        return lines

    @override
    def instrument_metadata(self) -> dict[str, Any]:
        # TODO @DMD: what can we do here?
        # https://github.com/orgs/rosalindfranklininstitute/projects/19/views/1?pane=issue&itemId=213504413
        return {}

    @override
    def user_metadata(self) -> dict[str, Any]:
        return {}

    @override
    def experiment_metadata(self) -> dict[str, Any]:
        return self.experiment_data

    @override
    def shape(self) -> DataShape:
        return DataShape(self.total_shape, worst_case_density=1.0)

    @override
    def signal_definition(self) -> Signal:
        return Signal("counts", np.int32, units="counts")

    @override
    def output_chunk_priorities(self) -> dict[str, Shape]:
        return dict(images=(1, 1, 2), spectra=(2, 2, 1))

    @override
    def chunk_read_count(self, memory_chunk: Shape) -> int:
        return int(np.prod(memory_chunk[0:2]))

    @override
    def axis_definitions(self) -> list[Axis]:
        return [
            Axis(
                name="x",
                primary_axis=0,
                units="um",
                dtype=np.float64,
            ),
            Axis(
                name="time",
                primary_axis=0,
                units="s",
                dtype=np.float64,
            ),
            Axis(
                name="y",
                primary_axis=1,
                units="um",
                dtype=np.float64,
            ),
            Axis(
                name="mz",
                primary_axis=2,
                units="mz",
                dtype=np.float64,
            ),
        ]

    @override
    def axis_values(self, axis: Axis) -> np.ndarray:
        """
        Returns the values for the specified continuous axis.
        """
        match axis.name:
            case "time":
                return self.time_values
            case "x":
                return self.x_values
            case "y":
                return self.y_values
            case "mz":
                return self.mz_edges
            case _:
                raise ValueError(f"Unknown continuous axis requested: {axis.name}")

    @override
    def output_accumulations(self) -> list[AbstractWriter]:
        return [
            SumAndMax("total_image", axes=["mz"]),
            SumAndMax("total_spectra", axes=["x", "y"]),
        ]

    @override
    def fill_chunk(
        self,
        memory_chunk: Chunk,
        update: Callable[[int], None],
    ) -> np.ndarray | MultiCOO:
        min_time = self.time_values[memory_chunk[0].start]
        max_time = self.time_values[memory_chunk[0].stop - 1]

        coords: list[Int3D32] = []
        data: list[Int1D32] = []
        mz_data: list[Float1D32] = []

        for yy in memory_chunk.range(1):
            line_data, line_times = self.data_lines.get_line_data(yy, 1)
            line_chunk_indices = np.array(
                [inx for inx, t in enumerate(line_times) if min_time <= t <= max_time]
            )
            time_axis_indices = np.searchsorted(
                self.time_values, line_times[line_chunk_indices], side="left"
            )
            for inx, scan_inx in enumerate(line_chunk_indices):
                tt = time_axis_indices[inx]
                scan_number = line_data.get_scan_number(scan_inx)
                mass, spec = line_data.get_segmented_scan(scan_number)

                count = len(spec)
                scan_coords = np.tile(
                    np.array([tt, yy, 0]).reshape(3, 1),
                    (1, count),
                )
                coords.append(scan_coords)
                data.append(spec)
                mz_data.append(mass)
            update(memory_chunk.shape[0])

        axis = np.concatenate(mz_data)
        labels = np.searchsorted(self.mz_edges[1:], axis)
        labels[labels == self.total_shape[-1]] = self.total_shape[-1] - 1
        final_coords = np.concatenate(coords, axis=1)
        final_coords[2, :] = labels
        return MultiCOO(
            coords=final_coords, values={"signal": np.concatenate(data), "mz": axis}
        )
