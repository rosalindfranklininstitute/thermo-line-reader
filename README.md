<!--
SPDX-FileCopyrightText: 2026 Duncan McDougall <duncan.mcdougall@rfi.ac.uk>

SPDX-License-Identifier: LicenseRef-RFI-Apache-2.0-Commons-clause
-->

# Thermo Raw Reader NeXus Tools

This tool can be used to read a collection Thermo `.raw` files into a Mass
Spectrometry Image (MSI) file. Each `.raw` file represents a horizontal line of
the image ,travelling left to right. So that the number of $y$ values is equal
to the number of files, and the $x$ values are related to the time of each
scan.

This repo uses
[RawFileReader](https://github.com/thermofisherlsms/RawFileReader/) to read
mass spectrometry data. The data is converted to an HDF5 file structured with
[NeXus](https://www.nexusformat.org/). The conversion is done with the
[`Xartus`](https://github.com/rosalindfranklininstitute/Xartus) library.

## License

This tool is released under a dual license. The core code is under the Apache-2.0
license with a Commons Clause. The included dlls are provided licensed
under the "RawFileReader" license: "RawFileReader reading tool. Copyright ©
2016 by Thermo Fisher Scientific, Inc. All rights reserved." However,
`thermo_line_reader` is not an official product of Themo Fischer, nor is it
associated with them in any way, beyond using the library they provide in
"RawFileReader".

## Example

Given '/input/data/', a directory with the thermo `.raw` files for the lines,
named 'some_prefix-{line_number}.raw'. The following will write the converted
data to '/output/data.nxs'. If the output should have pixels of $200\mu m \times 200\mu m$, 
and a line speed of $20 \fram{\mu m}{s}$. You could convert
this wither using a pixel width in time or in space. As shown below.

### Using distance

The $x$ distance is specified as distance $200\mu m$. This will use the line
speed ($20 \mu m/s$) to determine how to bin the input files. It can be
done using only the command line with:

```bash
thermo_line_reader --input '/input/data/' --output '/output/data.nxs' --micron-per-line 200 --micron-per-second 20 --pixel-width 200 --pixel-metric 'distance'
```

Or if using a config file, run:

```bash
thermo_line_reader --config '/input/data/config.toml'
```

With the content of '/input/data/config.toml' as:

```toml
[thermo_line_reader]
input="/input/data/" 
output="/output/data.nxs" 
micron-per-line=200 
micron-per-second=20 
pixel-width=200
pixel-metric='distance'
```

### Using time

The $x$ distance is specified as $10s$. It can be done using only the command line with:
The $x$ distance is specified as time with $10s$. This will bin the input, and use the line
speed ($20 \mu m/s$) convert this to distance. It can be done using only the command line with:

```bash
thermo_line_reader --input '/input/data/' --output '/output/data.nxs' --micron-per-line 200 --micron-per-second 20 --pixel-width 10 --pixel-metric 'time'
```

Or if using a config file, run:

```bash
thermo_line_reader --config '/input/data/config.toml'
```

With the content of '/input/data/config.toml' as:

```toml
[thermo_line_reader]
input="/input/data/" 
output="/output/data.nxs" 
micron-per-line=200 
micron-per-second=20 
pixel-width=10
pixel-metric='time'
```


## Usage

The Thermo Line Reader can be used on the command line. The only required
inpurs are `--input` and `--output`, the rest have defaults. However specifying
the pixel diemnsions is important for representing the data correctly. In
addition, there are three extra output formats that can be written. The details
of these options are given below. All these options can also be specified in a
config file that removes the need for repeating them. In addition, this config
file can be used as a kind of metadata for a particular set of line scans.

A config file can be specified with `--config '/path/to/config.toml'`. In
addition, a simple graphical interface can be shown using `--interactive`. Note
that if cli options are specified (including the config) those values will be
pre-loaded into the GUI.

### Pixel Dimensions

The pixel dimensions are specified in $y$ as the line height
(`--micron-per-line`). For the $x$ direction the scan speed
(`--micron-per-second`) and the size of the pixels. The size of the pixels is
specified with a width (`--pixel-width`) in either seconds (with
`--pixel-metric 'time'`) or in distance (with `--pixel-metric 'distance'`/`pixel-metric='distance'`).

### Additional outputs

By default this code writes out to NeXus files using [`Xartus`](https://github.com/rosalindfranklininstitute/Xartus). However, three additional output formats are supported.

1. `--write-imzml`/`write-imzml=true`: 'processed' imzML files can be written ( '{output}.imzML' and '{output_path}.ibd'). This takes more time than the core
   conversion. For moderate datasizes ($30 \times 30\times 4000$) this should
   not be a problem. But for bigger datasets this could be disabled to speed up
   conversion.

2. `--write-unidec`/`write-unidec=true`: `unidec` files can be written. This outputs two files that can be read by `unidec`:
   1. A '{output}.total_spectrum.txt' file which contains the total spectrum of all the input data.
   2. A '{output}.unidec.hdf5' file which writes the total spectrum for each line ($y$ value in the output). 

3. `--write-tic`/`write-tic=true`: a Total Ion Count (TIC) image. This writes a TIC image to disc for the data
   ('{output}.total_ion_count.png'). Useful for checking that the pixel sizes
   are as expected. 


### List of all outputs

```bash
options:
  -h, --help            show this help message and exit
  --int, --interactive  If present will present the arguments interactively, instead of on the command line.
  -c, --config CONFIG   The path to a configuration file. This config file can be used instead of passing files into the command line. It is a TOML formatted file. The arguments will be   
                        read out of the name of the program. Default 'None'.
  -o, --output OUT_PATH
                        The directory to place the requested images and spectra. Default 'None'.
  --chunk-bytes, --chunk-max-byte-count CHUNK_MAX_BYTE_COUNT
                        The maximum number of bytes of a chunk in the output file. Default '8388608'.
  --memory-bytes, --memory-max-byte-count MEMORY_MAX_BYTE_COUNT
                        The maximum number of bytes to use as the memory buffer. Each thread uses this much memory. Default '4294967296'.
  --show-progress, --no-show-progress
                        If true (if not specified) tqdm is used to show the progress. Otherwise no progress is shown. Default 'True'.
  --dry-run             If true (if specified) no conversion will be done, but all the metadata will be read and displayed.
  --input, --in-path IN_PATH
                        The name of the folder to import data from. Default 'None'.
  --filename-prefix FILENAME_PREFIX
                        A prefix to ignore on the filenames. Default 'None'.
  --pixel-width PIXEL_WIDTH
                        The amount of time or distance per pixel. The units are determined by --width-metric Default '0.025'.
  --pixel-metric {PixelMetric.TIME,PixelMetric.DISTANCE}
                        The type metric used to define the width of a pixel. Default 'PixelMetric.TIME'.
  --time-bounds {TimeBounds.SUBSET,TimeBounds.SUPERSET}
                        What timing to use for aligning the lines. Subset: data between the largest start time and the smallest end time will be used, and the rest discarded. Superset:    
                        all data will be included from the smallest start time to the largest end time. Lines which do not fully span this will have their data padded with zeros. Default  
                        'TimeBounds.SUBSET'.
  --micron-per-second MICRON_PER_SECOND
                        The line speed of the scan. Used to convert the line times into spacial dimensions. Each line is stored along the X axis. Default '1.0'.
  --micron-per-line MICRON_PER_LINE
                        The distance between each line. Used to convert the separate files into the spacial dimensions. The lines are concatenated align the Y axis. Default '1.0'.
  --write-unidec, --no-write-unidec
                        If present will not write out the total spectra and line spectra to a format UniDec can read. Default 'False'.
  --write-tic, --no-write-tic
                        If present will not write out the tic image file. Default 'False'.
  --write-imzml, --no-write-imzml
                        If present will write out an imzML file in processed mode, removing zeros formthe spectra. Default 'False'.
  -n, --down-sampling, --downsampling, --down-sampling DOWN_SAMPLING
                        The number of 'mass resolutions' to include in each mass bucket. For instance if the mass is between 2000 and 5000 with a mass resolution of 0.5, there should be   
                        6000 mass bins. If '--down-sample=10', then there will be 600 mass bins in the output file. Default '1'.
```

## Config


All the settings can be specified in a config file. If a setting is not
specified, it will take its default value. In addition, values specified on the
command line are facoured over those specified in the config file. 

Particular parameters to note are the `write-` parameters. These use TOML
`true` and `false` not the python syntax. In addition, the `pixel-metric`
which should be either 'time' or 'distance'.

```toml
[thermo_line_reader]
input="path/to/directory/with/.raw/"
output="path/to/output.nxs"
micron-per-second=20
micron-per-line=200
pixel-width=10
pixel-metric='time'
show-progress=true
write-unidec=false
write-tic=true
write-imzml=false
down-sampling=1
```

