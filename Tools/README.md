# Tools for generating imgCIF files

## `imgcif_mapper.py`

## `cbf_scan_extractor`

`cbf_scan_extractor` analyses a sequence of miniCBF files contained in 
the indicated directory and makes a best-effort attempt to divide them
into scans and to extract the axis settings for each scan. This has not
been widely tested on the many possible miniCBF formats and file naming
schemes, so please check the final output.

The final output is a _fragment_ of a complete imgCIF file. To produce a
complete file from this fragment, the following steps are necessary:

1. the axis and detector geometric information (such as those
that have been collected together in this repository) should be 
prepended.

2. Actual values for any missing geometric items (`?` characters in
the above file) should be inserted. Currently this is the beam centre
in mm when all axes are at their zero settings and the dimensions of
the detector in pixels.

3. A data block header of the form `data_<someheadername>` prepended.

## Limitations

1. All frames must be in a single local directory
2. Frame filenames must be of the form `<xxxxx>_ss(0|_)ffff.cbf`
where `<xxxxx>` is arbitrary, `ss` is an arbitrary-length integer
labelling the scan, `ffff` is an arbitrary-length integer
sequentially numbering the frames within the scan starting from
1, and the scan label is followed either by the digit `0` or
an underscore.
3. Many variants of axis names may not be recognised. Please
report any that you think should be recognised.

### Installation

1. [Install Julia](https://julialang.org/downloads)

2. Run Julia. At the prompt, type `]` and then `add ArgParse`. This
will download and install the `ArgParse` Julia package.

3. Press `backspace` and then exit Julia (`exit`)

### Running

Type `julia cbf_scan_extractor.jl --help` at a command prompt to get
a list of options.

#### Examples

```julia cbf_scan_extractor.jl /home/me/CBF_crystal_2```

Analyse frames in directory `/home/me/CBF_crystal_2` and output an imgCIF fragment
to the terminal, using a local URL to specify the frame locations.

```julia cbf_scan_extractor.jl -o stoe_test.cif -a Detector_distance trans -a detector_2theta two_theta -l https://zenodo.org/record/123456/files/cbf_crystal2.tar.bz2 -i /home/me/CBF_crystal_2```

Analyse frames in directory `/home/me/CBF_crystal_2` and then output an imgCIF fragment to
`stoe_test.cif` assuming that the actual network location of an archive file containing
the frames is 
`https://zenodo.org/record/123456/files/cbfcrystal2.tar.bz2`, prepending the directory
name (`-i` option) to the location specification within the archive. Axes named
`detector_2theta` and `Detector_distance` in the frames are renamed to `two_theta`
and `trans` in the imgCIF file.

