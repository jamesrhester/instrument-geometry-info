
# WARNING: UNDER CONSTRUCTION
# Instrument layouts

## Introduction

This repository collects information on the geometrical layout of specific X-ray and neutron diffraction instruments and
beamlines. Such information is essential for reducing raw diffraction data. Information in the raw data files,
when available, should be used in preference to the information here. Use this information as a default or backup option 
where the geometrical information is not otherwise available. 

## Contributing

The following options are available for contributing a new layout:

1. Raise an issue using the appropriate issue template: click on 'Issues' above, then 'New Issue' and submit a new beamline layout.
2. Create a pull request on the repository. See `Beamlines/beamline_geometries.cif` for the way the information should be presented.

Note that by contributing you are implicitly licensing your contribution as CC0 (public domain).
  
## Format of information

The geometric information is stored using imgCIF data names. Each beamline configuration is contained in a separate data block
of files `Beamlines/beamline_geometries.cif` (for beamlines) or `Laboratory/laboratory_geometries.cif` (for off-the-shelf 
laboratory instruments).

## Accuracy

While we try our best, unfortunately we can make no claims as to the accuracy of this information as it has often been provided
by third parties. Treat this information as default or backup information, always preferring information in the raw data files
themselves. Mistakes in geometry specifications will usually lead to problems in data reduction, so do check that data
are being correctly reduced, and please advise of any errors by raising an issue.

## Tools

The `Tools` directory contains some useful programs for generating imgCIF files. See the `README` in
that directory for further information.
