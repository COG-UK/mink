# mink

**M**utation **I**ntegrated **N**umerical **K**inetics (tbc)

**mink** is a report generation tool for tracking mutations of interest through time and space in the UK. 

It was created by Verity Hill, and uses code from Áine O'Toole and Ben Jackson, all in Rambaut Group at the University of Edinburgh.

## Installation and requirements

mink is compatible with Linux, MacOS and Windows subsystem for linux. 
Python version 3.6 or higher is required to run mink.

To install:
1. ``git clone https://github.com/COG-UK/mink.git`` and ``cd mink``
2. ``python setup.py install``

There is currently no environment required to run mink.

### Check the install worked

Type:

```
mink
```
and you should see the help menu of mink printed

### Update mink

1. ``git pull``
2. ``python setup.py install``


## Running mink

usage: 
`mink [options]`

options:

``` 
--metadata-file METADATA
                CSV file containing name, date and location information for each sequence
--snp-file SNP-FILE
                CSV file containing sequence name and all snps found in that sequence
--snp-list SNP-LIST
                Comma separated string containing SNPs of interest that the report will be run on
--snps-for-matrix SNPs
                Comma separated string containing SNPs for pairwise co-occurrence matrix figure. Default is the full SNP list.
--date-data DATE
                Date that report was run, printed at the top of the report.
--figdir FIGDIR
                Path to where figures will be made. Can be relative or absolute. Default is "figures" in the output directory.
--outdir OUTDIR
                Path to where output files will be produced. Default is "mink_results" in the current working directory.
--date-start DATE
                Earliest date for inclusion in report, in format YYYY-MM-DD. Default is 2020-1-1.
--date-end DATE
                Latest date for includsion in report, in format YYYY-MM-DD. Default is one week before the current date.
```


