# mink

**M**utation **I**dentification **N** **K**inetics (tbc)

**mink** is a report generation tool for tracking mutations of interest through time and space in the UK. 

It was created by Verity Hill, and uses code from Áine O'Toole and Ben Jackson, all in Rambaut Group at the University of Edinburgh.

NB mink is non-phylogenetic and relies on SNPs already having been called by other software (for example COG-UK's phylopipe). 

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
--snp-csv SNP-CSV
                CSV file as alternative to SNP-LIST containing groups of SNPs of interest. 
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
--find-fastest
                Flag to find the ten fastest growing SNPs in the last month.
```


Format of SNP input:

To use the snp-csv input, the csv should contain three columns headed: group, snps, and description.
The group column allows mink to analyse the snps in groups, and the group name is used to structure the report.
The snps contains a pipe-separated list of amino acid changes.
The description is a free text line which will  be added to the appropriate section of the report.

The snps, whether given as a list on the command line, or in the right column of the csv must be in the following format:

gene:{referenceaminoacid}POSITION{newaminoacid}

eg to look for the D614G mutation in the spike protein, the input will be S:D614G

To look for all possible mutations from the reference in the same gene, input S:D614
mink will find all of the mutations at that position in the input file.

Stop codons are represented by a "*",eg S:D614\*

mink does not currently find deletions or insertions. This will be added at a later date.