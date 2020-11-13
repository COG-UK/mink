#!/usr/bin/env python3
import sys
import argparse
import os
import datetime as dt
import pkg_resources

import mapping as map_funks
import report_writer as r_writer


def main(sysargs = sys.argv[1:]):

    cwd = os.getcwd()

    parser = argparse.ArgumentParser(description='Run mutation report')

    parser.add_argument("--metadata-file", dest="metadata_file", help="path to metadata file with date and time information for sequences.")
    parser.add_argument("--snp-file", dest="snp_file", help="path to csv file containing which snps are in which sequences")
    parser.add_argument("--snp-list",help="list of snps desired in the report", required=True,dest="snp_list")
    parser.add_argument("--snps-for-matrix", dest="snps_for_matrix", default=None, help="list of snps for co-occurence matrix. Default is all snps in snp-list")
    parser.add_argument('--date-data', help="When the report was run", dest="date_data")
    parser.add_argument("--figdir", default="figures", help="figure directory name")
    parser.add_argument("--outdir", default="mink_results", help="output directory name")
    parser.add_argument("--date-start", dest = "date_start", help="restrict analysis to this date at the earliest")
    parser.add_argument("--date-end", dest="date_end", help="restrict analysis to this date at the latest")

    """
    Exit with help menu if no args supplied
    """
    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(0)
    else:
        args = parser.parse_args(sysargs)
        if args.help:
            parser.print_help()
            sys.exit(0)

    args = parser.parse_args()

    date_data = args.date_data
    metadata_file = args.metadata_file
    snp_file = args.snp_file
    snp_list = args.snp_list
    figdir = args.figdir
    snps_for_matrix = args.snps_for_matrix
    outdir = args.outdir
    date_start = args.date_start
    date_end = args.date_end

    if not os.path.exists(outdir):
        os.mkdir(os.path.join(cwd,outdir))

    if outdir not in figdir:
        figdir_writing = os.path.join(outdir,figdir)
    else:
        figdir_writing = figdir
        figdir = figdir_writing.replace(outdir,".")

    if not os.path.exists(figdir_writing):
        os.mkdir(figdir_writing)

    if type(snp_list) == str:
        snp_list = snp_list.split(",")
    if type(snps_for_matrix) == str:
        snps_for_matrix = snps_for_matrix.split(",")

    if not snps_for_matrix:
        snps_for_matrix = snp_list

    if args.date_start:
        date_start = dt.datetime.strptime(args.date_start,"%Y-%m-%d").date()
    else:
        date_start = dt.date(2020,1,1)

    if args.date_end:
        date_end = dt.datetime.strptime(args.date_end,"%Y-%m-%d").date()
    else:
        date_end = dt.date.today()


    uk = pkg_resources.resource_filename('mink', 'data/mapping_files/gadm36_GBR_2.json')
    ni = pkg_resources.resource_filename('mink', 'data/mapping_files/NI_counties.geojson')
    channels = pkg_resources.resource_filename('mink', 'data/mapping_files/channel_islands.json')

    map_files = [uk,channels,ni]

    all_uk = map_funks.generate_all_uk_dataframe(map_files)

    
    r_writer.generate_report(metadata_file, date_data, snp_file, snp_list, date_start, date_end, snps_for_matrix, figdir_writing, figdir, outdir, all_uk)

if __name__ == '__main__':
    main()