#!/usr/bin/env python3
import sys
import argparse
import os
import datetime as dt

import mapping as map_funks
import report_writer as r_writer


def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(description='Run mutation report')

    parser.add_argument('--date-data', help="When the report was run", dest="date_data")
    parser.add_argument("--metadata-file", dest="metadata_file")
    parser.add_argument("--snp-file", dest="snp_file")
    parser.add_argument("--snp-list",help="list of snps desired in the report", required=True,dest="snp_list")
    parser.add_argument("--figdir")
    parser.add_argument("--snps-for-matrix", dest="snps_for_matrix", default=None)
    parser.add_argument("--uk-map", dest="uk_map")
    parser.add_argument("--ni-counties", dest="ni_counties")
    parser.add_argument("--channel-islands", dest="channel_islands")
    parser.add_argument("--outdir")
    parser.add_argument("--date-range-start", dest = "date_start")
    parser.add_argument("--date-range-end", dest="date_end")


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

    uk = args.uk_map
    ni = args.ni_counties
    channels = args.channel_islands

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if outdir not in figdir:
        figdir_writing = f'{outdir}/{figdir}'
    else:
        figdir_writing = figdir
        figdir = figdir_writing.replace(outdir,".")

    if not os.path.exists(figdir_writing):
        os.mkdir(figdir_writing)

    map_files = [uk,channels,ni]

    all_uk = map_funks.generate_all_uk_dataframe(map_files)

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

    
    r_writer.generate_report(metadata_file, date_data, snp_file, snp_list, date_start, date_end, snps_for_matrix, figdir_writing, figdir, outdir, all_uk)

if __name__ == '__main__':
    main()