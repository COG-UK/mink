#!/usr/bin/env python3
import sys
import argparse
import os
import datetime as dt
import pkg_resources
import csv
from collections import defaultdict

import mapping as map_funks
import report_writer as r_writer
import mutation_funcs as mfunks

import logging
logging.getLogger("imported_module").setLevel(logging.ERROR)


def main(sysargs = sys.argv[1:]):

    cwd = os.getcwd()

    parser = argparse.ArgumentParser(add_help=False, description='Run mutation report')

    parser.add_argument("--metadata-file", dest="metadata_file", help="path to metadata file with date and time information for sequences.")
    parser.add_argument("--snp-file", dest="snp_file", help="path to csv file containing which snps are in which sequences")
    parser.add_argument("--snp-list",help="list of snps desired in the report", dest="snp_list")
    parser.add_argument("--snp-csv", help="list of snps in a csv desired in the report", dest="snp_csv")
    parser.add_argument("--snps-for-matrix", dest="snps_for_matrix", default=None, help="list of snps for co-occurence matrix. Default is all snps in snp-list")
    parser.add_argument('--date-data', help="When the report was run", dest="date_data")
    parser.add_argument("--figdir", default="figures", help="figure directory name")
    parser.add_argument("--outdir", default="mink_results", help="output directory name")
    parser.add_argument("--raw-data-dir", default="raw_data", help="Where the raw data will be written to", dest="raw_data_dir")
    parser.add_argument("--date-start", dest = "date_start", help="restrict analysis to this date at the earliest")
    parser.add_argument("--date-end", dest="date_end", help="restrict analysis to this date at the latest")
    parser.add_argument("--title", dest="title", default="Mutation report", help="report title")
    parser.add_argument("--flag-fastest", dest="flag_fastest", action="store_true", help="find the ten fastest growing genetic changes")
    
    parser.add_argument("-h","--help", action="store_true", dest="help")

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
    figdir = args.figdir
    raw_data_dir = args.raw_data_dir
    outdir = args.outdir
    date_start = args.date_start
    date_end = args.date_end
    title = args.title
    flag_fastest = args.flag_fastest

    snp_list = args.snp_list
    snp_csv = args.snp_csv

    if not os.path.exists(outdir):
        os.mkdir(os.path.join(cwd,outdir))

    if outdir not in figdir:
        figdir_writing = os.path.join(outdir,figdir)
    else:
        figdir_writing = figdir
        figdir = figdir_writing.replace(outdir,".")

    if not os.path.exists(figdir_writing):
        os.mkdir(figdir_writing)

    if outdir not in raw_data_dir:
        raw_data_dir = os.path.join(outdir,raw_data_dir)

    if not os.path.exists(raw_data_dir):
        os.mkdir(raw_data_dir)

    if args.date_start:
        date_start = dt.datetime.strptime(args.date_start,"%Y-%m-%d").date()
    else:
        date_start = dt.date(2020,1,1)

    if args.date_end:
        date_end = dt.datetime.strptime(args.date_end,"%Y-%m-%d").date()
    else:
        date_end = mfunks.find_mrd(snp_file) #by default find the most recent date in the dataset


    uk = pkg_resources.resource_filename('mink', 'data/mapping_files/gadm36_GBR_2.json')
    ni = pkg_resources.resource_filename('mink', 'data/mapping_files/NI_counties.geojson')
    channels = pkg_resources.resource_filename('mink', 'data/mapping_files/channel_islands.json')

    map_files = [uk,channels,ni]

    all_uk = map_funks.generate_all_uk_dataframe(map_files)

    group_descriptions = {}
    snps = defaultdict(list)
    snps_for_matrix = defaultdict(list)
    if snp_csv:
        print("Detected input csv")
        with open(snp_csv) as f:
            r = csv.DictReader(f)
            data = [i for i in r]
            for line in data:
                snp_list = line["snps"].split(";")
                group = line["group"]

                new_snp_list = []
                for i in snp_list:
                    if i[-1].isdigit():
                        i += "."
                    new_snp_list.append(i)

                snps[group] = new_snp_list 
                
                try:
                    description = line["description"]
                except:
                    description = ""
                
                try:
                    snps_for_matrix_list = line["snps_for_matrix"].split(";")
                except:
                    snps_for_matrix_list = None

                group_descriptions[group] = description
                snps_for_matrix[group] = snps_for_matrix_list
    
    elif snp_list: 
        print("Detected command line input")
        snp_list = snp_list.split(",")
        new_snp_list = []
        for i in snp_list:
            if i[-1].isdigit():
                i += "."
            new_snp_list.append(i)
        
        if type(args.snps_for_matrix) == str:
            snps_for_matrix_list = snps_for_matrix.split(",")
        else:
            snps_for_matrix_list = None

        snps["genetic_changes_of_interest"] = new_snp_list
        snps_for_matrix["genetic_changes_of_interest"] = snps_for_matrix_list
        group_descriptions["genetic_changes_of_interest"] = ""


    if flag_fastest:
        print("Finding fastest growing changes")
        fastest = mfunks.find_fastest_growing(snp_file, date_end)
        snps["Fastest_growers"] = fastest
        snps_for_matrix["Fastest_growers"] = fastest
        group_descriptions["Fastest_growers"] = "Amino acid changes with the highest growth rate in the thirty days preceeding the most recent sample"

    r_writer.generate_report(metadata_file, date_data, snp_file, snps, snps_for_matrix, date_start, date_end, figdir_writing, figdir, outdir, raw_data_dir, all_uk, group_descriptions, title)

if __name__ == '__main__':
    main()