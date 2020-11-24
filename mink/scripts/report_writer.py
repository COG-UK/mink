#!/usr/bin/env python3

from collections import defaultdict
from collections import Counter
from collections import OrderedDict
import matplotlib.pyplot as plt
import pandas as pd
import os
import datetime as dt
import csv
import re

import mutation_funcs as mfunk
import mapping as map_funks

import logging
logging.getLogger("imported_module").setLevel(logging.ERROR)


def process_data(metadata_file, snp_file, snp_list, date_start, date_end, snps_for_matrix, figdir_writing, figdir, outdir, raw_data_dir, all_uk, group):
    
    print(f"Running report for {group}")

    ##Data procesing and manipulation
    print("Parsing metadata")
    taxon_dict = mfunk.parse_metadata(metadata_file)
    query_to_snps, snp_to_queries, snp_list, regex_to_query = mfunk.parse_snp_data(snp_file, snp_list, taxon_dict, date_start, date_end)

    print("Making heatmap")
    if not snps_for_matrix:
        snps_for_matrix = list(snp_to_queries.keys())
    mfunk.make_heatmap(snps_for_matrix, query_to_snps, snp_to_queries, figdir_writing, group)

    print("Making maps and sorting adm2 out")
    adm2_perc_dict = defaultdict(dict)
    adm2_count_dict = defaultdict(dict)

    for snp in snp_list:
        relevant_taxa = {}
        relevant_taxa_names = snp_to_queries[snp]
        for tax in relevant_taxa_names:
            if tax in taxon_dict:
                relevant_taxa[tax] = taxon_dict[tax]
        output = map_funks.map_adm2(relevant_taxa, all_uk, figdir_writing, snp)
        if type(output) != bool:
            adm2_counter, adm2_percentages = output
            adm2_perc_dict[snp] = adm2_percentages
            adm2_count_dict[snp] = adm2_counter
        else:
            adm2_perc_dict[snp] = False

    print("Making tables")
    df, snp_to_dates, snp_last_date = mfunk.make_snp_table(snp_to_queries, taxon_dict, adm2_count_dict)
    df.to_csv(f"{outdir}/SNP_summary_table_{group}.csv")

    print("Making line figure")
    mfunk.make_overall_lines(taxon_dict,snp_to_dates,snp_last_date,figdir_writing, raw_data_dir, None, group)

    return df, adm2_perc_dict, adm2_count_dict, snp_to_queries, snp_to_dates, snp_last_date, taxon_dict, snp_list, regex_to_query

def quick_parse(snps, snp_file):

    group_to_snps = defaultdict(set)
    with open(snp_file) as f:
        reader = csv.DictReader(f)
        data = [r for r in reader]
        for line in data:
            seq_name = line["query"]
            snps_in_file = line["variants"]        

            if type(snps) == defaultdict:
                for group, snp_list in snps.items():
                    for snp in snp_list:
                        identified_snps = re.findall(snp, snps_in_file)
                        for ide in identified_snps:
                            group_to_snps[group].add(ide)

            if type(snps) == list:
                for snp in snp_list:
                    identified_snps = re.findall(snp, snps_in_file)
                    for ide in identified_snps:
                        group_to_snps["all_snps"].add(ide)

    #I'm sure I can fix this and make it neater
    if type(snps) == defaultdict:
        for group, snp_list in snps.items():
            for i in snp_list:
                if i not in group_to_snps[group] and "." not in i:
                    group_to_snps[group].add(i)
    else:
        for i in snps:
            if i not in group_to_snps["all_snps"] and "." not in i:
                group_to_snps["all_snps"].add(i)

    group_to_snps = OrderedDict(sorted(group_to_snps.items()))

    return group_to_snps

def write_start_report(title, date_data, snps, snp_file, outdir):

    print("Writing the start of the report")

    file_name = f"{outdir}/mutation_report.md"
    fw = open(file_name, 'w')
    fw.write(f"# {title}\n")
    fw.write(f'Date report run: {date_data} \n\n')
    fw.write("### Contents\n\n")

    print("finding regexed snps")
    group_to_snps = quick_parse(snps, snp_file)
    for group, snp_list in group_to_snps.items():
        if group != "all_snps":
            nice_group = group.replace("_"," ").title()
            link_group = group.replace("_","-").lower()
            fw.write(f" - [{nice_group}](#{link_group})\n")
        for i in snp_list:
            link_i = i.lower().replace(":","-")
            fw.write(f"    - [{i}](#{link_i})\n")

    fw.write("\n")

    return fw

def write_snp_sections(fw, figdir, figdir_writing, raw_data_dir, snp_df, snp_list, adm2_perc_dict, adm2_count_dict, snp_to_queries, snp_to_dates, snp_last_date,taxon_dict, description, group):
    
    if group == "all_snps":
        fw.write("## SNP summaries\n\n")
    else:
        nice_group = group.replace("_"," ").title()
        fw.write(f"## {nice_group}\n\n")
        fw.write("### Summaries\n\n")
    
    if description != "":
        nice_des = description.strip('"')
        fw.write(f"{nice_des}\n\n")

    fw.write(snp_df.to_markdown())
    fw.write("\n\n")

    fw.write("> Rolling seven day average of SNP frequency over time\n\n")
    fw.write(f'![]({figdir}/{group}_frequencies.svg)')
    fw.write("\n\n")
    
    fw.write("> Rolling seven day average of SNP counts over time\n\n")
    fw.write(f'![]({figdir}/{group}_counts.svg)')
    fw.write("\n\n")

    fw.write("> Co-occurence matrix\n\n")

    fw.write(f'![]({figdir}/pairwise_cooccurance_{group}.svg)')
    fw.write("\n\n")

    for snp in snp_list:
        
        adm2_in_map = adm2_perc_dict[snp]
        adm2_counts = adm2_count_dict[snp]
        nice_snp = snp.replace(":"," ")
        fw.write(f"### {nice_snp}\n\n")

        if len(snp_to_queries[snp]) > 10:
            small_snp_dict = defaultdict(list)
            small_snp_dict[snp] = snp_to_dates[snp]
            mfunk.make_overall_lines(taxon_dict,small_snp_dict, snp_last_date, figdir_writing, raw_data_dir, snp, group)

            fw.write(f'![]({figdir}/{snp}_frequencies.svg)')
            fw.write("\n\n")
            fw.write(f'![]({figdir}/{snp}_counts.svg)')
            fw.write("\n\n")
        else:
            if len(snp_to_queries[snp]) > 0:
                fw.write("There are fewer than ten sequences with this SNP in this time period, so count and frequency plots are not displayed.\n\n")
                fw.write("The COG IDs and dates are:\n")
                for query in snp_to_queries[snp]:
                    date = taxon_dict[query].date
                    if date != "NA":
                        fw.write(f" - {query}\t{date}\n")
                    else:
                        fw.write(f" - {query}\n")
                fw.write("\n")

        if len(snp_to_queries[snp]) > 0:
            if adm2_in_map:
                sorted_adm2_in_map =  {k: v for k, v in sorted(adm2_in_map.items(), key=lambda item: item[1], reverse=True)}
                fw.write(f"There are sequences with {snp} from {str(len(adm2_in_map))} admin2 regions\n")
                if len(adm2_in_map) > 5:
                    fw.write("The top five are:\n")
                else:
                    fw.write("This is divided into:\n")
                
                count = 0
                for adm2,percentage in sorted_adm2_in_map.items():
                    if count < 5:
                        fw.write(f'- {percentage}% ({adm2_counts[adm2]}) in {adm2}\n')
                        count += 1
                
                fw.write("\n")
                fw.write(f'![]({figdir}/{snp}_map.svg)')
                fw.write("\n\n")
            else:
                fw.write("There is no geographical data for this SNP\n\n")
        else:
            fw.write("This snp was not observed in this time period\n\n")

    return fw


def generate_report(metadata_file, date_data, snp_file, snps, snps_for_matrix, date_start, date_end, figdir_writing, figdir, outdir, raw_data_dir, all_uk, group_descriptions, title):

    #probably make snps a dictionary always and key the list by "all_snps" or something 

    fw = write_start_report(title, date_data, snps, snp_file, outdir)     

    if type(snps) == list:
        print("detected list input for snps")
        group = "all_snps"
        description = ""
        snp_df, adm2_perc_dict, adm2_count_dict, snp_to_queries, snp_to_dates, snp_last_date, taxon_dict, snp_list, regex_to_query = process_data(metadata_file, snp_file, snp_list, date_start, date_end, snps_for_matrix, figdir_writing, figdir, outdir, raw_data_dir, all_uk, group)
        fw = write_snp_sections(fw, figdir, figdir_writing, raw_data_dir, snp_df, snp_list, adm2_perc_dict, adm2_count_dict, snp_to_queries, snp_to_dates,snp_last_date, taxon_dict, description, group)

    elif type(snps) ==  defaultdict:
        print("detected input csv")
        snps = OrderedDict(sorted(snps.items()))
        for group, snp_list in snps.items():
            snps_for_matrix_actual = snps_for_matrix[group]
            group = group.replace(" ","_")
            description = group_descriptions[group]
            snp_df, adm2_perc_dict, adm2_count_dict, snp_to_queries, snp_to_dates, snp_last_date, taxon_dict, snp_list, regex_to_query = process_data(metadata_file, snp_file, snp_list, date_start, date_end, snps_for_matrix_actual, figdir_writing, figdir, outdir, raw_data_dir, all_uk, group)
            fw = write_snp_sections(fw, figdir, figdir_writing, raw_data_dir, snp_df, snp_list, adm2_perc_dict, adm2_count_dict, snp_to_queries, snp_to_dates,snp_last_date, taxon_dict, description, group)

    fw.close()
    
