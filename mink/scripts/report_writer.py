#!/usr/bin/env python3

from collections import defaultdict
from collections import Counter
from collections import OrderedDict
import matplotlib.pyplot as plt
import pandas as pd
import os
import datetime as dt

import mutation_funcs as mfunk
import mapping as map_funks


def process_data(metadata_file, snp_file, snp_list, date_start, date_end, snps_for_matrix, figdir_writing, figdir, outdir, raw_data_dir, all_uk):

    ##Data procesing and manipulation
    print("Parsing metadata")
    taxon_dict = mfunk.parse_metadata(metadata_file)
    query_to_snps, snp_to_queries = mfunk.parse_snp_data(snp_file, snp_list, taxon_dict, date_start, date_end)

    print("Making heatmap")
    mfunk.make_heatmap(snps_for_matrix, query_to_snps, figdir_writing)

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
    df.to_csv(f"{outdir}/SNP_summary_table.csv")

    print("Making line figure")
    mfunk.make_overall_lines(taxon_dict,snp_to_dates,snp_last_date,figdir_writing, raw_data_dir, None)

    return df, adm2_perc_dict, adm2_count_dict, snp_to_queries, snp_to_dates, snp_last_date, taxon_dict


def write_report(outdir, date_data, figdir, figdir_writing, raw_data_dir, snp_df, snp_list, adm2_perc_dict, adm2_count_dict, snp_to_queries, snp_to_dates, snp_last_date,taxon_dict):
   
    ## Writing the report ##
    print("Writing the report")
    fw = open(f"{outdir}/mutation_report.md", 'w')

    fw.write("# Mutation report\n")

    fw.write(f'Date report run: {date_data} \n\n')

    fw.write("Please note that adm2s can be joined together if there are sequences with ambiguous location data, so the adm2 number should only be used as an estimate.\n\n")

    fw.write(snp_df.to_markdown())
    fw.write("\n\n")

    fw.write("## SNP summaries\n")
    fw.write("### COG-UK amino acid variants by date\n\n")

    fw.write(f'![]({figdir}/all_snps_line.svg)')
    fw.write("\n\n")

    fw.write("### Co-occurence matrix\n\n")

    fw.write(f'![]({figdir}/pairwise_cooccurance.svg)')
    fw.write("\n\n")

    for snp in snp_list:
        
        adm2_in_map = adm2_perc_dict[snp]
        adm2_counts = adm2_count_dict[snp]
        fw.write(f"### {snp}\n\n")

        if len(snp_to_queries[snp]) > 0:
            small_snp_dict = defaultdict(list)
            small_snp_dict[snp] = snp_to_dates[snp]
            mfunk.make_overall_lines(taxon_dict,small_snp_dict, snp_last_date, figdir_writing, raw_data_dir, snp)

            fw.write(f'![]({figdir}/{snp}_line.svg)')
            fw.write("\n\n")

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

    fw.close()

def generate_report(metadata_file, date_data, snp_file, snp_list, date_start, date_end, snps_for_matrix, figdir_writing, figdir, outdir, raw_data_dir, all_uk):

    snp_df, adm2_perc_dict, adm2_count_dict, snp_to_queries, snp_to_dates, snp_last_date, taxon_dict = process_data(metadata_file, snp_file, snp_list, date_start, date_end, snps_for_matrix, figdir_writing, figdir, outdir, raw_data_dir, all_uk)
    write_report(outdir, date_data, figdir, figdir_writing, raw_data_dir, snp_df, snp_list, adm2_perc_dict, adm2_count_dict, snp_to_queries, snp_to_dates,snp_last_date, taxon_dict)
