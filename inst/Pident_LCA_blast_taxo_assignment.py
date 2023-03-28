#!/usr/bin/env python


'''
This scripts uses python 3 and the following libraries need to be installed 'pandas', 'ete3' and 'argparse'.

The blastn file needs no modification. As long as it is the blastn format 6 output with options 
"query.id", "query.length", "pident", "subject.id", "subject.GBid", "evalue", "bit.score","staxids", "sscinames", "sblastnames", "qcovs", "qcovhsp" 

If you use the taxonly option because you do not want to or do not have the ASV/OTU table, the output taxonomy file may be missing some entries compared to the number of ASVs/OTUs you are expecting. This is because some ASVs/OTUs will either have no hits with blastn or that no taxonomy could be returned (e.g. minimun coverage not respected, conflicting taxonomy at Kingdom level, etc.).

'''

import pandas as pd
import csv
from ete3 import NCBITaxa
import re
import argparse
import numpy as np
import sys
import time
from tqdm import tqdm


# Parser program
def handle_program_options():
    """Parses the given options passed in at the command line."""

    parser = argparse.ArgumentParser(
        description='Takes blastn multiple hits, trim uncultured or unidentified hits or environmental samples, assign taxo to feature based on percent similarity (for species) and Last Common Ancestor method')
    parser.add_argument('-b', "--btbl", required=True,
                        help='blastn format 6 output with options "query.id", "query.length", "pident", "subject.id", "subject.GBid", "evalue", "bit.score","staxids", "sscinames", "sblastnames", "qcovs", "qcovhsp", [REQUIRED]')
    parser.add_argument("-f", "--ftbl", required=False,
                        help="Feature id table in txt format [REQUIRED]")
    parser.add_argument("-o", "--output_ftbl", required=True,
                        help="Output path for the transformed feature_id table. [REQUIRED]")
    parser.add_argument('--minSim', default=97, type=int,
                        help='Minimum similarity to assign species hits (default: 97)')
    parser.add_argument('--minCov', default=80, type=int,
                        help='Minimum coverage to keep hits (default: 80)')
    parser.add_argument('--update', default="False", type=str,
                        help='Should the taxonomy database be updated')
    parser.add_argument('--pident', default='no', type=str,
                        help='To reduce taxonomy assingment according to default percent identity thresholds. Options are: before or after LCA assingment')
    parser.add_argument('--pgenus', default=95, type=int,
                        help='Minimum similarity to assign genus (default: 95)')
    parser.add_argument('--pfamily', default=87, type=int,
                        help='Minimum similarity to assign family (default: 87)')
    parser.add_argument('--porder', default=83, type=int,
                        help='Minimum similarity to assign order (default: 83)')
    parser.add_argument('--pclass', default=81, type=int,
                        help='Minimum similarity to assign class (default: 81)')
    parser.add_argument('--pphylum', default=79, type=int,
                        help='Minimum similarity to assign phylum (default: 79)')
    parser.add_argument('--pkingdom', default=71, type=int,
                        help='Minimum similarity to assign kingdom (default: 71)')                        
    parser.add_argument('--taxonly', default="True", type=str,
                        help='Do not require the ASV/OTU table')
    parser.add_argument('-v', '--verbose', action='store_true')

    return parser.parse_args()

# Progress bar
def progress(count, total, status=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))
    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush()

# Function to retrieve taxonomy
def get_desired_ranks(ncbi, taxid, desired_ranks):
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = dict((rank, taxid)
                         for (taxid, rank) in lineage2ranks.items())
    taxo = []
    for rank in desired_ranks:
        if rank in ranks2lineage:
            taxo.append(names[ranks2lineage[rank]])
        else:
            taxo.append("NA")
    return ";".join(taxo)

# Insert taxonomy into blast table
def taxo_assignment(tabl, dict_blast):
    print('\nAssigning full taxonomy to blast table')
    taxo = []
    for index, row in tabl.iterrows():
        taxo.append(dict_blast[tabl.loc[index]['staxids']])
    print('\nAssigning full taxonomy to blast table is completed')
    return taxo

# PIDENT
# Function to reduce taxonomy assignment depending on pident value
def pidentThresholds(row, minSim,pkingdom,pphylum,pclass,porder,pfamily,pgenus):
    if row['pident'] < pkingdom:
        row[["kingdom", 'phylum', 'class', 'order', 'family', 'genus',
             'species']] = "NA", "NA", "NA", "NA", "NA", "NA", "NA"
        return row
    elif row['pident'] < pphylum:
        row[["kingdom", 'phylum', 'class', 'order', 'family', 'genus',
             'species']] = row["kingdom"], "NA", "NA", "NA", "NA", "NA", "NA"
        return row
    elif row['pident'] < pclass:
        row[["kingdom", 'phylum', 'class', 'order', 'family', 'genus', 'species']
            ] = row["kingdom"], row['phylum'], "NA", "NA", "NA", "NA", "NA"
        return row
    elif row['pident'] < porder:
        row[["kingdom", 'phylum', 'class', 'order', 'family', 'genus', 'species']
            ] = row["kingdom"], row['phylum'], row['class'], "NA", "NA", "NA", "NA"
        return row
    elif row['pident'] < pfamily:
        row[["kingdom", 'phylum', 'class', 'order', 'family', 'genus', 'species']
            ] = row["kingdom"], row['phylum'], row['class'], row['order'], "NA", "NA", "NA"
        return row
    elif row['pident'] < pgenus:
        row[["kingdom", 'phylum', 'class', 'order', 'family', 'genus', 'species']
            ] = row["kingdom"], row['phylum'], row['class'], row['order'], row['family'], "NA", "NA"
        return row
    elif row['pident'] < minSim:
        row[["kingdom", 'phylum', 'class', 'order', 'family', 'genus', 'species']
            ] = row["kingdom"], row['phylum'], row['class'], row['order'], row['family'], row['genus'], "NA"
        return row
    else:
        row[["kingdom", 'phylum', 'class', 'order', 'family', 'genus', 'species']
            ] = row["kingdom"], row['phylum'], row['class'], row['order'], row['family'], row['genus'], row['species']
        return row

# LCA
# If similarity of best hit => minSim%, assign to species level, otherwise assign to last common ancestor
def taxo_consensus(tabl, tabl2, minSim):
    new = tabl
    #new['species'] = ["" if new[new.index == ind]["pident"].iat[0] < minSim else new[new.index == ind]["species"].iat[0] for ind in new.index]

    def Remove(sets):
        sets.discard("")
        return(sets)

    rankLevel = 0
    listRanks = ['species', 'genus', 'family', 'order',
                 'class', 'phylum', 'kingdom', 'superkingdom']
    #t0 = time.time()
    for i in tqdm(range(len(listRanks))):
        #t1 = time.time()
        #progress(rankLevel,len(listRanks[i]), status = "Progress")
        for query, row in tqdm(new.iterrows()):
            setTaxo = set(tabl2[tabl2['query.id'] == query][listRanks[i]])
            setTaxo = Remove(setTaxo)
            if row['pident'] < minSim and len(setTaxo) > 1:
                new.loc[query, listRanks[i]] = ""
                x = rankLevel
                while x > 0:
                    new.loc[query, listRanks[i-x]] = ""
                    x -= 1
            elif row['pident'] < minSim:
                s = list(setTaxo)
                s = ['' if v is None else v for v in s]
                s = ''.join(str(s))
                new.loc[query, listRanks[i]] = s
        rankLevel += 1
        #t2 = time.time()
        #print(" {}s (estimated remaining time: {}m,s)".format(t2-t1, t2-t0))

    for query, row in new.iterrows():
        a = new.columns.get_loc('superkingdom')
        b = new.columns.get_loc('species')
        c = str(row[a:b + 1].str.cat(sep=';'))
        c = c.replace("{", "").replace("}", "").replace("[", "").replace("]", "").replace("'", "").replace(
            " ,", "").replace(", ", "").replace(",", "").replace("NA", "").replace("nan", "")
        c = re.sub(' sp\..*', ' sp.', c)
        # need to correspond to number of ranks...
        c = re.sub('^;;;;;;;', 'Unknown', c)
        new.loc[query, 'taxonomy'] = c
    return new


# Functions for taxo assignment based on any of the 3 options
def pident_bef_LCA(b_trimmed, minSim, pkingdom, pphylum, pclass ,porder, pfamily, pgenus):
    b_trimmed = b_trimmed[b_trimmed.taxonomy != "NA"]
    b_trimmed = b_trimmed.apply(pidentThresholds, args = (minSim,pkingdom,pphylum,pclass,porder,pfamily,pgenus), axis=1)
    b_trimmed = b_trimmed.replace(r'NA', "", regex=True)
    dummy2 = b_trimmed.groupby('query.id', group_keys=False).apply(
        lambda x: x.loc[x.evalue.idxmin()])
    f_btbl = taxo_consensus(dummy2, b_trimmed, minSim)
    print('\nPident trimming and LCA completed')
    return f_btbl


def LCA_bef_pident(b_trimmed, minSim, pkingdom, pphylum, pclass, porder, pfamily, pgenus):
    b_trimmed = b_trimmed.replace(r'NA', np.nan, regex=True)
    dummy2 = b_trimmed.groupby('query.id', group_keys=False).apply(
        lambda x: x.loc[x.evalue.idxmin()])
    f_btbl = taxo_consensus(dummy2, b_trimmed, minSim)
    # Pident thresholds
    f_btbl2 = f_btbl[f_btbl.taxonomy != "NA"]
    f_btbl = f_btbl2.apply(pidentThresholds, args = (minSim, pkingdom, pphylum, pclass, porder, pfamily, pgenus), axis=1)
    #f_btbl = f_btbl.replace([None], ['NA'], regex=True)
    f_btbl = f_btbl.replace(r'^\s*$', 'NA', regex = True)
    f_btbl = f_btbl.replace(np.nan, 'NA', regex=True)
    f_btbl['taxonomy'] = f_btbl[['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']].apply(lambda x: ';'.join(x), axis=1)
    f_btbl['taxonomy'] = f_btbl.taxonomy.str.replace("[", "", regex=False).str.replace("]", "", regex=False).str.replace("'", "", regex=False).str.replace(" ,", "", regex=False).str.replace(", ", "", regex=False).str.replace(",", "", regex=False).str.replace("nan", "", regex=False).str.replace("NA", "", regex=False).str.replace("^;;;;;;;", "Unknown", regex=True)
    print('\nLCA and Pident trimming completed')
    return f_btbl

def LCA_only(b_trimmed, minSim):
    b_trimmed = b_trimmed.replace(r'NA', np.nan, regex=True)
    dummy2 = b_trimmed.groupby('query.id', group_keys=False).apply(
        lambda x: x.loc[x.evalue.idxmin()])
    f_btbl = taxo_consensus(dummy2, b_trimmed, minSim)
    print('\nLCA completed')
    return f_btbl


# Assign taxonomy to each feature-id
def blast_to_feature_tbl(btbl, ftbl):
    print('\nAssign full taxonomy to each feature-id')
    new_ftbl = ftbl
    new_ftbl['taxonomy'] = ""
    new_ftbl["Percent_Identity"] = ""
    new_ftbl["Sequence_coverage"] = ""
    for query, row in new_ftbl.iterrows():
        if row['#OTU ID'] in list(btbl['query.id']):
            otu_x = row['#OTU ID']
            a = btbl[btbl['query.id'] == otu_x]['taxonomy'].iat[0]
            b = btbl[btbl['query.id'] == otu_x]['pident'].iat[0]
            c = btbl[btbl['query.id'] == otu_x]['qcovs'].iat[0]
            new_ftbl.at[query, 'taxonomy'] = a
            new_ftbl.at[query, 'Percent_Identity'] = round(b)
            new_ftbl.at[query, 'Sequence_coverage'] = round(c)

        else:
            new_ftbl.at[query, 'taxonomy'] = "No hit"

    # To return only taxo assingment information
    new_ftbl = new_ftbl[["#OTU ID", "taxonomy",
                         "Percent_Identity", "Sequence_coverage"]]
    return new_ftbl

# Create taxonomy table only
def blast_to_taxonomy_tbl(btbl, ftbl):
    print('\nAssign full taxonomy to each feature-id')
    #new_ftbl = ftbl["ASVs"].to_frame(name=None)
    new_ftbl = ftbl["ASVs"].to_frame()   
    new_ftbl['taxonomy'] = ""
    new_ftbl["Percent_Identity"] = ""
    new_ftbl["Sequence_coverage"] = ""
    for query, row in new_ftbl.iterrows():
        if row['ASVs'] in list(btbl['query.id']):
            otu_x = row['ASVs']
            a = btbl[btbl['query.id'] == otu_x]['taxonomy'].iat[0]
            b = btbl[btbl['query.id'] == otu_x]['pident'].iat[0]
            c = btbl[btbl['query.id'] == otu_x]['qcovs'].iat[0]
            new_ftbl.at[query, 'taxonomy'] = a
            new_ftbl.at[query, 'Percent_Identity'] = round(b)
            new_ftbl.at[query, 'Sequence_coverage'] = round(c)

        else:
            new_ftbl.at[query, 'taxonomy'] = "No hit"
    # To return only taxo assingment information
    new_ftbl = new_ftbl[["ASVs", "taxonomy",
                         "Percent_Identity", "Sequence_coverage"]]
    return new_ftbl




##############################        MAIN        #################################################
def main():
    args = handle_program_options()
    taxonly = args.taxonly.lower()
    update = args.update.lower()
    blast = pd.read_table(args.btbl,
                          names=["query.id", "query.length", "pident", "subject.id", "subject.GBid", "evalue",
                                 "bit.score", "staxids", "sscinames", "sblastnames", "qcovs", "qcovhsp"])

    ncbi = NCBITaxa()
    if update == "true":
        ncbi.update_taxonomy_database()
    else:
        pass

    # 1-  Remove uncultured and unidentified, and remove hit with coverage less than 80%
    blast_trimmed = blast[~blast.sscinames.str.contains(
        'uncultured', na=False)]  # May have to remove
    blast_trimmed = blast_trimmed[~blast_trimmed.sscinames.str.contains(
        'unidentified', na=False)]
    blast_trimmed = blast_trimmed[~blast_trimmed.sscinames.str.contains(
        'environmental sample', na=False)]
    blast_trimmed = blast_trimmed[blast_trimmed.qcovs >= args.minCov]

    # 2- Make a list of all unique staxids
    blast_trimmed.staxids = blast_trimmed.staxids.apply(str)
    blast_trimmed.staxids = blast_trimmed.staxids.str.split(';').str[0]
    l_staxids = list(set(blast_trimmed['staxids']))

    # 3- Create a dictionary with staxids as Key and taxonomy as value
    desired_ranks = ['superkingdom', "kingdom", 'phylum',
                     'class', 'order', 'family', 'genus', 'species']
    dict_blast = {}
    for i in l_staxids:
        try:
            dict_blast[i] = get_desired_ranks(ncbi, i, desired_ranks)
        except ValueError:
            dict_blast[i] = "NA"

    # 4- Insert taxonomy into blast table
    blast_trimmed['taxonomy'] = taxo_assignment(blast_trimmed, dict_blast)
    dummy1 = blast_trimmed['taxonomy'].str.split(';', expand=True)
    dummy1.columns = ['superkingdom', "kingdom", 'phylum',
                      'class', 'order', 'family', 'genus', 'species']
    blast_trimmed['superkingdom'] = dummy1['superkingdom']
    blast_trimmed['kingdom'] = dummy1['kingdom']
    blast_trimmed['phylum'] = dummy1['phylum']
    blast_trimmed['class'] = dummy1['class']
    blast_trimmed['order'] = dummy1['order']
    blast_trimmed['family'] = dummy1['family']
    blast_trimmed['genus'] = dummy1['genus']
    blast_trimmed['species'] = dummy1['species']

    # 5-6 If pidentThresholds is True, reduce taxo assingment using pident threshold values
    # Insert taxonomy into blast table
    if args.pident == "before":
        print('\nReducing taxonomy resolution based on percent identity thresholds and then LCA')
        final_btbl = pident_bef_LCA(b_trimmed = blast_trimmed, minSim = args.minSim, pkingdom = args.pkingdom, pphylum = args.pphylum, pclass = args.pclass, porder = args.porder, pfamily = args.pfamily, pgenus = args.pgenus)
    elif args.pident == "after":
        print('\nReducing taxonomy resolution based on LCA and then percent identity thresholds')
        final_btbl = LCA_bef_pident(b_trimmed = blast_trimmed, minSim = args.minSim, pkingdom = args.pkingdom, pphylum = args.pphylum, pclass = args.pclass, porder = args.porder, pfamily = args.pfamily, pgenus = args.pgenus)

    else:
        # LCA assingment. If similarity of best hit => 97%, assign to species level, otherwise assign to last common ancestor
        print('\nReducing taxonomy resolution based on LCA only')
        final_btbl = LCA_only(blast_trimmed, minSim = args.minSim)

    # 7- Assign taxonomy to each feature-id
    if taxonly == "false":
        try:
            # feature_table = pd.read_table(args.ftbl, skiprows=[0]) # Will skip the 'Constructed from biom file'
            # Will skip the 'Constructed from biom file'
            feature_table = pd.read_table(args.ftbl)

        except SimulationException as sim_exc:
            print("No feature table provided", sim_exc)

        feature_table = blast_to_feature_tbl(final_btbl, feature_table)
    else:
        feature_table = final_btbl['query.id'].unique()
        feature_table = pd.DataFrame(data=feature_table, columns=["ASVs"])
        feature_table = blast_to_taxonomy_tbl(final_btbl, feature_table)

    # 8- Write output to csv
    pd.DataFrame.to_csv(feature_table, args.output_ftbl, index=False)
    print("\nProgram completed with success!!!\n")
    if args.verbose:
        print("\nTransformed table written to: {}".format(args.output_ftbl))


if __name__ == '__main__':
    main()
