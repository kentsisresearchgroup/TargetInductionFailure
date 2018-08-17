# mergeSamples v3
# by John Philip
# 20170810
#
# This takes a directory full of files from NYGC and merges them into one master file.
# Uses updated Kentsis Suggestions:
# SNVs: group by genes, sort by recurrence in descending order.
# CVs: extract out all genes names: contained, left, right. Make a row for each, keeping annotations. Group by Genes
#      and sort by recurrence in descending order
# Fusions: extract out all genes, sort by recurrence in descending order
# Indels: extract out all genes, sort by recurrence in descending order
#
# do all this with and without SNEPP_EFFECT filter
# added pivot tables
#

import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pdb

def flatten(l, ltypes=(list, tuple)):
    # code from http://rightfootin.blogspot.com/2006/09/more-on-python-flatten.html
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)


def calc_sample_name_from_filename(filename):
    # get sample name from File name

    sample_name = filename[filename.find('TARGET-21-'):filename.find('--TARGET-21')]
    sample_name = sample_name[:sample_name.find('_')]
    # print(sample_name)
    return sample_name


def getINDELS(reportlist):

    maintable_df = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'CALLED_BY', 'SNPEFF_GENE_NAME',
                                         'COSMIC_GENE', 'COSMIC_AA_CHANGE', 'COSMIC_CDS', 'COSMIC_CNT', '1000G_AF',
                                         'ExAC_AF', 'SNPEFF_IMPACT', 'SNPEFF_EFFECT', 'SNPEFF_FUNCTIONAL_CLASS',
                                         'SNPEFF_AA_CHANGE', 'SNPEFF_CDS_CHANGE', 'SNPEFF_CODON_CHANGE',
                                         'SNPEFF_EXON_ID', 'SNPEFF_GENE_BIOTYPE', 'SNPEFF_TRANSCRIPT_ID',
                                         'TARGET_SCAN_miR_TARGET', 'MutationAssessor_score', 'MutationAssessor_pred',
                                         'FATHMM_SOMATIC_score', 'FATHMM_SOMATIC_pred', 'CHASM_score', 'CHASM_pred',
                                         'Actionable', 't_alt_count', 't_ref_count', 'n_alt_count', 'n_ref_count'])

    for report in reportlist:
        t_df = pd.read_csv(report, delimiter='\t')
        t_df['Sample'] = calc_sample_name_from_filename(report)
        # print(t_df.shape)
        maintable_df = maintable_df.append(t_df, ignore_index=True)

    # print(F_list)
    # print(V_list)
    print(maintable_df.shape)
    maintable_df['mutation'] = 'indels'

    return maintable_df


def getSNVs(reportlist):

    maintable_df = pd.DataFrame(columns=['CHROM','POS','ID','REF','ALT','CALLED_BY','SNPEFF_GENE_NAME','COSMIC_GENE',
                                         'COSMIC_AA_CHANGE','COSMIC_CDS','COSMIC_CNT','1000G_AF','ExAC_AF',
                                         'SNPEFF_IMPACT','SNPEFF_EFFECT','SNPEFF_FUNCTIONAL_CLASS','SNPEFF_AA_CHANGE',
                                         'SNPEFF_CDS_CHANGE','SNPEFF_CODON_CHANGE','SNPEFF_EXON_ID',
                                         'SNPEFF_GENE_BIOTYPE','SNPEFF_TRANSCRIPT_ID','TARGET_SCAN_miR_TARGET',
                                         'MutationAssessor_score','MutationAssessor_pred','FATHMM_SOMATIC_score',
                                         'FATHMM_SOMATIC_pred','CHASM_score','CHASM_pred','Actionable','t_alt_count',
                                         't_ref_count','n_alt_count','n_ref_count'])

    for report in reportlist:
        t_df = pd.read_csv(report, delimiter='\t')

        # get sample name from File name
        t_df['Sample'] = calc_sample_name_from_filename(report)

        # print(t_df.shape)
        maintable_df = maintable_df.append(t_df, ignore_index=True)

    print(maintable_df.shape)
    maintable_df['mutation'] = 'SNVs'

    return maintable_df


def getSVs(reportlist):

    maintable_df = pd.DataFrame()

    for report in reportlist:
        t_df = pd.read_csv(report, delimiter='\t', header=None)

        # get sample name from File name
        t_df['Sample'] = calc_sample_name_from_filename(report)

        #print(t_df.shape)
        maintable_df = maintable_df.append(t_df, ignore_index=True)

    print(maintable_df.shape)
    maintable_df['mutation'] = 'SVs'

    return maintable_df


def getCNVs(reportlist):

    maintable_df = pd.DataFrame()

    for report in reportlist:
        t_df = pd.read_csv(report, delimiter='\t', header=None)

        # get sample name from File name
        t_df['Sample'] = calc_sample_name_from_filename(report)

        # print(t_df.shape)
        maintable_df = maintable_df.append(t_df, ignore_index=True)

    print(maintable_df.shape)
    maintable_df['mutation'] = 'CNVs'

    return maintable_df


def getFusions(reportlist):

    maintable_df = pd.DataFrame(columns=['Sample', 'Gene1', 'Gene2', 'Description', 'Spanning_pairs',
                                         'Spanning_unique_reads', 'Breakpoint1', 'Breakpoint2', 'GeneID1', 'GeneID2',
                                         'Fusion_predicted_effect'])

    for report in reportlist:
        # print(report)
        t_df = pd.read_csv(report, delimiter='\t')
        # print(t_df.shape)
        maintable_df = maintable_df.append(t_df, ignore_index=True)

    maintable_df['mutation'] = 'fusions'

    return maintable_df


def process_SNVs(SNVs):

    genes_to_look_df = pd.read_csv('GenestoLook.txt', header=None)
    print(genes_to_look_df.shape)
    genes_to_look = genes_to_look_df[0].tolist()

    SNV_df = SNVs

    SNV_df['myGeneName'] = SNV_df['SNPEFF_GENE_NAME']
    SNV_df['SampleName'] = SNV_df.Sample.str[10:]
    SNV_df['mutation_Name'] = 'SNV_' + SNV_df['SNPEFF_GENE_NAME'].map(str) + '_AA:' + \
                              SNV_df['SNPEFF_AA_CHANGE'].map(str) + '_CDS:' + SNV_df['SNPEFF_CDS_CHANGE'].map(str) + \
                              '_' + SNV_df['SNPEFF_EFFECT'].map(str)

    SNV_df['ToLook'] = SNV_df.apply(lambda row: 1 if (row.SNPEFF_GENE_NAME in genes_to_look) else 0, axis=1)
    maintable_writer = pd.ExcelWriter('SNVs_unfiltered.xlsx', engine='xlsxwriter')
    SNV_df.to_excel(maintable_writer)
    print(SNV_df.shape)
    maintable_writer.close()

    # filter for only mutations called by more than one caller.
    # SNV_df['keep'] = SNV_df.apply(lambda row: 1 if (row.CALLED_BY == 'Intersection') or
    #                                                (row.CALLED_BY == 'mutect-strelka') or
    #                                                (row.CALLED_BY == 'mutect-lofreq') or
    #                                                (row.CALLED_BY == 'lofreq-strelka') else 0, axis=1)

    #SNV_df = SNV_df[SNV_df.keep == 1]
    #SNV_df = SNV_df.drop('keep',axis=1)

    SNV_df['keep'] = SNV_df.apply(lambda row: 1 if (row.SNPEFF_IMPACT == 'HIGH') or
                                                   (row.SNPEFF_IMPACT == 'MODERATE') or
                                                   (row.SNPEFF_IMPACT == 'LOW') else 0, axis=1)
    SNV_df = SNV_df[SNV_df.keep == 1]
    SNV_df = SNV_df.drop('keep',axis=1)


    # filter SNVs further
    filter_terms = set(['NON_SYNONYMOUS_CODING', 'NON_SYNONYMOUS_CODING+SPLICE_SITE_REGION',
                        'SPLICE_SITE_ACCEPTOR+INTRON', 'SPLICE_SITE_DONOR+INTRON', 'SPLICE_SITE_REGION+INTRON',
                        'SPLICE_SITE_REGION+NON_CODING_EXON_VARIANT', 'SPLICE_SITE_REGION+START_GAINED+UTR_5_PRIME',
                        'SPLICE_SITE_REGION+SYNONYMOUS_CODING', 'SPLICE_SITE_REGION+UTR_3_PRIME',
                        'SPLICE_SITE_REGION+UTR_5_PRIME', 'START_GAINED+UTR_5_PRIME', 'START_LOST', 'STOP_GAINED',
                        'STOP_GAINED+SPLICE_SITE_REGION', 'STOP_LOST'])

    SNV_df = SNV_df.dropna(subset=['SNPEFF_EFFECT'])

    SNV_df['mark4deletion'] = SNV_df['SNPEFF_EFFECT'].apply(lambda x:len(set(x.split(',')).intersection(filter_terms)))
    SNV_df = SNV_df[SNV_df.mark4deletion != 0]
    SNV_df = SNV_df.drop('mark4deletion', axis=1)

    SNV_df = SNV_df[SNV_df.t_alt_count >= 8]

    SNV_df['n_AF'] = SNV_df['n_alt_count']/(SNV_df['n_alt_count'] + SNV_df['n_ref_count'])
    SNV_df['t_AF'] = SNV_df['t_alt_count']/(SNV_df['t_alt_count'] + SNV_df['t_ref_count'])

    maintable_writer = pd.ExcelWriter('SNVs_filtered.xlsx', engine='xlsxwriter')
    SNV_df.to_excel(maintable_writer)
    print(SNV_df.shape)
    maintable_writer.close()


def process_indels(indels):

    genes_to_look_df = pd.read_csv('GenestoLook.txt', header=None)
    print("genes_to_look_df.shape: ", genes_to_look_df.shape)
    genes_to_look = genes_to_look_df[0].tolist()

    indel_df = indels
    indel_df['myGeneName'] = indel_df['SNPEFF_GENE_NAME']

    indel_df['ToLook'] = indel_df.apply(lambda row: 1 if (row.SNPEFF_GENE_NAME in genes_to_look) else 0, axis=1)

    indel_df['SampleName'] = indel_df.Sample.str[10:]
    indel_df['mutation_Name'] = 'indel_' + indel_df['SNPEFF_GENE_NAME'].map(str) + '_AA:' + \
                                indel_df['SNPEFF_AA_CHANGE'].map(str) + '_CDS:' + \
                                indel_df['SNPEFF_CDS_CHANGE'].map(str) + '_' + indel_df['SNPEFF_EFFECT'].map(str)
    # Keep calls with only high read counts
    indel_df = indel_df[indel_df.t_alt_count > 7]

    maintable_writer = pd.ExcelWriter('indel_unfiltered.xlsx', engine='xlsxwriter')
    indel_df.to_excel(maintable_writer)
    maintable_writer.close()

    # filter for only mutations that are moderate or high impact
    indel_df['keep'] = indel_df.apply(lambda row: 1 if (row.SNPEFF_IMPACT == 'HIGH') or
                                                   (row.SNPEFF_IMPACT == 'MODERATE') else 0, axis=1)

    indel_df = indel_df[indel_df.keep == 1]
    indel_df = indel_df.drop('keep', axis=1)

    indel_df['t_AF'] = indel_df['t_alt_count']/(indel_df['t_alt_count'] + indel_df['t_ref_count'])
    indel_df['n_AF'] = indel_df['n_alt_count']/(indel_df['n_alt_count'] + indel_df['n_ref_count'])

    maintable_writer = pd.ExcelWriter('indel_filtered.xlsx', engine='xlsxwriter')
    indel_df.to_excel(maintable_writer)
    maintable_writer.close()


def process_CNV(CNVs):

    CNV_df = CNVs
    CNV_df['mutation_Name'] = CNV_df[0].map(str) + '_' + CNV_df[1].map(str) + '_' + CNV_df[2].map(str) + '_' \
                              + CNV_df[3]

    CNV_df['myGeneName'] = ''

    numrows = CNV_df.shape[0]
    print("numrows: ", numrows)

    tempDF = pd.DataFrame()

    for i in range (0,numrows):
        geneStr = CNV_df.loc[i,6]

        genelist = ''
        disrupt_l = ''
        disrupt_r = ''

        str_start = geneStr.find('Cancer_census=')
        if (str_start > -1):  # Cancer Census found
            str_start = str_start + 14
            str_end = geneStr.find(';', str_start)
            genelist = geneStr[str_start:str_end]

        if (genelist != ''):
            genelist = genelist.split(',')
        else:
            genelist = []

        str_start = geneStr.find('DisruptL=')
        if (str_start > -1):
            str_start = str_start + 9
            str_end = geneStr.find(';', str_start)
            disrupt_l = geneStr[str_start:str_end]

        if (disrupt_l != ''):
            disrupt_l = disrupt_l.split(',')
            genelist.append(disrupt_l)

        str_start = geneStr.find('DisruptR=')
        if (str_start > -1):
            str_start = str_start + 9
            str_end = geneStr.find(';', str_start)
            disrupt_r = geneStr[str_start:str_end]

        if (disrupt_r != ''):
            disrupt_r = disrupt_r.split(',')
            genelist.append(disrupt_r)

        if len(genelist)>0:
            genelist = flatten(genelist)
        #print(genelist)
        #print ('i = ' + str(i))

        genelist_len = len(genelist)
        current_row = CNV_df.loc[i]

        t_DF = pd.DataFrame([current_row] * genelist_len)
        t_DF = t_DF.reset_index()

        for j in range(0,genelist_len):
            t_DF.loc[j,'myGeneName'] = genelist[j]
            #print(tempDF.shape)

        #print(genelist)
        #print(genelist_len)
        #print(t_DF)

        tempDF = tempDF.append(t_DF, ignore_index=True)

    CNV_df = CNV_df.append(tempDF, ignore_index=True)

    CNV_df['SampleName'] = CNV_df.Sample.str[10:]
    CNV_df['mutation_Name'] = 'CNV_' + CNV_df['myGeneName'].map(str) + ':' + CNV_df['mutation_Name'].map(str)

    print("CNV_df.shape before filtering: ", CNV_df.shape)
    maintable_writer = pd.ExcelWriter('CNVs_unfiltered.xlsx', engine='xlsxwriter')
    CNV_df.to_excel(maintable_writer)
    maintable_writer.close()

    CNV_df = CNV_df[CNV_df[3] != 'NEU']
    CNV_df.reset_index(inplace=True)
    CNV_df.drop('level_0',axis=1, inplace=True)
    CNV_df = CNV_df[CNV_df['myGeneName'] != '']
    CNV_df = CNV_df.reset_index()

    print("CNV_df.shape after filtering: ", CNV_df.shape)
    maintable_writer = pd.ExcelWriter('CNVs_filtered.xlsx', engine='xlsxwriter')
    CNV_df.to_excel(maintable_writer)
    maintable_writer.close()


def process_SVs(SVs):

    SV_df = SVs
    SV_df['mutation_Name'] = 'SV_' + SV_df[0].map(str) + ':' + SV_df[1].map(str) + ':' + SV_df[2].map(str) + ':' \
                             + SV_df[6]

    SV_df['myGeneName'] = ''

    SV_df.reset_index(inplace=True)
    numrows = SV_df.shape[0]
    print("numrows: ", numrows)

    tempDF = pd.DataFrame()

    for i in range(0, numrows):
        geneStr = SV_df.loc[i, 11]

        genelist = ''
        disrupt_l = ''
        disrupt_r = ''

        str_start = geneStr.find('Cancer_census=')
        if (str_start > -1):  # Cancer Census found
            str_start = str_start + 14
            str_end = geneStr.find(';', str_start)
            genelist = geneStr[str_start:str_end]

        if (genelist != ''):
            genelist = genelist.split(',')
        else:
            genelist = []

        str_start = geneStr.find('DisruptL=')
        if (str_start > -1):
            str_start = str_start + 9
            str_end = geneStr.find(';', str_start)
            disrupt_l = geneStr[str_start:str_end]

        if (disrupt_l != ''):
            disrupt_l = disrupt_l.split(',')
            genelist.append(disrupt_l)

        str_start = geneStr.find('DisruptR=')
        if (str_start > -1):
            str_start = str_start + 9
            str_end = geneStr.find(';', str_start)
            disrupt_r = geneStr[str_start:str_end]

        if (disrupt_r != ''):
            disrupt_r = disrupt_r.split(',')
            genelist.append(disrupt_r)

        if len(genelist) > 0:
            genelist = flatten(genelist)

        genelist_len = len(genelist)
        current_row = SV_df.loc[i]

        t_DF = pd.DataFrame([current_row] * genelist_len)
        t_DF = t_DF.reset_index()

        for j in range(0, genelist_len):
            t_DF.loc[j, 'myGeneName'] = genelist[j]

        tempDF = tempDF.append(t_DF, ignore_index=True)

    SV_df = SV_df.append(tempDF, ignore_index=True)
    SV_df = SV_df.drop('level_0', axis=1)
    SV_df = SV_df.reset_index()
    SV_df = SV_df.drop('level_0', axis=1)

    SV_df['mutation_Name'] = SV_df['myGeneName'].map(str) + ':' + SV_df['mutation_Name'].map(str)
    maintable_writer = pd.ExcelWriter('SVs_unfiltered.xlsx', engine='xlsxwriter')
    SV_df.to_excel(maintable_writer)
    maintable_writer.close()

    SV_df = SV_df[SV_df['myGeneName'] != '']

    print("SV_df.shape: ", SV_df.shape)

    SV_df['SampleName'] = SV_df.Sample.str[10:]

    maintable_writer = pd.ExcelWriter('SVs_unfiltered_expanded.xlsx', engine='xlsxwriter')
    SV_df.to_excel(maintable_writer)
    maintable_writer.close()


def process_fusions(fusions):

    fusion_df = fusions
    fusion_df['mutation_Name'] = 'Fusion_' + fusion_df['Gene1'].map(str) + ':' + fusion_df['Gene2'].map(str) + ':' + \
                                 fusion_df['Breakpoint1'].map(str) + fusion_df['Breakpoint2'].map(str) + '_' + \
                                 fusion_df['Fusion_predicted_effect'].map(str)

    maintable_writer = pd.ExcelWriter('fusions_unfiltered.xlsx', engine='xlsxwriter')
    fusion_df.to_excel(maintable_writer)
    print("fusion_df.shape: ", fusion_df.shape)
    maintable_writer.close()

    # filter fusions
    filter_terms = set(['banned', 'bodymap2', 'cacg', 'cta_gene', 'ctb_gene', 'ctc_gene', 'ctd_gene', 'distance1000bp',
                    'distance100kbp', 'distance10kbp', 'duplicates', 'ensembl_fully_overlapping',
                    'ensembl_partially_overlapping', 'ensembl_same_strand_overlapping', 'gtex', 'hpa', 'mt', 'paralogs',
                    'readthrough', 'refseq_fully_overlapping', 'refseq_partially_overlapping',
                    'refseq_same_strand_overlapping', 'rp11_gene', 'rp_gene', 'rrna', 'short_distance',
                    'similar_reads', 'similar_symbols', 'ucsc_fully_overlapping', 'ucsc_partially_overlapping',
                    'ucsc_same_strand_overlapping'])

    fusion_df = fusion_df.dropna(subset=['Description'])

    fusion_df['mark4deletion'] = fusion_df['Description'].apply(lambda
                                                                    x:len(set(x.split(',')).intersection(filter_terms)))
    fusion_df = fusion_df[fusion_df.mark4deletion == 0]
    fusion_df = fusion_df.drop('mark4deletion', axis=1)

    fusion_df.reset_index(inplace=True)

    maintable_writer = pd.ExcelWriter('fusions_NYGCfilter.xlsx', engine='xlsxwriter')
    fusion_df.to_excel(maintable_writer)
    print("fusion_df.shape post NYGC filter: ", fusion_df.shape)
    maintable_writer.close()

    fusion_df['myGeneName'] = ''

    numrows = fusion_df.shape[0]
    print("numrows: ", numrows)

    tempDF = pd.DataFrame()

    for i in range(0, numrows):
        current_row = fusion_df.loc[i]

        t_DF = pd.DataFrame([current_row] * 2)
        t_DF = t_DF.reset_index()

        t_DF.loc[0, 'myGeneName'] = current_row['Gene1']
        t_DF.loc[1, 'myGeneName'] = current_row['Gene2']

        tempDF = tempDF.append(t_DF, ignore_index=True)

    fusion_df = fusion_df.append(tempDF, ignore_index=True)
    fusion_df = fusion_df[fusion_df['myGeneName'] != '']
    fusion_df = fusion_df[fusion_df['Fusion_predicted_effect'] == 'in-frame']
    fusion_df = fusion_df[fusion_df['Spanning_pairs'] >= 5]

    print("fusion_df.shape post expansion: ", fusion_df.shape)

    fusion_df['SampleName'] = fusion_df.Sample.str[10:]

    maintable_writer = pd.ExcelWriter('fusions_MSKfiltered.xlsx', engine='xlsxwriter')
    fusion_df.to_excel(maintable_writer)
    print("fusion_df.shape post MSK filter: ", fusion_df.shape)
    maintable_writer.close()


def generate_mutation_table():

    doSNVs = 0
    doCNVs = 0
    doSVs = 0
    doIndels = 1
    doFusions = 1

    do_pre_treatement = 1
    do_normal = 0
    do_post_treatment = 1

    post_treatment_samples = ['TARGET-21-PAMYMA-41A-01D', 'TARGET-21-PAMYMA-41A-01R', 'TARGET-21-PANZLR-41A-02D',
                              'TARGET-21-PANZLR-41A-02R', 'TARGET-21-PARBTV-41A-02D', 'TARGET-21-PARBTV-41A-02R',
                              'TARGET-21-PARHRS-41A-01D', 'TARGET-21-PARLSL-41A-02D', 'TARGET-21-PARLSL-41A-02R',
                              'TARGET-21-PARNAW-41A-01D', 'TARGET-21-PARNAW-41A-01R', 'TARGET-21-PARXYR-41A-02D',
                              'TARGET-21-PARXYR-41A-02R', 'TARGET-21-PARZIA-41A-02D', 'TARGET-21-PARZIA-41A-02R',
                              'TARGET-21-PASDKZ-41A-02D', 'TARGET-21-PASDKZ-41A-03R', 'TARGET-21-PASDXR-41A-03D',
                              'TARGET-21-PASDXR-41A-03R', 'TARGET-21-PASFHK-41A-01D', 'TARGET-21-PASFHK-41A-01R',
                              'TARGET-21-PASFJJ-41A-01D', 'TARGET-21-PASFJJ-41A-01R', 'TARGET-21-PASFLG-41A-02D',
                              'TARGET-21-PASFLG-41A-02R', 'TARGET-21-PASFLG-41A-03R', 'TARGET-21-PASIGA-41A-01D',
                              'TARGET-21-PASIGA-41A-01R', 'TARGET-21-PASLZE-41A-01D', 'TARGET-21-PASLZE-41A-02R',
                              'TARGET-21-PASNKZ-41A-02D', 'TARGET-21-PASNKZ-41A-02R', 'TARGET-21-PASSLT-41A-01D',
                              'TARGET-21-PASSLT-41A-02R', 'TARGET-21-PASTZK-41A-03D', 'TARGET-21-PASTZK-41A-03R',
                              'TARGET-21-PASVJS-41A-02D', 'TARGET-21-PASVJS-41A-02R', 'TARGET-21-PASYEJ-41A-02D',
                              'TARGET-21-PASYEJ-41A-03R', 'TARGET-21-PASYWA-41A-01D', 'TARGET-21-PASYWA-41A-01R',
                              'TARGET-21-PATHIU-41A-01D', 'TARGET-21-PATISD-41A-01D', 'TARGET-21-PATISD-41A-01R',
                              'TARGET-21-PATJMY-41A-01D', 'TARGET-21-PATJMY-41A-01R', 'TARGET-21-PATKBK-41A-02D',
                              'TARGET-21-PATKBK-41A-02R', 'TARGET-21-PATKKJ-41A-01D', 'TARGET-21-PATKKJ-41A-01R',
                              'TARGET-21-PATKWH-41A-02D', 'TARGET-21-PATKWH-41A-02R', 'TARGET-21-PATAIJ-42A-01D',
                              'TARGET-21-PATAIJ-42A-01R']
    
    normal_samples = ['TARGET-21-PAMXZY-15A-01D', 'TARGET-21-PAMYMA-15A-01D', 'TARGET-21-PANVPB-15A-01D',
                      'TARGET-21-PANZLR-15A-01D', 'TARGET-21-PARBTV-15A-01D', 'TARGET-21-PARHRS-15A-01D',
                      'TARGET-21-PARLSL-15A-01D', 'TARGET-21-PARNAW-15A-01D', 'TARGET-21-PARXYR-15A-01D',
                      'TARGET-21-PARZIA-15A-01D', 'TARGET-21-PASDKZ-15A-01D', 'TARGET-21-PASDXR-15A-01D',
                      'TARGET-21-PASFHK-15A-01D', 'TARGET-21-PASFJJ-15A-01D', 'TARGET-21-PASFLG-15A-01D',
                      'TARGET-21-PASIGA-15A-01D', 'TARGET-21-PASLZE-15A-01D', 'TARGET-21-PASNKZ-15A-01D',
                      'TARGET-21-PASSLT-15A-01D', 'TARGET-21-PASTZK-15A-01D', 'TARGET-21-PASVJS-15A-01D',
                      'TARGET-21-PASYEJ-15A-01D', 'TARGET-21-PASYWA-15A-01D', 'TARGET-21-PATAIJ-15A-01D',
                      'TARGET-21-PATHIU-15A-01D', 'TARGET-21-PATISD-15A-01D', 'TARGET-21-PATJMY-15A-01D',
                      'TARGET-21-PATKBK-15A-01D', 'TARGET-21-PATKKJ-15A-02D', 'TARGET-21-PATKWH-15A-01D']
    
    pre_treatment_samples = ['TARGET-21-PAMXZY-09A-03D', 'TARGET-21-PAMXZY-09A-03R', 'TARGET-21-PAMYMA-09A-02D',
                             'TARGET-21-PAMYMA-09A-02R', 'TARGET-21-PANVPB-09A-02D', 'TARGET-21-PANVPB-09A-02R',
                             'TARGET-21-PANZLR-09A-03D', 'TARGET-21-PANZLR-09A-03R', 'TARGET-21-PARBTV-09A-03D',
                             'TARGET-21-PARBTV-09A-03R', 'TARGET-21-PARHRS-09A-01D', 'TARGET-21-PARHRS-09A-02R',
                             'TARGET-21-PARLSL-09A-02D', 'TARGET-21-PARLSL-09A-02R', 'TARGET-21-PARNAW-09A-01D',
                             'TARGET-21-PARNAW-09A-01R', 'TARGET-21-PARZIA-09A-01D', 'TARGET-21-PARZIA-09A-01R',
                             'TARGET-21-PASDXR-09A-02D', 'TARGET-21-PASDXR-09A-02R', 'TARGET-21-PASFHK-09A-01D',
                             'TARGET-21-PASFHK-09A-01R', 'TARGET-21-PASFJJ-09A-01D', 'TARGET-21-PASFJJ-09A-01R',
                             'TARGET-21-PASFLG-09A-04D', 'TARGET-21-PASIGA-09A-01D', 'TARGET-21-PASIGA-09A-01R',
                             'TARGET-21-PASLZE-09A-01D', 'TARGET-21-PASLZE-09A-01R', 'TARGET-21-PASNKZ-09A-01D',
                             'TARGET-21-PASNKZ-09A-01R', 'TARGET-21-PASSLT-09A-01D', 'TARGET-21-PASSLT-09A-01R',
                             'TARGET-21-PASTZK-09A-01D', 'TARGET-21-PASTZK-09A-01R', 'TARGET-21-PASVJS-09A-01D',
                             'TARGET-21-PASVJS-09A-01R', 'TARGET-21-PASYEJ-09A-01D', 'TARGET-21-PASYEJ-09A-01R',
                             'TARGET-21-PASYWA-09A-01D', 'TARGET-21-PASYWA-09A-01R', 'TARGET-21-PATAIJ-09A-01D',
                             'TARGET-21-PATAIJ-09A-01R', 'TARGET-21-PATHIU-09A-01D', 'TARGET-21-PATHIU-09A-01R',
                             'TARGET-21-PATISD-09A-01D', 'TARGET-21-PATISD-09A-01R', 'TARGET-21-PATJMY-09A-01D',
                             'TARGET-21-PATJMY-09A-01R', 'TARGET-21-PATKKJ-09A-01D', 'TARGET-21-PATKKJ-09A-01R',
                             'TARGET-21-PATKWH-09A-01D', 'TARGET-21-PATKWH-09A-01R', 'TARGET-21-PARXYR-03A-03D',
                             'TARGET-21-PARXYR-03A-03R', 'TARGET-21-PASDKZ-03A-02D', 'TARGET-21-PASDKZ-03A-02R',
                             'TARGET-21-PATKBK-03A-01D', 'TARGET-21-PATKBK-03A-01R']

    samples_to_keep = []
    if do_pre_treatement:
        samples_to_keep.append(pre_treatment_samples)
    if do_normal:
        samples_to_keep.append(normal_samples)
    if do_post_treatment:
        samples_to_keep.append(post_treatment_samples)

    samples_to_keep = flatten(samples_to_keep)
    if doSNVs:
        reportlist = glob.glob('SNVs/*.txt')
        SNVfiles_to_process = []
        for SNVfile in reportlist:
            # get sample name from File name
            SNVsample = calc_sample_name_from_filename(SNVfile)

            if SNVsample in samples_to_keep:
                SNVfiles_to_process.append(SNVfile)

        print("SNVs: ", SNVfiles_to_process)
        SNVs = getSNVs(SNVfiles_to_process)
        process_SNVs(SNVs)

    if doCNVs:
        reportlist = glob.glob('CNVs/*.bed')
        CNVfiles_to_process = []
        for CNVfile in reportlist:
            # get sample name from File name
            CNVsample = calc_sample_name_from_filename(CNVfile)

            if CNVsample in samples_to_keep:
                CNVfiles_to_process.append(CNVfile)

        print("CNVs: ", CNVfiles_to_process)
        CNVs = getCNVs(CNVfiles_to_process)
        process_CNV(CNVs)
    
    if doSVs:
        reportlist = glob.glob('SVs/*.bedpe')
        SVfiles_to_process = []
        for SVfile in reportlist:
            # get sample name from File name
            SVsample = calc_sample_name_from_filename(SVfile)

            if SVsample in samples_to_keep:
                SVfiles_to_process.append(SVfile)

        print("SVs: ", SVfiles_to_process)
        SVs = getSVs(SVfiles_to_process)
        process_SVs(SVs)
    
    if doIndels:
        reportlist = glob.glob('indels/*.txt')
        indelfiles_to_process = []
        for indelfile in reportlist:
            # get sample name from File name
            indelsample = calc_sample_name_from_filename(indelfile)

            if indelsample in samples_to_keep:
                indelfiles_to_process.append(indelfile)

        print("indels: ", indelfiles_to_process)
        indels = getINDELS(indelfiles_to_process)
        process_indels(indels)

    if doFusions:
        print('fusion')
        reportlist = glob.glob('./fusions/*.txt')
        fusions = getFusions(reportlist)

        # This is the file that Nico sent that has supported fusions. We can add it to the list
        fusionfile2 = pd.read_excel('fusions_unfiltered.subset.SVannotated.xlsx')
        print("fusionfile2.shape: ", fusionfile2.shape)
        fusions = fusions.append(fusionfile2, ignore_index=True)

        fusions['toKeep'] = fusions.apply(lambda row:1 if (row.Sample in samples_to_keep) else 0, axis=1)
        fusions = fusions[fusions['toKeep'] == 1]
        fusions.drop('toKeep', axis=1, inplace=True)
        process_fusions(fusions)

    SNVfile = pd.read_excel('SNVs_filtered.xlsx')
    print("SNVfile.shape: ", SNVfile.shape)

    SVfile = pd.read_excel('SVs_unfiltered_expanded.xlsx')
    print("SVfile.shape: ", SVfile.shape)

    CNVfile = pd.read_excel('CNVs_filtered.xlsx')
    print("CNVfile.shape: ", CNVfile.shape)

    indelfile = pd.read_excel('indel_unfiltered.xlsx')
    print("indelfile.shape: ", indelfile.shape)

    fusionfile = pd.read_excel('fusions_MSKfiltered.xlsx')
    print("fusionfile.shape: ", fusionfile.shape)

    mergefile = pd.DataFrame()
    mergefile = mergefile.append(SNVfile, ignore_index=True)
    mergefile = mergefile.append(SVfile, ignore_index=True)
    mergefile = mergefile.append(CNVfile, ignore_index=True)
    mergefile = mergefile.append(indelfile, ignore_index=True)
    mergefile = mergefile.append(fusionfile, ignore_index=True)
    mergefile['Patient'] = mergefile.SampleName.str[:6]

    census_df = pd.read_csv('cancer_gene_census.v84.txt', header=None)
    print("census_df.shape: ", census_df.shape)
    census = census_df[0].tolist()
    # print(census)

    mergefile['census'] = mergefile.apply(lambda row: 1 if (row.myGeneName in census) else 0, axis=1)

    ig_genes = ['IGFN1', 'IGHA1', 'IGHA2', 'IGHD', 'IGHD1-1', 'IGHD1-14', 'IGHD1-20', 'IGHD1-26', 'IGHD1-7', 'IGHD2-15',
                'IGHD2-2', 'IGHD2-21', 'IGHD2-8', 'IGHD3-10', 'IGHD3-16', 'IGHD3-22', 'IGHD3-3', 'IGHD3-9', 'IGHD4-11',
                'IGHD4-17', 'IGHD4-23', 'IGHD4-4', 'IGHD5-12', 'IGHD5-18', 'IGHD5-24', 'IGHD5-5', 'IGHD6-13',
                'IGHD6-19', 'IGHD6-25', 'IGHD6-6', 'IGHD7-27', 'IGHE', 'IGHEP1', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4',
                'IGHGP', 'IGHJ1', 'IGHJ1P', 'IGHJ2', 'IGHJ2P', 'IGHJ3', 'IGHJ3P', 'IGHJ4', 'IGHJ5', 'IGHJ6', 'IGHM',
                'IGHMBP2', 'IGHV1-12', 'IGHV1-14', 'IGHV1-17', 'IGHV1-18', 'IGHV1-2', 'IGHV1-24', 'IGHV1-3', 'IGHV1-45',
                'IGHV1-46', 'IGHV1-58', 'IGHV1-67', 'IGHV1-68', 'IGHV1-69', 'IGHV1-8', 'IGHV2-10', 'IGHV2-26',
                'IGHV2-5', 'IGHV2-70', 'IGHV3-11', 'IGHV3-13', 'IGHV3-15', 'IGHV3-16', 'IGHV3-19', 'IGHV3-20',
                'IGHV3-21', 'IGHV3-22', 'IGHV3-23', 'IGHV3-25', 'IGHV3-29', 'IGHV3-30', 'IGHV3-30-2', 'IGHV3-32',
                'IGHV3-33', 'IGHV3-33-2', 'IGHV3-35', 'IGHV3-36', 'IGHV3-37', 'IGHV3-38', 'IGHV3-41', 'IGHV3-42',
                'IGHV3-43', 'IGHV3-47', 'IGHV3-48', 'IGHV3-49', 'IGHV3-50', 'IGHV3-52', 'IGHV3-53', 'IGHV3-54',
                'IGHV3-57', 'IGHV3-6', 'IGHV3-60', 'IGHV3-62', 'IGHV3-63', 'IGHV3-64', 'IGHV3-65', 'IGHV3-66',
                'IGHV3-7', 'IGHV3-71', 'IGHV3-72', 'IGHV3-73', 'IGHV3-74', 'IGHV3-75', 'IGHV3-76', 'IGHV3-79',
                'IGHV3-9', 'IGHV3OR16-8', 'IGHV4-28', 'IGHV4-31', 'IGHV4-34', 'IGHV4-39', 'IGHV4-4', 'IGHV4-55',
                'IGHV4-59', 'IGHV4-61', 'IGHV4-80', 'IGHV5-51', 'IGHV5-78', 'IGHV6-1', 'IGHV7-27', 'IGHV7-34-1',
                'IGHV7-40', 'IGHV7-56', 'IGHV7-81', 'IGHVII-1-1', 'IGHVII-15-1', 'IGHVII-20-1', 'IGHVII-22-1',
                'IGHVII-26-2', 'IGHVII-28-1', 'IGHVII-30-1', 'IGHVII-31-1', 'IGHVII-33-1', 'IGHVII-40-1',
                'IGHVII-43-1', 'IGHVII-44-2', 'IGHVII-46-1', 'IGHVII-49-1', 'IGHVII-51-2', 'IGHVII-53-1', 'IGHVII-60-1',
                'IGHVII-62-1', 'IGHVII-65-1', 'IGHVII-67-1', 'IGHVII-74-1', 'IGHVII-78-1', 'IGHVIII-11-1',
                'IGHVIII-13-1', 'IGHVIII-16-1', 'IGHVIII-2-1', 'IGHVIII-22-2', 'IGHVIII-25-1', 'IGHVIII-26-1',
                'IGHVIII-38-1', 'IGHVIII-44', 'IGHVIII-47-1', 'IGHVIII-5-1', 'IGHVIII-5-2', 'IGHVIII-51-1',
                'IGHVIII-67-2', 'IGHVIII-67-3', 'IGHVIII-67-4', 'IGHVIII-76-1', 'IGHVIII-82', 'IGHVIV-44-1',
                'IGKC', 'IGKJ1', 'IGKJ2', 'IGKJ3', 'IGKJ4', 'IGKJ5', 'IGKV1-12', 'IGKV1-13', 'IGKV1-16', 'IGKV1-17',
                'IGKV1-22', 'IGKV1-27', 'IGKV1-32', 'IGKV1-33', 'IGKV1-35', 'IGKV1-37', 'IGKV1-39', 'IGKV1-5',
                'IGKV1-6', 'IGKV1-8', 'IGKV1-9', 'IGKV1D-12', 'IGKV1D-13', 'IGKV1D-16', 'IGKV1D-17', 'IGKV1D-22',
                'IGKV1D-27', 'IGKV1D-32', 'IGKV1D-33', 'IGKV1D-35', 'IGKV1D-37', 'IGKV1D-39', 'IGKV1D-42', 'IGKV1D-43',
                'IGKV1D-8', 'IGKV2-10', 'IGKV2-14', 'IGKV2-18', 'IGKV2-19', 'IGKV2-23', 'IGKV2-24', 'IGKV2-26',
                'IGKV2-28', 'IGKV2-29', 'IGKV2-30', 'IGKV2-36', 'IGKV2-38', 'IGKV2-4', 'IGKV2-40', 'IGKV2D-10',
                'IGKV2D-14', 'IGKV2D-18', 'IGKV2D-19', 'IGKV2D-23', 'IGKV2D-24', 'IGKV2D-26', 'IGKV2D-28', 'IGKV2D-29',
                'IGKV2D-30', 'IGKV2D-36', 'IGKV2D-38', 'IGKV2D-40', 'IGKV3-11', 'IGKV3-15', 'IGKV3-20', 'IGKV3-25',
                'IGKV3-31', 'IGKV3-34', 'IGKV3-7', 'IGKV3D-11', 'IGKV3D-15', 'IGKV3D-20', 'IGKV3D-25', 'IGKV3D-31',
                'IGKV3D-34', 'IGKV3D-7', 'IGKV4-1', 'IGKV5-2', 'IGKV6-21', 'IGKV6D-21', 'IGKV6D-41', 'IGKV7-3', 'IGLC1',
                'IGLC2', 'IGLC3', 'IGLC4', 'IGLC5', 'IGLC6', 'IGLC7', 'IGLJ1', 'IGLJ2', 'IGLJ3', 'IGLJ4', 'IGLJ5',
                'IGLJ6', 'IGLJ7', 'IGLL5', 'IGLV1-36', 'IGLV1-40', 'IGLV1-41', 'IGLV1-44', 'IGLV1-47', 'IGLV1-50',
                'IGLV1-51', 'IGLV1-62', 'IGLV10-54', 'IGLV10-67', 'IGLV11-55', 'IGLV2-11', 'IGLV2-14', 'IGLV2-18',
                'IGLV2-23', 'IGLV2-28', 'IGLV2-33', 'IGLV2-34', 'IGLV2-5', 'IGLV2-8', 'IGLV3-1', 'IGLV3-10', 'IGLV3-12',
                'IGLV3-13', 'IGLV3-15', 'IGLV3-16', 'IGLV3-17', 'IGLV3-19', 'IGLV3-2', 'IGLV3-21', 'IGLV3-22',
                'IGLV3-24', 'IGLV3-25', 'IGLV3-26', 'IGLV3-27', 'IGLV3-29', 'IGLV3-30', 'IGLV3-31', 'IGLV3-32',
                'IGLV3-4', 'IGLV3-6', 'IGLV3-7', 'IGLV3-9', 'IGLV4-3', 'IGLV4-60', 'IGLV4-69', 'IGLV5-37', 'IGLV5-45',
                'IGLV5-48', 'IGLV5-52', 'IGLV6-57', 'IGLV7-35', 'IGLV7-43', 'IGLV7-46', 'IGLV8-61', 'IGLV9-49',
                'IGLVI-20', 'IGLVI-38', 'IGLVI-42', 'IGLVI-56', 'IGLVI-63', 'IGLVI-68', 'IGLVI-70', 'IGLVIV-53',
                'IGLVIV-59', 'IGLVIV-64', 'IGLVIV-65', 'IGLVIV-66-1', 'IGLVV-58', 'IGLVV-66', 'IGLVVI-22-1',
                'IGLVVI-25-1', 'IGLVVII-41-1', 'IGSF1', 'IGSF3']

    mergefile['ig_locus'] = mergefile.apply(lambda row: 1 if (row.myGeneName in ig_genes) else 0, axis=1)

    tcr_genes = ['TRAJ17', 'TRAJ18', 'TRAJ19', 'TRAJ20', 'TRAJ21', 'TRAJ22', 'TRAJ23', 'TRAJ24', 'TRAJ25', 'TRAJ26',
                 'TRAJ27', 'TRAJ28', 'TRAJ29', 'TRAJ30', 'TRAJ31', 'TRAJ32', 'TRAJ33', 'TRAJ34', 'TRAJ35', 'TRAJ36',
                 'TRAJ37', 'TRAJ38', 'TRAJ39', 'TRAJ40', 'TRAJ41', 'TRAJ42', 'TRAJ43', 'TRAJ44', 'TRAJ45', 'TRAJ46',
                 'TRAJ47', 'TRAJ48', 'TRAJ49', 'TRAJ50', 'TRAJ51', 'TRAJ52', 'TRAJ53', 'TRAJ54', 'TRAJ55', 'TRAJ56',
                 'TRAJ57', 'TRAJ58', 'TRAJ59', 'TRAJ60', 'TRAJ61', 'TRAK1', 'TRAPPC8', 'TRAV12-1', 'TRAV12-2',
                 'TRAV12-3', 'TRAV13-1', 'TRAV13-2', 'TRAV14DV4', 'TRAV15', 'TRAV16', 'TRAV17', 'TRAV18', 'TRAV19',
                 'TRAV20', 'TRAV21', 'TRAV22', 'TRAV23DV6', 'TRAV24', 'TRAV25', 'TRAV26-1', 'TRAV26-2', 'TRAV27',
                 'TRAV28', 'TRAV29DV5', 'TRAV30', 'TRAV31', 'TRAV32', 'TRAV33', 'TRAV34', 'TRAV35', 'TRAV36DV7',
                 'TRAV37', 'TRAV38-1', 'TRAV38-2DV8', 'TRAV39', 'TRAV40', 'TRAV41', 'TRAV7', 'TRAV8-2', 'TRAV8-3',
                 'TRAV8-4', 'TRAV8-5', 'TRAV8-6', 'TRAV8-7', 'TRAV9-2', 'TRBC2', 'TRBJ2-1', 'TRBJ2-2', 'TRBJ2-2P',
                 'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-6', 'TRBJ2-7', 'TRBV1', 'TRBV10-1', 'TRBV10-2', 'TRBV11-1',
                 'TRBV11-2', 'TRBV12-1', 'TRBV12-2', 'TRBV19', 'TRBV2', 'TRBV20-1', 'TRBV21-1', 'TRBV22-1', 'TRBV23-1',
                 'TRBV24-1', 'TRBV25-1', 'TRBV26', 'TRBV27', 'TRBV28', 'TRBV29-1', 'TRBV3-1', 'TRBV30', 'TRBV4-1',
                 'TRBV4-2', 'TRBV5-1', 'TRBV5-2', 'TRBV5-3', 'TRBV5-4', 'TRBV5-5', 'TRBV5-6', 'TRBV5-7', 'TRBV6-1',
                 'TRBV6-4', 'TRBV6-5', 'TRBV6-6', 'TRBV6-7', 'TRBV6-8', 'TRBV6-9', 'TRBV7-1', 'TRBV7-3', 'TRBV7-4',
                 'TRBV7-5', 'TRBV7-6', 'TRBV7-7', 'TRBV7-8', 'TRBV8-1', 'TRBV8-2', 'TRBV9', 'TRBVA', 'TRBVB', 'TRDC',
                 'TRDD1', 'TRDD2', 'TRDD3', 'TRDJ1', 'TRDJ2', 'TRDJ3', 'TRDJ4', 'TRDV1', 'TRDV2', 'TRDV3', 'TRGC1',
                 'TRGC2', 'TRGJ2', 'TRGJP1', 'TRGJP2', 'TRGV1', 'TRGV11', 'TRGV4', 'TRGV8', 'TRGV9']

    mergefile['tcr'] = mergefile.apply(lambda row: 1 if (row.myGeneName in tcr_genes) else 0, axis=1)

    return mergefile


def generate_mutation_pivot(mergefile):
        # pivot table and draw plot
        merge_pivot = pd.pivot_table(mergefile, index='mutation_Name', columns='SampleName',
                                     values='myGeneName', aggfunc='count')

        return merge_pivot


def generate_gene_pivot(mergefile):
    # pivot table and draw plot
    merge_pivot = pd.pivot_table(mergefile, index='myGeneName', columns='SampleName',
                                 values='mutation_Name', aggfunc='count')

    return merge_pivot


def list_filtered():
    # read the source data
    SNVfile = pd.read_excel('SNVs_unfiltered.xlsx')
    print("SNVfile.shape from list_filtered: ", SNVfile.shape)

    SVfile = pd.read_excel('SVs_unfiltered_expanded.xlsx')
    print("SVfile.shape from list_filtered: ", SVfile.shape)

    CNVfile = pd.read_excel('CNVs_unfiltered.xlsx')
    print("CNVfile.shape from list_filtered: ", CNVfile.shape)

    indelfile = pd.read_excel('indel_unfiltered.xlsx')
    print("indelfile.shape from list_filtered: ", indelfile.shape)

    fusionfile = pd.read_excel('fusions_MSKfiltered.xlsx')
    print("fusionfile.shape from list_filtered: ", fusionfile.shape)

    #create the unfiltered merged data
    unfiltered_mergefile = pd.DataFrame()
    unfiltered_mergefile = unfiltered_mergefile.append(SNVfile, ignore_index=True)
    unfiltered_mergefile = unfiltered_mergefile.append(SVfile, ignore_index=True)
    unfiltered_mergefile = unfiltered_mergefile.append(CNVfile, ignore_index=True)
    unfiltered_mergefile = unfiltered_mergefile.append(indelfile, ignore_index=True)
    unfiltered_mergefile = unfiltered_mergefile.append(fusionfile, ignore_index=True)
    unfiltered_mergefile['Patient'] = unfiltered_mergefile.SampleName.str[:6]
    print("unfiltered_mergefile from list_filtered: ", unfiltered_mergefile.shape)
    unfiltered_mergefile.to_excel("unfiltered_mergefile.xlsx")

    # read the filtered merged data
    mergefile = pd.read_excel('merged.xlsx')
    print("mergefile.shape from list_filtered: ", mergefile.shape)

    # take the genelist from each file and then subtract the filtered genes from the unfiltered
    # this leaves behind a list of the genes that were removed by the filters
    filtered_genes = mergefile['myGeneName']
    unfiltered_genes = unfiltered_mergefile['myGeneName']
    filtered_genes = set(filtered_genes)
    unfiltered_genes = set(unfiltered_genes)
    removed_genes = list(unfiltered_genes - filtered_genes)
    print("removed genes from list_filtered: ", len(removed_genes))

    removed_genes_df = unfiltered_mergefile[unfiltered_mergefile['myGeneName'].isin(removed_genes)]
    print("removed_genes_df.shape from list_filtered: ", removed_genes_df.shape)
    removed_genes_df = removed_genes_df.drop('level_0',axis=1)
    removed_genes_df = removed_genes_df.reset_index()

    print("removed_genes_df.shape after reset_index: ", removed_genes_df.shape)
    removed_genes_df.to_excel('removed_genes_df_beforePivot.xlsx')
    removed_genes_pivot = pd.pivot_table(removed_genes_df, index=['myGeneName','mutation_Name'], columns='Patient',
                                 values='Sample', aggfunc=pd.Series.nunique)
    print("removed_genes_pivot.shape from list_filtered: ", removed_genes_pivot.shape)

    removed_genes_pivot['mysum'] = removed_genes_pivot[:].count(axis=1)
    # removed_genes_pivot = removed_genes_pivot.sort_values(by=['mysum'], ascending=False)
    removed_genes_pivot = removed_genes_pivot[removed_genes_pivot['mysum'] > 1]
    print("removed_genes_pivot.shape from list_filtered after removing singletons: ", removed_genes_pivot.shape)
    removed_genes_pivot = removed_genes_pivot.reset_index()

    census_df = pd.read_csv('cancer_gene_census.v84.txt', header=None)
    print("census_df.shape from list_filtered: ", census_df.shape)
    census = census_df[0].tolist()

    removed_genes_pivot['census'] = removed_genes_pivot.apply(lambda row: 1 if (row.myGeneName in census) else 0,
                                                              axis=1)
    maintable_writer = pd.ExcelWriter('removed_genes_pivot.xlsx', engine='xlsxwriter')
    removed_genes_pivot.to_excel(maintable_writer)
    maintable_writer.close()

    return removed_genes_pivot


def main():
    make_merge = 1

    if make_merge == 1:
        mergefile = generate_mutation_table()

        maintable_writer = pd.ExcelWriter('merged_pre_post.xlsx', engine='xlsxwriter')

        # clean up table
        # remove patients for whom we don't have all sequence samples
        filter_patients = set(['PAMXZY', 'PANVPB'])

        mergefile['mark4deletion'] = mergefile['Patient'].apply(lambda x:
                                                                len(set(x.split(',')).intersection(filter_patients)))
        mergefile = mergefile[mergefile.mark4deletion == 0]
        #mergefile = mergefile.drop('mark4deletion', axis=1)

        # remove samples that failed the QC review for CNV.
        filter_samples = set(['PASNKZ-09A-01D', 'PASFLG-09A-04D'])

        mergefile['deletion1'] = mergefile['SampleName'].apply(
            lambda x: len(set(x.split(',')).intersection(filter_samples)))

        filter_samples = set(['CNVs'])

        mergefile['deletion2'] = mergefile['mutation'].apply(
            lambda x: len(set(x.split(',')).intersection(filter_samples)))

        mergefile['deletion3'] = mergefile.apply(lambda row: 1 if (row.deletion1 and row.deletion2) else 0, axis=1)
        mergefile = mergefile[mergefile.deletion3 == 0]

        # remove samples that failed the QC review for SV.
        filter_samples = set(['PASNKZ-09A-01D'])

        mergefile['deletion1'] = mergefile['SampleName'].apply(
            lambda x: len(set(x.split(',')).intersection(filter_samples)))

        filter_samples = set(['SVs'])

        mergefile['deletion2'] = mergefile['mutation'].apply(
            lambda x: len(set(x.split(',')).intersection(filter_samples)))

        mergefile['deletion3'] = mergefile.apply(lambda row: 1 if (row.deletion1 and row.deletion2) else 0, axis=1)
        mergefile = mergefile[mergefile.deletion3 == 0]

        #drop all the genes that are a part of tcr_rearrangement and immuno_globulin
        mergefile = mergefile[mergefile.tcr == 0]
        mergefile = mergefile[mergefile.ig_locus == 0]
        mergefile = mergefile.drop(['deletion1', 'deletion2', 'deletion3', 'ig_locus', 'tcr',
                                    'mark4deletion'], axis=1)

        # mergefile = mergefile[(mergefile.mark4deletion == 0) and (mergefile.mutation == 'CNVs')]
        # mergefile = mergefile.drop([(mergefile.mark4deletion == 0) and (mergefile.mutation == 'CNVs')])

        # mark which samples are pre and post treatment
        mergefile['treatment'] = mergefile.apply(lambda row: 'pre' if ('-03A-' in row.SampleName) or
                                                                ('-09A-' in row.SampleName) else 'post', axis=1)

        mergefile['myGeneName'] = mergefile.apply(lambda row: (row.Gene1 + ":" + row.Gene2) if ('fusions' in row.mutation) else row.myGeneName, axis=1)

        #remove rows with no gene name calculated
        mergefile['myGeneName'].replace('',np.nan, inplace=True)
        mergefile.dropna(subset=['myGeneName'], inplace=True)

        mergefile.to_excel(maintable_writer)
        print("mergefile.shape from make_merge=1: ", mergefile.shape)
        maintable_writer.close()
    else:
        mergefile = pd.read_excel('merged.xlsx')

    # do it on the mutation level
    print("mergefile.shape from mutation: ", mergefile.shape)

    merge_pivot = generate_mutation_pivot(mergefile)
    table_writer = pd.ExcelWriter('merge_mutation.xlsx', engine='xlsxwriter')
    merge_pivot.to_excel(table_writer)
    table_writer.close()

    merge_pivot['mysum'] = merge_pivot[:].sum(axis=1)
    merge_pivot = merge_pivot.sort_values(by=['mysum'], ascending=False)
    merge_pivot = merge_pivot[merge_pivot['mysum'] > 3]

    merge_pivot = merge_pivot.drop('mysum', axis=1)
    # pdb.set_trace()
    merge_pivot.fillna(0, inplace=True)
    table_writer = pd.ExcelWriter('merge_mutation_pivottable.xlsx', engine='xlsxwriter')
    merge_pivot.to_excel(table_writer)
    table_writer.close()

    mycol = list(merge_pivot.columns.values)
    mycol = [e[:6] for e in mycol]  # read first 6 characters since that is what represents each patient
    merge_pivot.columns = mycol
    merge_pivot = merge_pivot.groupby(by=merge_pivot.columns, axis=1).sum()

    table_writer = pd.ExcelWriter('merge_mutation_pivottable_patient.xlsx', engine='xlsxwriter')
    merge_pivot.to_excel(table_writer)
    table_writer.close()

    # do it on the gene level
    print("mergefile.shape from genes: ", mergefile.shape)
    merge_pivot = generate_gene_pivot(mergefile)
    merge_pivot['mysum'] = merge_pivot[:].sum(axis=1)
    merge_pivot = merge_pivot.sort_values(by=['mysum'], ascending=False)
    merge_pivot = merge_pivot[merge_pivot['mysum'] > 0]

    merge_pivot = merge_pivot.drop('mysum', axis=1)
    #pdb.set_trace()
    merge_pivot.fillna(0, inplace=True)
    table_writer = pd.ExcelWriter('merge_gene_pivottable.xlsx', engine='xlsxwriter')
    merge_pivot.to_excel(table_writer)
    table_writer.close()

    mycol = list(merge_pivot.columns.values)
    mycol = [e[:6] for e in mycol] # read first 6 characters since that is what represents each patient
    merge_pivot.columns = mycol
    merge_pivot = merge_pivot.groupby(by=merge_pivot.columns, axis=1).sum()
    table_writer = pd.ExcelWriter('merge_gene_pivottable_patient.xlsx', engine='xlsxwriter')
    merge_pivot.to_excel(table_writer)
    table_writer.close()
    return 1


if __name__ == "__main__":
    main()
    #list_filtered()

# after merge_pre_post.xls is created, make two files, one for pre and one for post.
# remove the duplicate entries from exploded fusion table.
# remove all the RPLP and RPS proteins fusions where it's with another Ribosomal protein
# remove indels and SNVs with < 8 t_alt_count
# delete Malat1 since it's not real
# The following fusions were removed since they were only in one sample and weren't found to have STAR-Fusion support.
# BAZ1B--BCL7B not found.
# BCL7B--BAZ1B not found.
# SPI1--AZU1 not found.
# SPI1--DGAT1 not found.
# TMEM41B--DENND5A not found.
