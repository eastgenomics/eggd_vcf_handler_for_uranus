"""
Script takes mutect VCF and and cgppindel vcf to:
    1. Filter and annotate both VCFs against defined transcripts
    2. Create vcf for import into BSVI
    3. Create excel workbook with sheets for each panel

BSVI vcf requires multiallelic sites are decomposed and that all genotype
fields are simplified to 0/1 (from, for example, 0/0/1/0).

Excel workbook for scientists to look at variants have sheets, one for each
panel, with annotations and preformatted report text for Epic.

Inputs are vcf files after VEP annotation and filtering:
    - *_allgenesvep.vcf (specified with -a argument)
    - *_[panel-name]vep.vcf (any number, specified with -v argument)
        e.g. *_lymphoidvep.vcf *_myeloidvep.vcf
    - pindel vcf (specified with -p)

Outputs:
    - annotated and filtered mutect2 and pindel vcf
    - annotated and filtered panel vcfs
    - bsvi vcf
    - tsv with variants
    - excel workbook with variants, one panel per sheet
"""

import argparse
from pathlib import Path
import re
import sys
import numpy as np
import pandas as pd


def parse_args():
    """
    Parse command line args

    Args: None

    Returns:
        - args (Namespace): object containing parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            'Generates a tsv of variants, a decomposed VCF for BSVI,'
            'and an excel workbook of variants for each panel'
        )
    )

    parser.add_argument(
        '-a', '--allgenes',
        help='allgenes VCF from which to generate tsv file and BSVI vcf'
    )

    parser.add_argument(
        '-v', '--vcfs', nargs='*',
        help='Panel filtered VCF(s) from which to generate excel workbook'
    )
    parser.add_argument(
        '-p', '--pindel',
        help='Output VCF from cgppindel'
    )

    args = parser.parse_args()

    return args


def read_vcf(input_vcf):
    """
    Reads vcf into pandas df, returns header as a list for output (bsvi) vcf

    Args:
        - input_vcf (file): vcf file to read in

    Returns:
        - vcf_df (df): df of variants from vcf
        - vcf_header (list): header from vcf, for writing output (bsvi) vcf
    """
    # read in vcf to get header lines
    with open(input_vcf, 'r') as fh:
        vcf_data = fh.readlines()

    vcf_header = []

    for line in vcf_data:
        if line.startswith('#'):
            # dump out header to list to write back for output vcf
            vcf_header.append(line)
        else:
            break

    if "TUMOUR" in vcf_header[-1]:
        # pindel vcf has NORMAL & TUMOUR instead of SAMPLE
        cols = [
            "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", "NORMAL", "TUMOUR"
        ]
    else:
        cols = [
            "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", "SAMPLE"
        ]

    # read vcf records into df
    vcf_df = pd.read_csv(
        input_vcf, sep="\t", comment='#', names=cols, compression='infer'
    )

    return vcf_df, vcf_header


def mod_genotype(vcf_df):
    """
    Overwrite GT with 0/1 when it's derived from a multiallelic
    e.g. 0/0/1/0 -> 0/1 and 0/0/0/1 -> 0/1

    Args:
        - input_vcf (file): vcf file passed at cmd line

    Returns:
        - vcf_df (df): df of variants with modified genotype
    """
    for index, row in vcf_df.iterrows():
        # loop over rows, change genotype if contains greater than 2 fields
        sample = row['SAMPLE'].split(':')

        # sense check that genotype field doesn't have anything funky,
        # if it does then it can be reviewed manually
        assert len(sample[0]) >= 3, \
            f'Genotype field has < 3 characters: {sample[0]} for sample: {row}'

        if len(sample[0]) > 3:
            # >3 => not 0/1 => modify
            sample[0] = '0/1'
            sample = ':'.join(sample)

            # write new entry back to row
            vcf_df.at[index, 'SAMPLE'] = sample

    return vcf_df


def get_field_index(column, field):
    """
    Returns the index of a given field in a column with values separated by ':'

    Args:
        - column (pd.Series): df column to get field index
        - field (str): field to get index of

    Returns: index of given field
    """
    return list(column.apply(lambda x: x.split(':').index(field)))[0]


def get_field_value(column, index):
    """
    Given a column with ':' separated values and index, return a
    series of values with the value split by the index

    Args:
        - column (pd.Series): column to split and return from
        - index (int): index to select

    Returns: series of values selected by index
    """
    return column.apply(lambda x: int(x.split(':')[index]))


def filter_common(vcf_df):
    """
    Filters for common variants by prev_count > 50% & synonymous variants
    EXCEPT in TP53 and GATA2

    Args: vcf_df (df): df of variants

    Returns:
        - vcf (df): df of unfiltered variants
        - filter_vcf (df): df of filtered out common variants
    """
    # prev_count formatted as prev_ac/prev_samples (i.e. 45/205)
    filter_idxs = np.where((
        vcf_df['Prev_Count'].apply(
            lambda x: (int(x.split("/")[0]) / int(x.split("/")[1])) > 0.5
        ) | (
            vcf_df['CONSEQ'] == 'synonymous_variant'
        )) & (
            np.logical_not(vcf_df['GENE'].isin(['TP53', 'GATA2']))
    ))

    # filter df by indixes of filter conditions
    filtered_df = vcf_df.loc[filter_idxs]
    vcf_df = vcf_df.drop(filter_idxs[0])

    return vcf_df, filtered_df


def df_report_formatting(panel, vcf_df):
    """
    Formats df of vcf records for report with INFO column split out to
    individual columns. Expects min. 15 '|' separated fields in INFO
    column to split out.

    Args:
        - panel (str): panel name of vcf
        - vcf_df (df): df of variants from vcf

    Returns:
        - vcf_df (df): df of variants with fun formatting
    """
    # Get DP field from INFO column, use join instead of
    # using index in case of missing and being empty => index error
    vcf_df['Read_Depth'] = vcf_df['INFO'].str.split(';').apply(
        lambda x: ''.join((y for y in x if y.startswith('DP=')))
    )
    vcf_df['Read_Depth'] = vcf_df['Read_Depth'].apply(
        lambda x: x.replace('DP=', '')
    )

    # get CSQ field from INFO column, use join instead of
    # using index in case of missing and being empty => index error
    vcf_df['INFO'] = vcf_df['INFO'].str.split(';').apply(
        lambda x: ''.join((y for y in x if y.startswith('CSQ=')))
    )

    info_cols = [
        'GENE', 'VARIANT_CLASS', 'CONSEQ', 'EXON', 'HGVSc', 'HGVSp',
        'gnomAD_AF', 'CADD_PHRED', 'DB', 'ClinVar', 'ClinVar_CLNDN',
        'ClinVar_CLNSIG', 'COSMIC', 'Prev_AC', 'Prev_NS', 'Feature'
    ]

    # splits info column to cols defined in info_cols
    vcf_df[info_cols] = vcf_df['INFO'].str.split('|', -1, expand=True)

    # drop unndeeded Feature column (transcript field from vep filter)
    vcf_df.drop(columns=['Feature'], inplace=True)

    # remove info id from gene
    vcf_df['GENE'] = vcf_df['GENE'].apply(lambda x: x.replace('CSQ=', ''))

    # calculate Prev_count, first adjust those not previously seen that have
    # empty strings for prev_ac and prev_ns

    # first get total number of samples across all, should return single value
    uniq_prev_ns = list(set(filter(None, vcf_df['Prev_NS'])))

    assert len(uniq_prev_ns) == 1, \
        f"Differing total previous samples identified: {uniq_prev_ns}"

    uniq_prev_ns = uniq_prev_ns[0]

    # those not previously seen will have empty string, fill appropriately to
    # display as 0/{total}
    vcf_df['Prev_AC'] = vcf_df['Prev_AC'].apply(
        lambda x: str(0) if x == "" else x)
    vcf_df['Prev_NS'] = vcf_df['Prev_NS'].apply(
        lambda x: uniq_prev_ns if x == "" else x)

    vcf_df['Prev_Count'] = vcf_df['Prev_AC'] + '/' + vcf_df['Prev_NS']

    if panel == "pindel":
        # handle pindel vcf, AF to be calculated from TUMOUR field
        # this is calculated as (PU + NU) / (PR + NR)
        # values are described in table 15.7.3 here:
        # # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6097606/
        caller = "Pindel"

        # get indices of required fields
        pu_index = get_field_index(vcf_df['FORMAT'], 'PU')
        nu_index = get_field_index(vcf_df['FORMAT'], 'NU')
        pr_index = get_field_index(vcf_df['FORMAT'], 'PR')
        nr_index = get_field_index(vcf_df['FORMAT'], 'NR')

        # get values for given field from TUMOUR column
        pu_values = get_field_value(vcf_df['TUMOUR'], pu_index)
        nu_values = get_field_value(vcf_df['TUMOUR'], nu_index)
        pr_values = get_field_value(vcf_df['TUMOUR'], pr_index)
        nr_values = get_field_value(vcf_df['TUMOUR'], nr_index)

        # calculate af & format as pct to 1dp
        af_values = (pu_values + nu_values) / (pr_values + nr_values)
        af_pcts = af_values.apply(lambda x: '{:.1f}'.format(float(x * 100)))

        vcf_df.insert(16, 'Pindel_AF%', af_pcts)

        # depth is also not in the INFO field so will be blank for pindel,
        # calulating this as PR + NR
        vcf_df['Read_Depth'] = pr_values + nr_values
    else:
        # mutect2 vcf
        caller = "Mutect2"

        # get index of AF in format column, should all be same and have a
        # list with 1 value, used to get AF from the sample column
        af_index = list(set(vcf_df['FORMAT'].apply(
            lambda x: x.split(':').index('AF'))
        ))

        # sense check all AF at same index
        assert len(af_index) == 1, \
            'Error in FORMAT column, AF not all at same index.'

        # get AF values from sample column add to new AF column, convert to %
        af_index = af_index[0]
        af_values = vcf_df["SAMPLE"].apply(
            lambda x: x.split(':')[af_index]
        ).apply(
            lambda x: format(float(x) * 100)
        ).apply(
            lambda x: '{:.1f}'.format(float(x))
        )

        vcf_df.insert(16, 'Mutect2_AF%', af_values)

    # cosmic annotation returns duplicates for each record in cosmic vcf
    # turn to set to be unique, join in case there is more than one
    vcf_df['COSMIC'] = vcf_df['COSMIC'].apply(
        lambda x: (','.join(set(x.split('&')))) if x else x
    )

    # waiting to hear if the haemonc team want hgmd or not
    # vcf_df['HGMD'] = vcf_df['DB'].str.split(r'\&|\||,').apply(
    #     lambda x: ','.join((y for y in x if y.startswith('CM', 'CD')))
    # )
    vcf_df['dbSNP'] = vcf_df['DB'].str.split(r'\&|\||,').apply(
        lambda x: ','.join((y for y in x if y.startswith('rs')))
    )

    # Include word exon in exon field to overcome excel
    vcf_df['EXON'] = vcf_df['EXON'].apply(lambda x: 'exon ' + x if x else None)

    # scientists are picky and want NM_ and NP_ changing in HGVS
    vcf_df['Transcript_ID'] = vcf_df['HGVSc'].str.split(':').apply(
        lambda x: ','.join((y for y in x if y.startswith('NM')))
    )
    vcf_df['HGVSc'] = vcf_df['HGVSc'].str.split(':').apply(
        lambda x: ','.join((y for y in x if y.startswith('c')))
    )
    vcf_df['Protein_ID'] = vcf_df['HGVSp'].str.split(':').apply(
        lambda x: ','.join((y for y in x if y.startswith('NP')))
    )
    vcf_df['HGVSp'] = vcf_df['HGVSp'].str.split(':').apply(
        lambda x: ','.join((y for y in x if y.startswith('p')))
    )
    regex = re.compile(r'^[p\.]*')
    vcf_df['HGVSp'] = vcf_df['HGVSp'].apply(
        lambda x: regex.sub('p.(', x) + ')' if x else None)

    # add interestingly formatted report text column
    vcf_df['Report_text'] = vcf_df[vcf_df.columns.tolist()].apply(
        lambda x: (
            f"{x['GENE']} {x['CONSEQ']} "
            f"{'in ' + x['EXON'].split('/')[0] if x['EXON'] else ''} \n"
            f"HGVSc: {x['HGVSc'] if x['HGVSc'] else 'None'} \n"
            f"HGVSp: {x['HGVSp'] if x['HGVSp'] else 'None'} \n"
            f"COSMIC ID: {x['COSMIC'] if x['COSMIC'] else 'None'} \n"
            f"dbSNP: {x['dbSNP'] if x['dbSNP'] else 'None'} \n"
            f"""Allele Frequency (VAF): {
                str(x[f'{caller}_AF%']) + '%' if x[f'{caller}_AF%'] else 'None'
            }"""
        ), axis=1
    )

    # get sample name from filename (should always be >>20 characters)
    if len(fname) > 19:
        samplename = fname[0:20] + '...'
    else:
        samplename = fname
    vcf_df['samplename'] = samplename

    # select and re-order df columns
    vcf_df = vcf_df[[
        'samplename', 'CHROM', 'POS', 'GENE', 'Transcript_ID', 'EXON', 'HGVSc',
        'HGVSp', 'Protein_ID', 'CONSEQ', 'Read_Depth', f'{caller}_AF%',
        'FILTER', 'ClinVar', 'ClinVar_CLNSIG', 'ClinVar_CLNDN', 'COSMIC',
        'dbSNP', 'gnomAD_AF', 'CADD_PHRED', 'Prev_Count', 'Report_text'
    ]]

    return vcf_df


def to_report_formatting(col_names):
    """
    Build df of required placeholder text for the to report sheet

    Args: col_names (list): list of column names from a panel df

    Returns: report_df (df): formatted dataframe wth placeholder text as report
        template in specific cells
    """
    # add the column names for when they paste in reported variants
    report_df = pd.DataFrame(columns=col_names).astype('object')

    # add 12 empty rows as padding for visualise niceness
    report_df = report_df.append(
        [pd.Series([np.nan])] * 12).reindex(col_names, axis=1)

    col1_labels = [
        "Run QC", "250x", "Contamination", "Total reads M", "Fold 80",
        "Insert Size",
    ]


    for label in col1_labels:
        # set each of the col1 labels as a field in first column
        report_df = report_df.append(
            {report_df.columns[0]: label}, ignore_index=True
        )

    report_df.at[12, report_df.columns[3]] = "Sample QC"

    col6_labels = ["Analysed by", "Date", "Subpanel analysed"]
    col6_label_idx = 8

    for label in col6_labels:
        # start at row 8, add each of col6 labels as field in column 6
        report_df.at[col6_label_idx, report_df.columns[5]] = label
        col6_label_idx += 1

    return report_df


def set_column_widths(worksheet, df):
    """
    Sets column widths dyanmically based on cell content

    Args:
        - worksheet (xlsxwriter sheet object): sheet to modify
        - df (pd.DataFrame): dataframe for sheet to set widths from
    Returns:
        - - worksheet (xlsxwriter sheet object): xlsx sheet with new widths
    """
    for idx, col in enumerate(df):
        # specific formatting for columns with potential very long text
        if "ClinVar" in col:
            max_len = 15
        elif "Report_text" in col:
            max_len = 40
        else:
            series = df[col]
            max_len = max((
                series.astype(str).map(len).max(), len(str(series.name))
            )) + 2
        worksheet.set_column(idx, idx, max_len)  # set column width

    return worksheet


def write_bsvi_vcf(fname, bsvi_df, bsvi_vcf_header):
    """
    Write df of variants with modified genotype to vcf for BSVI

    Args:
        - fname (str): name of allgenes input vcf, used to name output file
        - bsvi_df (df): df of variants with modified genotype for bsvi
        - bsvi_vcf_header (list): vcf header from original vcf

    Returns: None

    Outputs: .vcf file of variants for importing into bsvi
    """
    vcf_fname = fname.replace('allgenesvep', 'allgenes_bsvi')

    with open(vcf_fname, 'w') as f:
        for line in bsvi_vcf_header:
            # write header to vcf
            f.write(line)

        # apend variants to vcf
        bsvi_df.to_csv(f, sep='\t', header=False, index=False)


def write_tsv(fname, all_genes_df):
    """
    Write tsv file of all genes from allgenes vcf

    Args:
        - fname (str): name of allgenes input vcf, used to name output file
        - all_genes_df (df): df of all genes from input vcf

    Returns: None

    Outputs: .tsv file with all gene variants
    """
    tsv_fname = fname.replace('allgenesvep.vcf', 'allgenes.tsv')

    # write tsv file
    with open(tsv_fname, 'w') as tsv:
        all_genes_df.to_csv(tsv, sep='\t', header=True, index=False)


def write_xlsx(fname, vcfs_dict):
    """
    Write xlsx file of panel filtered vcfs, each panel goes to a separate tab

    Args:
        - fname (str): name of allgenes input vcf, used to name output file
        - vcfs_dict (dict): dict of panel vcfs

    Returns: None

    Outputs: .xlsx file of panel vcfs in separate tabs
    """
    excel_fname = fname.replace('allgenesvep.vcf', 'panels.xlsx')

    # write excel
    writer = pd.ExcelWriter(excel_fname, engine="xlsxwriter")
    workbook = writer.book

    # get the column names of a df to write to report sheet
    header = vcfs_dict[next(iter(vcfs_dict))].columns
    report_df = to_report_formatting(col_names=header)
    report_df.to_excel(writer, sheet_name='to report', index=False)

    worksheet = writer.sheets['to report']

    # set specific cells to be bold in to report tab
    for cell in ['A14:A14', 'D14:D14', 'F10:F12']:
        header_format = workbook.add_format({'bold': True})
        worksheet.conditional_format(
            cell, {'type': 'no_errors', 'format': header_format}
        )

    # set dynamic column widths on cell content
    worksheet = set_column_widths(worksheet, report_df)
    worksheet.set_row(0, 12)

    for panel_name, vcf_df in vcfs_dict.items():
        # loop over panel dfs, write to sheet & apply formatting
        vcf_df.to_excel(writer, sheet_name=panel_name, index=False)
        # fun excel formatting
        worksheet = writer.sheets[panel_name]
        wrap_format = workbook.add_format({'text_wrap': True})
        worksheet.set_column(1, 21, 15)
        worksheet.set_column(22, 22, 70, wrap_format)
        worksheet.set_default_row(100)
        worksheet.set_row(0, 12)

        # set dynamic column widths on cell content
        worksheet = set_column_widths(worksheet, vcf_df)

    writer.save()


if __name__ == "__main__":
    args = parse_args()

    # read allgenes vcf into df, retain header for output vcf
    all_genes_df, all_genes_df_header = read_vcf(args.allgenes)

    # build dict of panel vcfs, allows to store panel name with df
    vcfs_dict = {}

    for vcf in args.vcfs:
        # loop over panel vcfs, read into df and add to dict
        # get panel name from vcf name
        panel = Path(vcf).stem.split('_')[-1].rstrip('vep')
        panel_df, _ = read_vcf(vcf)  # don't retain vcf header as not needed

        vcfs_dict[panel] = panel_df

    # read in pindel vcf and add to vcfs_dict to be included in excel
    pindel_df, _ = read_vcf(Path(args.pindel))
    vcfs_dict["pindel"] = pindel_df

    # modify genotype of bsvi vcf
    bsvi_vcf_df = all_genes_df.copy(deep=True)
    bsvi_vcf_df = mod_genotype(bsvi_vcf_df)

    # get filename, use name of allgenes vcf as prefix for all
    fname = str(Path(args.allgenes).name)

    # apply formatting to allgenes vcf df for tsv file
    all_genes_df = df_report_formatting(fname, all_genes_df)

    formatted_dfs = {}

    # apply formatting to each panel df for xlsx file
    for panel, vcf_df in vcfs_dict.items():
        if not vcf_df.empty:
            vcf_df = df_report_formatting(panel, vcf_df)
            formatted_dfs[panel] = vcf_df
        if panel == 'myeloid':
            # filtering myeloid panel for common variants
            vcf_df, filter_vcf_df = filter_common(vcf_df)
            formatted_dfs['myeloid'] = vcf_df
            formatted_dfs['myeloid_filtered'] = filter_vcf_df

    write_bsvi_vcf(fname, bsvi_vcf_df, all_genes_df_header)
    write_tsv(fname, all_genes_df)
    write_xlsx(fname, formatted_dfs)
