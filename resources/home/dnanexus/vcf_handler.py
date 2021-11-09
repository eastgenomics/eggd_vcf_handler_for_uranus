"""
Script takes mutect VCFs and:
    1. Create vcf for import into BSVI
    2. Create excel workbook with sheets for each panel

BSVI vcf requires multiallelic sites are decomposed and that all genotype
fields are simplified to 0/1 (from, for example, 0/0/1/0).

Excel workbook for scientists to look at variants have sheets, one for each
panel, with annotations and preformatted report text for Epic.

Inputs are vcf files after VEP annotation and filtering:
    - *_allgenesvep.vcf (specified with -a argument)
    - *_[panel-name]vep.vcf (any number, specified with -v argument)
        e.g. *_lymphoidvep.vcf *_myeloidvep.vcf

Outputs:
    - bsvi vcf
    - tsv with variants
    - excel workbook with variants, one panel per sheet
"""

import argparse
import io
from pathlib import Path
import re
import subprocess
import sys
import pandas as pd
import xlsxwriter


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
    # read in vcf
    process = subprocess.Popen(
        f"cat {input_vcf} ", shell=True, stdout=subprocess.PIPE
    )

    vcf_data = io.StringIO()

    vcf_header = []

    for line in process.stdout:
        line = line.decode()
        vcf_data.write(line)

        if line.startswith('#'):
            # dump out header to list to write back for output vcf
            vcf_header.append(line)

    vcf_data.seek(0)

    cols = [
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
        "FORMAT", "SAMPLE"
    ]

    # read vcf into df
    vcf_df = pd.read_csv(vcf_data, sep="\t", comment='#', names=cols)

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


def df_report_formatting(fname, vcf_df):
    """
    Formats df of vcf records for report with INFO column split out to
    individual columns. Expects min. 14 '|' separated fields in INFO
    column to split out.

    Args:
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

    # calc Prev_count
    vcf_df['Prev_Count'] = vcf_df['Prev_AC'] + '/' + vcf_df['Prev_NS']

    # get index of AF in format column, should all be same and have a
    # list with 1 value, used to get AF from the sample column
    af_index = list(set(
        vcf_df['FORMAT'].apply(lambda x: x.split(':').index('AF')).to_list()
    ))

    # sense check all AF at same index
    assert len(af_index) == 1, \
        'Error in FORMAT column, AF not all at same index.'

    # get AF values from sample column add to new AF column, convert to %
    af_index = af_index[0]
    af_values = vcf_df['SAMPLE'].apply(lambda x: x.split(':')[af_index]).apply(
        lambda x: format(float(x) * 100)).apply(
        lambda x: '{:.1f}'.format(float(x)))
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
            f"{'in ' + x['EXON'] if x['EXON'] else ''} \n"
            f"HGVSc: {x['HGVSc'] if x['HGVSc'] else 'None'} \n"
            f"HGVSp: {x['HGVSp'] if x['HGVSp'] else 'None'} \n"
            f"COSMIC ID: {x['COSMIC'] if x['COSMIC'] else 'None'} \n"
            f"dbSNP: {x['dbSNP'] if x['dbSNP'] else 'None'} \n"
            f"Allele Frequency (VAF): {x['Mutect2_AF%'] + '%' if x['Mutect2_AF%']else 'None'}"
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
        'HGVSp', 'Protein_ID', 'CONSEQ', 'Read_Depth', 'Mutect2_AF%', 'FILTER',
        'ClinVar', 'ClinVar_CLNSIG', 'ClinVar_CLNDN', 'COSMIC', 'dbSNP',
        'gnomAD_AF', 'CADD_PHRED', 'Prev_Count', 'Report_text'
    ]]

    return vcf_df


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

    for panel_name, vcf_df in vcfs_dict.items():
        # loop over panel dfs, write to sheet & apply formatting
        vcf_df.to_excel(writer, sheet_name=panel_name)
        # fun excel formatting
        worksheet = writer.sheets[panel_name]
        wrap_format = workbook.add_format({'text_wrap': True})
        worksheet.set_column(1, 21, 15)
        worksheet.set_column(22, 22, 70, wrap_format)
        worksheet.set_default_row(100)
        worksheet.set_row(0, 15)

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

    # modify genotype of bsvi vcf
    bsvi_vcf_df = all_genes_df.copy(deep=True)
    bsvi_vcf_df = mod_genotype(bsvi_vcf_df)

    # get filename, use name of allgenes vcf as prefix for all
    fname = str(Path(args.allgenes).name)

    # apply formatting to allgenes vcf df for tsv file
    all_genes_df = df_report_formatting(fname, all_genes_df)

    # apply formatting to each panel df for xlsx file
    for panel, vcf_df in vcfs_dict.items():
        if not vcf_df.empty:
            vcf_df = df_report_formatting(fname, vcf_df)
            vcfs_dict[panel] = vcf_df

    write_bsvi_vcf(fname, bsvi_vcf_df, all_genes_df_header)
    write_tsv(fname, all_genes_df)
    write_xlsx(fname, vcfs_dict)
