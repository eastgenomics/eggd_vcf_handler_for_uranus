"""
Script to modify VCFs from mutect2 for importing to BSVI.
Multiallelic sites need splitting first with bcftools norm and then the
genotype fields are split from 0/0/1/0 -> 0/1 for BSVI to handle.

Outputs vcf with split multiallelics and a tsv with formatted INFO
column to individual columns.

Jethro Rainford
210311
"""
import io
from pathlib import Path
import re
import subprocess
import sys

import pandas as pd


def mod_genotype(input_vcf):
    """
    Overwrite GT with 0/1 when it's derived from a multiallelic
    e.g. 0/0/1/0 -> 0/1 and 0/0/0/1 -> 0/1

    Args:
        - input_vcf (file): vcf file passed at cmd line

    Returns:
        - vcf_header (list): header lines read in from VCF
        - vcf_df (df): df of variants
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
            # dump out header to list to write back
            vcf_header.append(line)

    vcf_data.seek(0)

    cols = [
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
        "FORMAT", "SAMPLE"
    ]

    # read normalised vcf into df
    vcf_df = pd.read_csv(vcf_data, sep="\t", comment='#', names=cols)

    print('Adjusting multiallelic genotypes')
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

    return vcf_header, vcf_df


def generate_tsv(tsv_df):
    """
    Generates tsv file from modified vcf with INFO column split out to
    individual columns. Expects min. 9 '|' separated fields in INFO
    column to split out.

    Args:
        - vcf_df (df): df of variants from vcf

    Returns:
        - tsv_df (df): df of variants with split info column
    """
    # sense check correct annotation has been added to all rows else it
    # gives an unhelpful pandas error on trying to split
    assert all(tsv_df.INFO.str.count(r'\|') > 8), \
        "Incorrectly formatted INFO field, some records have < 9 fields."

    # keep just the CSQ field from INFO column, use join instead of
    # using index in case of missing and being empty => index error
    tsv_df['INFO'] = tsv_df['INFO'].str.split(';').apply(
        lambda x: ''.join((y for y in x if y.startswith('CSQ=')))
    )

    info_cols = [
        'GENE', 'VARIANT_CLASS', 'CONS', 'EXON', 'HGVSc',
        'HGVSp', 'gnomAD_AF', 'SIFT', 'POLYPHEN', 'DB'
    ]

    # splits info column to cols defined in info_cols
    tsv_df[info_cols] = tsv_df['INFO'].str.split('|', 9, expand=True)

    # remove info id from gene
    tsv_df['GENE'] = tsv_df['GENE'].apply(lambda x: x.replace('CSQ=', ''))

    # get index of AF in format column, should all be same and have a
    # list with 1 value, used to get AF from the sample column
    af_index = list(set(
        tsv_df['FORMAT'].apply(lambda x: x.split(':').index('AF')).to_list()
    ))

    # sense check all AF at same index
    assert len(af_index) == 1, \
        'Error in FORMAT column, AF not all at same index.'

    # get AF values from sample column add to new AF column
    af_index = af_index[0]
    af_values = tsv_df['SAMPLE'].apply(lambda x: x.split(':')[af_index])
    tsv_df.insert(11, 'AF', af_values)

    # split messy DB annotation column out to clinvar, cosmic & dbsnp
    # cols have multiple fields and diff delimeters then join with ','
    # in case of having more than one entry
    tsv_df['COSMIC'] = tsv_df['DB'].str.split(r'\&|\||,').apply(
        lambda x: ','.join((y for y in x if y.startswith('COS')))
    )
    tsv_df['CLINVAR'] = tsv_df['DB'].str.split(r'\&|\||,').apply(
        lambda x: ','.join((y for y in x if y.startswith('CM')))
    )
    tsv_df['dbSNP'] = tsv_df['DB'].str.split(r'\&|\||,').apply(
        lambda x: ','.join((y for y in x if y.startswith('rs')))
    )

    # scientists are picky and want NM_ and NP_ changing in HGVS
    regex = re.compile(r'^[A-Z]*_[0-9]*.[0-9]*:')
    tsv_df['HGVScshort'] = tsv_df['HGVSc'].apply(lambda x: regex.sub('', x))
    tsv_df['HGVSpshort'] = tsv_df['HGVSp'].apply(lambda x: regex.sub('', x))

    # add interestingly formatted report text column
    tsv_df['Report_text'] = tsv_df[tsv_df.columns.tolist()].apply(
        lambda x: (
            f"{x['GENE']} {x['VARIANT_CLASS']} variant "
            f"{'in exon ' + x['EXON'] if x['EXON'] else ''} \r\n"
            f"HGVSc: {x['HGVScshort'] if x['HGVScshort'] else 'None'} \r\n"
            f"HGVSp: {x['HGVSpshort'] if x['HGVSpshort'] else 'None'} \r\n"
            f"COSMIC ID : {x['COSMIC'] if x['COSMIC'] else 'None'} \r\n"
            f"Allele Frequency (VAF): {x['AF'] if x['AF'] else 'None'}"
        ), axis=1
    )

    # drop unneeded columns
    tsv_df = tsv_df.drop(['INFO', 'DB', 'HGVScshort', 'HGVSpshort'], axis=1)

    return tsv_df


def write_files(input_vcf, vcf_header, vcf_df, tsv_df):
    """
    Write modified vcf and tsv df with split info field to files

    Args:
        - input_vcf (file): vcf file passed at cmd line
        - vcf_header (list): header lines read in from VCF
        - vcf_df (df): df of variants
        - tsv_df (df): df of variants with split info column

    Outputs:
        - vcf file with modified multiallelic records
        - tsv file with modified multialleic records and split info field
    """
    vcf_fname = str(Path(input_vcf).name).replace('vepfilter', 'bsvi')
    tsv_fname = vcf_fname.replace('bsvi.vcf', 'variantlist.tsv')

    print(f'Writing vcf to outfile: {vcf_fname}')
    print(f'Writing tsv to outfile: {tsv_fname}')

    with open(vcf_fname, 'w') as f:
        for line in vcf_header:
            # write header to vcf
            f.write(line)

    # apend variants to vcf
    with open(vcf_fname, 'a') as f:
        vcf_df.to_csv(f, sep='\t', header=False, index=False)

    # write tsv file
    with open(tsv_fname, 'w') as tsv:
        tsv_df.to_csv(tsv, sep='\t', header=True, index=False)


if __name__ == "__main__":

    # check only one vcf passed
    assert len(sys.argv) == 2, 'Incorrect no. VCFs passed, requires one.'

    input_vcf = sys.argv[1]
    vcf_header, vcf_df = mod_genotype(input_vcf)

    # copy variant df and modify appropriately for tsv
    tsv_df = vcf_df.copy(deep=True)
    tsv_df = generate_tsv(tsv_df)

    write_files(input_vcf, vcf_header, vcf_df, tsv_df)
