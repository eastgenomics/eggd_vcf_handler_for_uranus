# eggd vcf handler for uranus

## What does this app do?
### This app uses bedtools, bcftools and VEP to take VCFs from sentieon mutect2 and pindel:
- Annotates and filters sentieon mutect2 and pindel VCFs
- Produces a excel workbook with sub-panels presented in separate sheets; also provides preformatted text to aid Epic data entry
- Produces a VCF that can be used as input for BSVI, as BSVI can't handle mutect2's (VCFv4.2 compliant) representation of multiallelic variants in VCF

## What are typical use cases for this app?
### This app uses the follow tools which are app assets:
* bcftools (v1.12)
* bedtools (v2.30.0)
* python_packages (numpy-1.20.1, pandas-1.2.3, pytz-2021.1, XlsxWriter-1.4.0)

### This app uses the following provided as inputs:

* VEP (v103.1) (docker image)
* VEP refseq (v103) annotation sources
* CADD (v1.6) which now includes splicing
* ClinVar VCF (20210501 release) modified to add chr prefix
* merged VCF containing counts of each variant detected known set of samples (provided as separate input from tarball)
    

### This app has access to the Internet

## What data are required for this app to run?
- VCF output from sentieon mutect2 as part of Uranus workflow
- VCF output (+ index) from cgppindel as part of Uranus workflow
- BED file that details ROIs for myeloid NGS panel (default specified)
- BED file that details ROIs for filtering cgppindel output (default specified)
- Genome FASTA and index that was used by to generate the VCF (default specified)
- VEP docker image (`vep_docker`, default specified)
- VEP plugins (`vep_plugins`, defaults specified)
- VEP reference files (`vep_refs`, defaults specified)
- VEP annotation sources (`vep_annotation`, defaults specified)
    - Cosmic Coding Variants VCF (v94)
    - Cosmic NonCoding Variants VCF (v94)
    - ClinVar VCF (20211002)
    - CADD (v1.6)
- MAF file created from known set of samples (default specified)
    - default provided is from **205** NovaSeq samples (as of 211007)
    - the process for generating this VCF and from what samples is documented [here](https://cuhbioinformatics.atlassian.net/wiki/spaces/URA/pages/2415591443/Creation+of+Myeloid+NovaSeq+samples+MAF+file)

## What does this app output?
- Excel workbook of annotated variants
- BSVI-friendly VCF
- Intermediate VCFs

## How does this app work?
- Filters mutect2 VCF with bedtools:
    - retain variants within ROI
- Filters mutect2 VCF with bcftools:
    - retain positions where at least one variant has AF > 0.03
    - retain positions where DP >99
    - split multiallelics using `--keep-sum AD` which changes the ref AD to be the sum of AD's
    - split multiallelics requires fixing AD and RPA number field in header from `.` to `R`
- Filters cgppindel VCF with bedtools:
    - only indels that intersect with the exons of interest bed file 
    - only insertions with length greater than 2. This will remove the 1 bp false positive insertions.
- Annotates mutect2 **and** cgppindel VCF with VEP:
    - Annotate against specified RefSeq transcripts with:
        - gene symbol
        - variant class
        - variant consequence
        - exon number
        - HGVS c. & p.
        - gnomAD AF
        - SIFT
        - PolyPhen
        - dbSNP
        - COSMIC
        - ClinVar
        - CADD
        - previous counts
- Filters VCFs with VEP:
    - Retain variants with gnomAD AF < 0.1
    - Remove synonymous variants
    - Filters by transcripts defined in `cpde.sh` for each panel
- Generates variant lists (one panel + pindel VCF per sheet) in an excel workbook
- Generates BSVI-friendly VCF

## Limitations
- Designed to be used as part of Uranus workflow for processing myeloid NGS panel
- When building the app, mark-section and mark-success must be executable for the app to build correctly - these two scripts are executable in this repo but it can be easy to miss if the executable bit is lost as it's a metadata change
