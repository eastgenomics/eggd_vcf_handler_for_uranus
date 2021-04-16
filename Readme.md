# eggd vcf handler for uranus


## What does this app do?

This app bundles up a set of assets to enable taking the VCFs from sentieon mutect2 and:
- Annotates and filters the sentieon mutect2 VCF
- Produces a tsv variant list that provides clinical scientists with a complete list of variants with some annotation and preformatted text for Epic
- Produces a VCF that can be used as input for BSVI, which is incompatible with mutect2's (correct) representation of multiallelic variants in VCF

## What are typical use cases for this app?
This app uses the follow tools which are app assets:
* bcftools (v1.12)
* bedtools (v2.30.0)
* VEP (v103.1) (docker image)

This app has access to the Internet

## What data are required for this app to run?
- VCF outputted by sentieon mutect2 as part of Uranus workflow
- BED file that details ROIs for myeloid NGS panel
- Genome FASTA and index that was used by to generate the VCF

## What does this app output?
- TSV variant list
- BSVI-friendly VCF

## How does this app work?
- Filters VCF with bedtools
- Filters VCF with bcftools
- Annotates VCF with VEP
- Filters VCF with VEP
- Generates TSV variant list
- Generate BSVI-friendly VCF

## Limitations
- Designed to be used as part of Uranus workflow for processing myeloid NGS panel
