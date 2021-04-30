#!/bin/bash
#
# Pulls together a bunch of proceedures to provide a workaround for BSVI
# Also generates a variantlist text tsv for clinical scientists

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

mark-section "downloading inputs"
dx-download-all-inputs --parallel

mark-section "filtering and splitting multiallelics"
# retain variants that are: # within ROIs (bed file),
#   have at least one allele >0.03 AF, and have DP >99
# fix AD and RPA number in header
# split multiallelics using --keep-sum AD which changes the ref AD to be a sum
#   of all other AD's rather than being ref AD alone
# note that --keep-sum AD is a one way conversion in bcftools 1.12.0 and can't
#   be undone with bcftools norm -m +any
# bedtools and bcftools are app assets
splitfile="${vcf_prefix}_split.vcf"
bedtools intersect -header -a "${vcf_path}" -b "${bed_path}" \
  | bcftools view -i "FORMAT/AF[*]>0.03" - \
  | bcftools view -i "DP>99" - \
  | sed 's/AD,Number=./AD,Number=R/g' \
  | sed 's/RPA,Number=./RPA,Number=R/g' \
  | bcftools norm -f "${mutect2_fasta_path}" -m -any --keep-sum AD - \
  -o ~/"${splitfile}"


mark-section "annotating and further filtering"
# vep needs permissions to write to /home/dnanexus
chmod a+rwx /home/dnanexus
# extract annotation tarball (asset) to /home/dnanexus
tar xf /homo_sapiens_refseq_vep_103_GRCh38.tar.gz -C /home/dnanexus
# place fasta and indexes for VEP (assets) in the annotation folder
mv /Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
  ~/homo_sapiens_refseq/103_GRCh38/
mv /Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai \
  ~/homo_sapiens_refseq/103_GRCh38/
mv /Homo_sapiens.GRCh38.dna.toplevel.fa.gz.gzi \
  ~/homo_sapiens_refseq/103_GRCh38/
# load vep docker (asset)
docker load -i /vep_v103.1_docker.tar.gz
# will run vep to annotate against specified transcripts for all, lymphoid
# and myeloid gene lists
# run vep for all genes list
allgenesfile="${vcf_prefix}_allgenes.vcf"
docker run -v /home/dnanexus:/opt/vep/.vep \
  ensemblorg/ensembl-vep:release_103.1 \
  ./vep -i /opt/vep/.vep/"${splitfile}" -o /opt/vep/.vep/"${allgenesfile}" \
  --vcf --cache --refseq --exclude_predicted --symbol --hgvs --af_gnomad \
  --check_existing --variant_class --numbers -sift s --polyphen s --fields \
  "SYMBOL,VARIANT_CLASS,Consequence,EXON,HGVSc,HGVSp,gnomAD_AF,SIFT,PolyPhen,Existing_variation" \
  --no_stats --transcript_filter \
  "stable_id match NM_002074 \
  or stable_id match NM_000760 \
  or stable_id match NM_005373 \
  or stable_id match NM_002227 \
  or stable_id match NM_002524 \
  or stable_id match NM_022552 \
  or stable_id match NM_012433 \
  or stable_id match NM_005896 \
  or stable_id match NM_002468 \
  or stable_id match NM_032638 \
  or stable_id match NM_000222 \
  or stable_id match NM_001127208 \
  or stable_id match NM_033632 \
  or stable_id match NM_002520 \
  or stable_id match NM_016222 \
  or stable_id match NM_006060 \
  or stable_id match NM_181500 \
  or stable_id match NM_004333 \
  or stable_id match NM_004456 \
  or stable_id match NM_170606 \
  or stable_id match NM_006265 \
  or stable_id match NM_004972 \
  or stable_id match NM_016734 \
  or stable_id match NM_017617 \
  or stable_id match NM_000314 \
  or stable_id match NM_005343 \
  or stable_id match NM_024426 \
  or stable_id match NM_001165 \
  or stable_id match NM_000051 \
  or stable_id match NM_001197104 \
  or stable_id match NM_005188 \
  or stable_id match NM_001987 \
  or stable_id match NM_018638 \
  or stable_id match NM_033360 \
  or stable_id match NM_001136023 \
  or stable_id match NM_005475 \
  or stable_id match NM_002834 \
  or stable_id match NM_004119 \
  or stable_id match NM_002168 \
  or stable_id match NM_004380 \
  or stable_id match NM_000546 \
  or stable_id match NM_001042492 \
  or stable_id match NM_012448 \
  or stable_id match NM_139276 \
  or stable_id match NM_003620 \
  or stable_id match NM_001195427 \
  or stable_id match NM_015559 \
  or stable_id match NM_004343 \
  or stable_id match NM_004364 \
  or stable_id match NM_015338 \
  or stable_id match NM_080425 \
  or stable_id match NM_001754 \
  or stable_id match NM_006758 \
  or stable_id match NM_007194 \
  or stable_id match NM_001429 \
  or stable_id match NM_005089 \
  or stable_id match NM_001123385 \
  or stable_id match NM_002049 \
  or stable_id match NM_001042750 \
  or stable_id match NM_001184772 \
  or stable_id match NM_001015877"
# run vep filtering to remove high AF and outside genelist
allgenesvepfile="${vcf_prefix}_allgenesvep.vcf"
docker run -v /home/dnanexus:/opt/vep/.vep \
  ensemblorg/ensembl-vep:release_103.1 \
  ./filter_vep -i /opt/vep/.vep/"${allgenesfile}" \
  -o /opt/vep/.vep/"${allgenesvepfile}" --only_matched --filter \
  "(gnomAD_AF < 0.10 or not gnomAD_AF)"
# run vep for lymphoid genes list
lymphoidfile="${vcf_prefix}_lymphoid.vcf"
docker run -v /home/dnanexus:/opt/vep/.vep \
  ensemblorg/ensembl-vep:release_103.1 \
  ./vep -i /opt/vep/.vep/"${splitfile}" -o /opt/vep/.vep/"${lymphoidfile}" \
  --vcf --cache --refseq --exclude_predicted --symbol --hgvs --af_gnomad \
  --check_existing --variant_class --numbers -sift s --polyphen s --fields \
  "SYMBOL,VARIANT_CLASS,Consequence,EXON,HGVSc,HGVSp,gnomAD_AF,SIFT,PolyPhen,Existing_variation" \
  --no_stats --transcript_filter \
  "stable_id match NM_000051 \
  or stable_id match NM_001165 \
  or stable_id match NM_004333 \
  or stable_id match NM_004380 \
  or stable_id match NM_001429 \
  or stable_id match NM_004456 \
  or stable_id match NM_033632 \
  or stable_id match NM_005343 \
  or stable_id match NM_033360 \
  or stable_id match NM_002468 \
  or stable_id match NM_017617 \
  or stable_id match NM_002524 \
  or stable_id match NM_016734 \
  or stable_id match NM_012433 \
  or stable_id match NM_139276 \
  or stable_id match NM_012448 \
  or stable_id match NM_000546"
# run vep filtering to remove high AF and outside genelist
lymphoidvepfile="${vcf_prefix}_lymphoidvep.vcf"
docker run -v /home/dnanexus:/opt/vep/.vep \
  ensemblorg/ensembl-vep:release_103.1 \
  ./filter_vep -i /opt/vep/.vep/"${lymphoidfile}" \
  -o /opt/vep/.vep/"${lymphoidvepfile}" --only_matched --filter \
  "(gnomAD_AF < 0.10 or not gnomAD_AF)" --filter "SYMBOL"
# run vep for myeloid genes list
myeloidfile="${vcf_prefix}_myeloid.vcf"
docker run -v /home/dnanexus:/opt/vep/.vep \
  ensemblorg/ensembl-vep:release_103.1 \
  ./vep -i /opt/vep/.vep/"${splitfile}" -o /opt/vep/.vep/"${myeloidfile}" \
  --vcf --cache --refseq --exclude_predicted --symbol --hgvs --af_gnomad \
  --check_existing --variant_class --numbers -sift s --polyphen s --fields \
  "SYMBOL,VARIANT_CLASS,Consequence,EXON,HGVSc,HGVSp,gnomAD_AF,SIFT,PolyPhen,Existing_variation" \
  --no_stats --transcript_filter \
  "stable_id match NM_015338 \
  or stable_id match NM_001123385 \
  or stable_id match NM_001184772 \
  or stable_id match NM_004333 \
  or stable_id match NM_004343 \
  or stable_id match NM_005188 \
  or stable_id match NM_004364 \
  or stable_id match NM_007194 \
  or stable_id match NM_000760 \
  or stable_id match NM_181500 \
  or stable_id match NM_016222 \
  or stable_id match NM_022552 \
  or stable_id match NM_018638 \
  or stable_id match NM_001987 \
  or stable_id match NM_004456 \
  or stable_id match NM_033632 \
  or stable_id match NM_004119 \
  or stable_id match NM_002049 \
  or stable_id match NM_032638 \
  or stable_id match NM_080425 \
  or stable_id match NM_002074 \
  or stable_id match NM_005343 \
  or stable_id match NM_005896 \
  or stable_id match NM_002168 \
  or stable_id match NM_006060 \
  or stable_id match NM_002227 \
  or stable_id match NM_004972 \
  or stable_id match NM_000222 \
  or stable_id match NM_001197104 \
  or stable_id match NM_033360 \
  or stable_id match NM_005373 \
  or stable_id match NM_001042492 \
  or stable_id match NM_001136023 \
  or stable_id match NM_017617 \
  or stable_id match NM_002520 \
  or stable_id match NM_002524 \
  or stable_id match NM_016734 \
  or stable_id match NM_001015877 \
  or stable_id match NM_003620 \
  or stable_id match NM_000314 \
  or stable_id match NM_002834 \
  or stable_id match NM_006265 \
  or stable_id match NM_001754 \
  or stable_id match NM_015559 \
  or stable_id match NM_012433 \
  or stable_id match NM_005475 \
  or stable_id match NM_001195427 \
  or stable_id match NM_001042750 \
  or stable_id match NM_139276 \
  or stable_id match NM_012448 \
  or stable_id match NM_001127208 \
  or stable_id match NM_000546 \
  or stable_id match NM_006758 \
  or stable_id match NM_024426 \
  or stable_id match NM_005089"
# run vep filtering to remove high AF and outside genelist
myeloidvepfile="${vcf_prefix}_myeloidvep.vcf"
docker run -v /home/dnanexus:/opt/vep/.vep \
  ensemblorg/ensembl-vep:release_103.1 \
  ./filter_vep -i /opt/vep/.vep/"${myeloidfile}" \
  -o /opt/vep/.vep/"${myeloidvepfile}" --only_matched --filter \
  "(gnomAD_AF < 0.10 or not gnomAD_AF)" --filter "SYMBOL"


mark-section "BSVI workaround (overwriting GT) and creating variant list"
# install required python packages (asset)
pip3 install /pytz-*.whl /numpy-*.whl /pandas-*.whl /XlsxWriter-*.whl
# call Python script (asset) to spit multiallelics
# and generate BSVI VCF & text report
python3 vcf_handler.py -a "${allgenesvepfile}" -v "${lymphoidfile}" "${myeloidfile}"


mark-section "uploading output"
mkdir -p ~/out/allgenes_filtered_vcf
mv ~/"${allgenesvepfile}" ~/out/allgenes_filtered_vcf/
mkdir -p ~/out/lymphoid_filtered_vcf
mv ~/"${lymphoidvepfile}" ~/out/lymphoid_filtered_vcf/
mkdir -p ~/out/myeloid_filtered_vcf
mv ~/"${myeloidvepfile}" ~/out/myeloid_filtered_vcf/
bsvivcf="${vcf_prefix}_allgenes_bsvi.vcf"
mkdir -p ~/out/bsvi_vcf
mv ~/"${bsvivcf}" ~/out/bsvi_vcf/
variantlist="${vcf_prefix}_allgenes.tsv"
mkdir -p ~/out/text_report
mv ~/"${variantlist}" ~/out/text_report/
excellist="${vcf_prefix}_panels.xlsx"
mkdir -p ~/out/excel_report
mv ~/"${excellist}" ~/out/excel_report/

dx-upload-all-outputs --parallel

mark-success
