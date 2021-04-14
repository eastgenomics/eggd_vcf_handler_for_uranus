#!/bin/bash
#

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

mark-section "downloading inputs"
dx-download-all-inputs --parallel

mark-section "filtering and splitting multiallelics"
# retain variants that are within ROIs (bed file), have at least one allele >0.03 AF, have DP >99
# fix AD and RPA number in header
# split multiallelics using --keep-sum AD which changes the ref AD to be a sum of all other AD's rather than being ref AD alone
# note that --keep-sum AD is a one way conversion in bcftools 1.12.0 and can't be undone with bcftools norm -m +any
bedtools intersect -header -a $vcf_path -b $bed_path | bcftools view -i "FORMAT/AF[*]>0.03" - | bcftools view -i "DP>99" - | sed 's/AD,Number=./AD,Number=R/g' | sed 's/RPA,Number=./RPA,Number=R/g' | bcftools norm -f $mutect2_fasta_path -m -any --keep-sum AD - -o ~/$vcf_prefix.split.vcf

mark-section "annotating and further filtering"
# vep needs permissions to write to /home/dnanexus
chmod a+rwx /home/dnanexus
# extract annotation tarball to /home/dnanexus
tar xf $vep_annotation_path -C /home/dnanexus
# place fasta and indexes for VEP in the annotation folder
mv $vep_fasta_path /home/dnanexus/homo_sapiens_refseq/103_GRCh38/
mv $vep_fai_path /home/dnanexus/homo_sapiens_refseq/103_GRCh38/
mv $vep_gzi_path /home/dnanexus/homo_sapiens_refseq/103_GRCh38/
# load vep docker
docker load -i $vep_docker_path
# run vep to annotate against specified transcripts
docker run -v /home/dnanexus:/opt/vep/.vep ensemblorg/ensembl-vep:release_103.1 ./vep -i /opt/vep/.vep/$vcf_prefix.split.vcf -o /opt/vep/.vep/$vcf_prefix.annotate.vcf --vcf --cache --refseq --exclude_predicted --symbol --hgvs --af_gnomad --check_existing --variant_class --numbers -sift s --polyphen s --fields "SYMBOL,VARIANT_CLASS,Consequence,EXON,HGVSc,HGVSp,gnomAD_AF,SIFT,PolyPhen,Existing_variation" --no_stats --transcript_filter "stable_id match NM_002074 or stable_id match NM_000760 or stable_id match NM_005373 or stable_id match NM_002227 or stable_id match NM_002524 or stable_id match NM_022552 or stable_id match NM_012433 or stable_id match NM_005896 or stable_id match NM_002468 or stable_id match NM_032638 or stable_id match NM_000222 or stable_id match NM_001127208 or stable_id match NM_033632 or stable_id match NM_002520 or stable_id match NM_016222 or stable_id match NM_006060 or stable_id match NM_181500 or stable_id match NM_004333 or stable_id match NM_004456 or stable_id match NM_170606 or stable_id match NM_006265 or stable_id match NM_004972 or stable_id match NM_016734 or stable_id match NM_017617 or stable_id match NM_000314 or stable_id match NM_005343 or stable_id match NM_024426 or stable_id match NM_001165 or stable_id match NM_000051 or stable_id match NM_001197104 or stable_id match NM_005188 or stable_id match NM_001987 or stable_id match NM_018638 or stable_id match NM_033360 or stable_id match NM_001136023 or stable_id match NM_005475 or stable_id match NM_002834 or stable_id match NM_004119 or stable_id match NM_002168 or stable_id match NM_004380 or stable_id match NM_000546 or stable_id match NM_001042492 or stable_id match NM_012448 or stable_id match NM_139276 or stable_id match NM_003620 or stable_id match NM_001195427 or stable_id match NM_015559 or stable_id match NM_004343 or stable_id match NM_004364 or stable_id match NM_015338 or stable_id match NM_080425 or stable_id match NM_001754 or stable_id match NM_006758 or stable_id match NM_007194 or stable_id match NM_001429 or stable_id match NM_005089 or stable_id match NM_001123385 or stable_id match NM_002049 or stable_id match NM_001042750 or stable_id match NM_001184772 or stable_id match NM_001015877"
# run vep filtering to remove high AF and synonymous
docker run -v /home/dnanexus:/opt/vep/.vep ensemblorg/ensembl-vep:release_103.1 ./filter_vep -i /opt/vep/.vep/$vcf_prefix.annotate.vcf -o /opt/vep/.vep/$vcf_prefix.vepfilter.vcf --only_matched --filter "(gnomAD_AF < 0.10 or not gnomAD_AF) and (Consequence != synonymous_variant)"

mark-section "BSVI workaround (overwriting GT) and creating text report"
# install required python packages from local tar
tar xf $python_packages_path -C /home/dnanexus
pip3 install pytz-*.whl numpy-*.whl pandas-*.whl
# call Python script to spit multiallelics and generate BSVI VCF & text report
python3 $python_script_path $vcf_prefix.vepfilter.vcf

mark-section "uploading output"
mkdir -p ~/out/filtered_annotated_vcf
mv ~/$vcf_prefix.vepfilter.vcf ~/out/filtered_annotated_vcf/
mkdir -p ~/out/bsvi_vcf
mv ~/$vcf_prefix.bsvi.vcf ~/out/bsvi_vcf
mkdir -p ~/out/text_report
mv ~/$vcf_prefix.variantlist.tsv ~/out/text_report

dx-upload-all-outputs --parallel

mark-success
