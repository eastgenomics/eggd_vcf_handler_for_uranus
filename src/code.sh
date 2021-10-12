#!/bin/bash
#
# Performs annotation and filtering of given VCF to produce multiple annotated VCFs
# Generates an excel workbook to aid variant interpretation
# Also generates a workaround for BSVI mis-handling multiallelics

function annotate_vep_vcf {
	# Function to run VEP for annotation on given VCF file

	# Inputs:
	# $1 -> input vcf (splitfile currently for all)
	# $2 -> name for output vcf
	# $3 -> list of comma separated transcripts to annotate and filter on

	input_vcf="$1"
	output_vcf="$2"
	transcript_list="$3"

	# # transcript list should be comma separated, format for passing to VEP
	transcript_list=$(echo "$transcript_list" | sed 's/[[:space:]]//g')  # remove any whitespace
	transcript_list=$(echo "$transcript_list" | sed 's/,/ or stable_id match /g')
	transcript_list="stable_id match $transcript_list"
	
	# fields to filter on
	# hard coded in function for now, can be made an input but all are the same
	filter_fields="SYMBOL,VARIANT_CLASS,Consequence,EXON,HGVSc,HGVSp,gnomAD_AF,CADD_PHRED,Existing_variation,ClinVar,ClinVar_CLNDN,ClinVar_CLNSIG,Prev_AC,Prev_NS"

	# find clinvar vcf, remove leading ./
	clinvar_vcf=$(find ./ -name "clinvar_*.vcf.gz" | sed s'/.\///')

	time docker run -v /home/dnanexus:/opt/vep/.vep \
	ensemblorg/ensembl-vep:release_103.1 \
	./vep -i /opt/vep/.vep/"${input_vcf}" -o /opt/vep/.vep/"${output_vcf}" \
	--vcf --cache --refseq --exclude_predicted --symbol --hgvs --af_gnomad \
	--check_existing --variant_class --numbers \
	--offline \
	--custom /opt/vep/.vep/"${clinvar_vcf}",ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
	--custom /opt/vep/.vep/"${maf_file_name}",Prev,vcf,exact,0,AC,NS \
	--plugin CADD,/opt/vep/.vep/whole_genome_SNVs.tsv.gz,/opt/vep/.vep/gnomad.genomes.r3.0.indel.tsv.gz \
	--fields "$filter_fields" \
	--no_stats \
	--transcript_filter "$transcript_list"
}

function filter_vep_vcf {
	# Function to filter annotated VCF with VEP to retain variants with AF < 0.10 in gnomAD

	# Inputs:
	# 	$1 -> input vcf (should be output vcf of annotation)
	# 	$2 -> name for output_vcf

	input_vcf="$1"
	output_vcf="$2"

	time docker run -v /home/dnanexus:/opt/vep/.vep \
	ensemblorg/ensembl-vep:release_103.1 \
	./filter_vep -i /opt/vep/.vep/"$input_vcf" \
	-o /opt/vep/.vep/"$output_vcf" --only_matched --filter \
	"(gnomAD_AF < 0.10 or not gnomAD_AF) and SYMBOL"
}

main() {
	set -e -x -v -o pipefail

	mark-section "downloading inputs"
	time dx-download-all-inputs --parallel

	# move maf file and index to home for vep to find
	mv "${maf_file_path}" /home/dnanexus/
	mv "${maf_file_tbi_path}" /home/dnanexus/

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

	time bedtools intersect -header -a "${vcf_path}" -b "${bed_path}" \
	| bcftools view -i "FORMAT/AF[*]>0.03" - \
	| bcftools view -i "DP>99" - \
	| sed 's/AD,Number=./AD,Number=R/g' \
	| sed 's/RPA,Number=./RPA,Number=R/g' \
	| bcftools norm -f "${mutect2_fasta_path}" -m -any --keep-sum AD - \
	-o ~/"${splitfile}"

	mark-section "annotating and further filtering"
	
	# vep needs permissions to write to /home/dnanexus
	chmod a+rwx /home/dnanexus
	
	# extract vep tarball (input) to /home/dnanexus
	# --transform used to ensure everything is extracted and not in a sub folder
	time tar xf "${vep_tarball_path}" -C /home/dnanexus --transform='s/.*\///'
	
	# extract annotation tarball to /home/dnanexus
	time tar xf ~/homo_sapiens_refseq_vep_103_GRCh38.tar.gz

	# place fasta and indexes for VEP in the annotation folder
	mv ~/Homo_sapiens.GRCh38.dna.toplevel.fa.gz ~/homo_sapiens_refseq/103_GRCh38/
	mv ~/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai ~/homo_sapiens_refseq/103_GRCh38/
	mv ~/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.gzi ~/homo_sapiens_refseq/103_GRCh38/

	# place plugins into plugins folder
	mkdir ~/Plugins
	mv ~/CADD.pm ~/Plugins/
	mv ~/plugin_config.txt ~/Plugins/

	# load vep docker (asset)
	docker load -i ~/vep_v103.1_docker.tar.gz

	# will run VEP to annotate against specified transcripts for all,
	# lymphoid and myeloid gene lists
	# for each, a transcript list is used to annotate on, as well as a *file (input vcf) and a
	# *vepfile (output VCF annotated and filtered by VEP)

	# full gene transcript list
	all_genes_transcripts="NM_002074,NM_000760,NM_005373,NM_002227,NM_002524,NM_022552,NM_012433,\
	NM_005896,NM_002468,NM_032638,NM_000222,NM_001127208,NM_033632,NM_002520,NM_016222,NM_006060,\
	NM_181500,NM_004333,NM_004456,NM_170606,NM_006265,NM_004972,NM_016734,NM_017617,NM_000314,\
	NM_005343,NM_024426,NM_001165,NM_000051,NM_001197104,NM_005188,NM_001987,NM_018638,NM_033360,\
	NM_001136023,NM_005475,NM_002834,NM_004119,NM_002168,NM_004380,NM_000546,NM_001042492,\
	NM_012448,NM_139276,NM_003620,NM_001195427,NM_015559,NM_004343,NM_004364,NM_015338,NM_080425,\
	NM_001754,NM_006758,NM_007194,NM_001429,NM_005089,NM_001123385,NM_002049,NM_001042750,\
	NM_001184772,NM_001015877"

	# annotate full VCF with VEP and all gene transcripts
	splitvepfile="${vcf_prefix}_split_filevep.vcf"
	annotate_vep_vcf "$splitfile" "$splitvepfile" "$all_genes_transcripts"


	# run vep for all gene transcripts
	allgenesfile="${vcf_prefix}_allgenes.vcf"
	allgenesvepfile="${vcf_prefix}_allgenesvep.vcf"

	annotate_vep_vcf "$splitfile" "$allgenesfile" "$all_genes_transcripts"
	filter_vep_vcf "$allgenesfile" "$allgenesvepfile"


	# run VEP for lymphoid genes list
	lymphoid_transcripts="NM_000051,NM_001165,NM_004333,NM_004380,NM_001429,NM_004456,\
	NM_033632,NM_005343,NM_033360,NM_002468,NM_017617,NM_002524,NM_016734,NM_012433,\
	NM_139276,NM_012448,NM_000546"

	lymphoidfile="${vcf_prefix}_lymphoid.vcf"
	lymphoidvepfile="${vcf_prefix}_pan-lymphoidvep.vcf"

	annotate_vep_vcf "$splitfile" "$lymphoidfile" "$lymphoid_transcripts"
	filter_vep_vcf "${lymphoidfile}" "$lymphoidvepfile"

	# run VEP for myeloid genes list
	myeloid_transcripts="NM_015338,NM_001123385,NM_001184772,NM_004333,NM_004343,NM_005188,\
	NM_004364,NM_007194,NM_000760,NM_181500,NM_016222,NM_022552,NM_018638,NM_001987,NM_004456,\
	NM_033632,NM_004119,NM_002049,NM_032638,NM_080425,NM_002074,NM_005343,NM_005896,NM_002168,\
	NM_006060,NM_002227,NM_004972,NM_000222,NM_001197104,NM_033360,NM_005373,NM_001042492,\
	NM_001136023,NM_017617,NM_002520,NM_002524,NM_016734,NM_001015877,NM_003620,NM_000314,\
	NM_002834,NM_006265,NM_001754,NM_015559,NM_012433,NM_005475,NM_001195427,NM_001042750,\
	NM_139276,NM_012448,NM_001127208,NM_000546,NM_006758,NM_024426,NM_005089"

	myeloidfile="${vcf_prefix}_myeloid.vcf"
	myeloidvepfile="${vcf_prefix}_myeloidvep.vcf"

	annotate_vep_vcf "$splitfile" "$myeloidfile" "$myeloid_transcripts"
	filter_vep_vcf "${myeloidfile}" "$myeloidvepfile"

	# run VEP for CLL_Extended genes list
	cll_transcripts="NM_001165,NM_004333,NM_033632,NM_005343,NM_033360,NM_002468,NM_017617,\
	NM_002524,NM_012433,NM_000546"

	cllfile="${vcf_prefix}_CLL-extended.vcf"
	cllvepfile="${vcf_prefix}_CLL-extendedvep.vcf"

	annotate_vep_vcf "$splitfile" "$cllfile" "$cll_transcripts"
	filter_vep_vcf "${cllfile}" "$cllvepfile"

	# run VEP for TP53
	tp53file="${vcf_prefix}_TP53.vcf"
	tp53vepfile="${vcf_prefix}_TP53vep.vcf"

	annotate_vep_vcf "$splitfile" "$tp53file" "NM_000546"
	filter_vep_vcf "${tp53file}" "$tp53vepfile"

	# run VEP for LGL
	lgl_transcripts="NM_139276,NM_012448"
	lglfile="${vcf_prefix}_LGL.vcf"
	lglvepfile="${vcf_prefix}_LGLvep.vcf"

	annotate_vep_vcf "$splitfile" "$lglfile" "$lgl_transcripts"
	filter_vep_vcf "${lglfile}" "$lglvepfile"

	# run vep for HCL
	hclfile="${vcf_prefix}_HCL.vcf"
	hclvepfile="${vcf_prefix}_HCLvep.vcf"

	annotate_vep_vcf "$splitfile" "$hclfile" "NM_004333"
	filter_vep_vcf "${hclfile}" "$hclvepfile"

	# run vep for LPL
	lplfile="${vcf_prefix}_LPL.vcf"
	lplvepfile=${vcf_prefix}_LPLvep.vcf

	annotate_vep_vcf "$splitfile" "$lplfile" "NM_002468"
	filter_vep_vcf "$lplfile" "$lplvepfile"


	mark-section "BSVI workaround (overwriting GT) and creating variant list"

	# install required python packages (asset)
	pip3 install /pytz-*.whl /numpy-*.whl /pandas-*.whl /XlsxWriter-*.whl

	# call Python script (asset) to spit multiallelics, generate BSVI VCF and excel report
	# note that order of VCFs passed to -v determines order of sheets in excel
	time python3 vcf_handler.py -a "${allgenesvepfile}" \
	-v "${myeloidvepfile}" "${cllvepfile}" "${tp53vepfile}" "${lglvepfile}" \
	"${hclvepfile}" "${lplvepfile}" "${lymphoidvepfile}"

	mark-section "uploading output"

	# make required output directories and move files
	mkdir -p ~/out/allgenes_filtered_vcf ~/out/lymphoid_filtered_vcf/ ~/out/myeloid_filtered_vcf/\
		~/out/cll_filtered_vcf ~/out/tp53_filtered_vcf ~/out/lgl_filtered_vcf ~/out/hcl_filtered_vcf\
		~/out/lpl_filtered_vcf ~/out/bsvi_vcf ~/out/text_report ~/out/excel_report

	bsvivcf="${vcf_prefix}_allgenes_bsvi.vcf"
	variantlist="${vcf_prefix}_allgenes.tsv"
	excellist="${vcf_prefix}_panels.xlsx"

	mv ~/"${allgenesvepfile}" ~/out/allgenes_filtered_vcf/
	mv ~/"${lymphoidvepfile}" ~/out/lymphoid_filtered_vcf/
	mv ~/"${myeloidvepfile}" ~/out/myeloid_filtered_vcf/
	mv ~/"${cllvepfile}" ~/out/cll_filtered_vcf/
	mv ~/"${tp53vepfile}" ~/out/tp53_filtered_vcf/
	mv ~/"${lglvepfile}" ~/out/lgl_filtered_vcf/
	mv ~/"${hclvepfile}" ~/out/hcl_filtered_vcf/
	mv ~/"${lplvepfile}" ~/out/lpl_filtered_vcf/
	mv ~/"${bsvivcf}" ~/out/bsvi_vcf/
	mv ~/"${variantlist}" ~/out/text_report/
	mv ~/"${excellist}" ~/out/excel_report/
	
	dx-upload-all-outputs --parallel
	mark-success
}
