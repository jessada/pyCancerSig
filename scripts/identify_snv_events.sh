#!/bin/bash 
set -e
set -u
set -o pipefail

source $PYCMM/bash/cmm_functions.sh

module load samtools/1.9
module load annovar/2017.07.16
module load vt/0.5772

#define default values
GT_FORMAT_DEFAULT="GTR"

gt_format=$GT_FORMAT_DEFAULT

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-i {file}           input VCF file (required)
-g {file}           genotype format (default=$GT_FORMAT_DEFAULT)
-r {file}           path to genome reference (required)
-o {file}           output events file (required)
-h                  this help
EOF
)

# parse option
while getopts ":i:g:o:r:h" OPTION; do
  case "$OPTION" in
    i)
      input_vcf="$OPTARG"
      ;;
    g)
      gt_format="$OPTARG"
      ;;
    r)
      reference="$OPTARG"
      ;;
    o)
      out_event_file="$OPTARG"
      ;;
    h)
      echo >&2 "$usage"
      ;;
    *)
      die "unrecognized option (-$OPTION) from executing: $0 $@"
      ;;
  esac
done

[ ! -z $out_event_file ] || die "output event file is required (-o)"
[ ! -z $reference ] || die "reference is required (-r)"
[ -f "$input_vcf" ] || die "$input_vcf is not found"
[ -f "$reference" ] || die "$reference is not found"


time_stamp=$( date )

working_dir=`mktemp -d`

# ****************************************  display configuration  ****************************************
new_section_txt "S T A R T <$script_name>"
info_msg
info_msg "description"
info_msg "  This application will count the following changes in all possible 3'5'"
info_msg "    - C > A"
info_msg "    - C > G"
info_msg "    - C > T"
info_msg "    - T > A"
info_msg "    - T > C"
info_msg "    - T > G"
info_msg
info_msg "version and script configuration"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "input VCF file (-i)" "$input_vcf"
display_param "genotype format (-g)" "$gt_format"
display_param "genome reference (-r)" "$reference"
display_param "output event file (-o)" "$out_event_file"
display_param "working directory" "$working_dir"

# ****************************************  executing  ****************************************

new_section_txt "Decompose variants"

cmd="vcf-query -l"
cmd+=" $input_vcf"
cmd+=" | tr \"\n\" \"\t\""
samples_list=`eval $cmd`
n_samples=`echo $samples_list | wc -w`

tmp_decompose="$working_dir/tmp_decompose"
cmd="gunzip -c"
cmd+=" $input_vcf"
cmd+=" | vt decompose -s -"
cmd+=" | grep -Pv \"\t\*\t\""
cmd+=" | grep -v \"\\x3b\""
cmd+=" > $tmp_decompose"
eval_cmd "$cmd"

# suppress QUAL and INFO column
printf_phrase="%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t%s"
param_phrase="\$1, \$2, \$3, \$4, \$5, \$7, \$9"

for ((n_sample=1; n_sample<=$n_samples; n_sample++));
do
    printf_phrase+="\t%s"
    param_phrase+=", \$$((n_sample+9))"
done;

new_section_txt "Removing non-titv variants and empty QUAL and INFO columns"

tmp_removed_non_titv="$working_dir/tmp_removed_non_titv"
cmd="awk -F\$'\t'"
cmd+=" '{ if (length(\$4) == 1 && length(\$5) == 1) print \$0 }'"
cmd+=" $tmp_decompose"
cmd+=" | awk -F '\t' '{ printf \"$printf_phrase\n\", $param_phrase}'"
cmd+=" > $tmp_removed_non_titv"
eval_cmd "$cmd"

tmp_removed_non_titv_vcf="$working_dir/tmp_removed_non_titv.vcf"
cmd="tabix -h $input_vcf 1:1-1 > $tmp_removed_non_titv_vcf" 
eval_cmd "$cmd"
cmd="cat $tmp_removed_non_titv >> $tmp_removed_non_titv_vcf" 
eval_cmd "$cmd"

cmd="bgzip -f $tmp_removed_non_titv_vcf"
eval_cmd "$cmd"

cmd="tabix -p vcf $tmp_removed_non_titv_vcf.gz"
eval_cmd "$cmd"

new_section_txt "Annotate variants with refGene information"

tmp_ta="$working_dir/ta"
table_annovar_cmd="table_annovar.pl"
table_annovar_cmd+=" $tmp_removed_non_titv_vcf.gz"
table_annovar_cmd+=" $ANNOVAR_HOME/humandb"
table_annovar_cmd+=" -buildver hg19"
table_annovar_cmd+=" -out $tmp_ta"
table_annovar_cmd+=" -remove"
table_annovar_cmd+=" -protocol refGene"
table_annovar_cmd+=" -operation g"
table_annovar_cmd+=" -nastring ."
table_annovar_cmd+=" -vcfinput"
eval_cmd "$table_annovar_cmd"

multianno_vcf="$tmp_ta.hg19_multianno.vcf"
tmp_ta_vcf="$tmp_ta.vcf"
cmd="sed 's/Func.refGene/Func_refGene/g'"
cmd+=" $multianno_vcf"
cmd+=" | sed 's/Gene.refGene/Gene_refGene/g'"
cmd+=" > $tmp_ta_vcf"
eval_cmd "$cmd"

cmd="bgzip -f $tmp_ta_vcf"
eval_cmd "$cmd"

cmd="tabix -p vcf $tmp_ta_vcf.gz"
eval_cmd "$cmd"

new_section_txt "generate list of all titv events in all samples"

info_msg
info_msg ">>> extracting variants information <<<"

TMP_VARIANTS_LIST_CHROM_COL_IDX=1
TMP_VARIANTS_LIST_POS_COL_IDX=2
TMP_VARIANTS_LIST_REF_COL_IDX=3
TMP_VARIANTS_LIST_ALT_COL_IDX=4
TMP_VARIANTS_LIST_GENE_REFGENE_COL_IDX=5
TMP_VARIANTS_LIST_FIRST_GT_COL_IDX=6

vcf_query_format="'"
vcf_query_format+="%CHROM"
vcf_query_format+="\t%POS"
vcf_query_format+="\t%REF"
vcf_query_format+="\t%ALT"
vcf_query_format+="\t%INFO/Gene_refGene"
vcf_query_format+="[\t%$gt_format]"
vcf_query_format+="\n"
vcf_query_format+="'"

tmp_filtered_variants_list="$working_dir/tmp_filtered_variants_list"
cmd="vcf-query"
cmd+=" -f $vcf_query_format"
cmd+=" $tmp_ta_vcf.gz"
cmd+=" > $tmp_filtered_variants_list"
eval_cmd "$cmd"

tmp_gtz="$working_dir/tmp_gtz"
info_msg
info_msg
info_msg ">>> extracting zygosities <<<"
raw_cmd=" cut"
raw_cmd+=" -f"
for col_idx in $(seq $TMP_VARIANTS_LIST_FIRST_GT_COL_IDX $((TMP_VARIANTS_LIST_FIRST_GT_COL_IDX+n_samples-1)))
do
    raw_cmd+=",$col_idx"
done
raw_cmd+=" $tmp_filtered_variants_list"
raw_cmd+=" > $tmp_gtz"
cmd="$( echo $raw_cmd | sed 's/-f,/-f/g' )"
eval_cmd "$cmd"

info_msg
info_msg
info_msg ">>> extracting snp_triplets <<<"

tmp_get_snp_triplet_cmds="$working_dir/tmp_get_snp_triplet_cmds"
cmd="awk -F\$'\t'"
cmd+=" '{ printf \"samtools faidx $reference %s:%s-%s | tail -1\n\", \$$TMP_VARIANTS_LIST_CHROM_COL_IDX, \$$TMP_VARIANTS_LIST_POS_COL_IDX-1, \$$TMP_VARIANTS_LIST_POS_COL_IDX+1 }'"
cmd+=" $tmp_filtered_variants_list"
cmd+=" > $tmp_get_snp_triplet_cmds"
eval_cmd "$cmd"

tmp_snp_triplets="$working_dir/tmp_snp_triplets"
:>$tmp_snp_triplets
while read cmd; do
    eval "$cmd >> $tmp_snp_triplets"
done <  $tmp_get_snp_triplet_cmds

info_msg
info_msg
info_msg ">>> extracting coordinates <<<"

tmp_coors="$working_dir/tmp_coors"
cmd="awk -F\$'\t'"
cmd+=" '{ printf \"%s\t%s\t%s\t%s\n\", \$$TMP_VARIANTS_LIST_CHROM_COL_IDX, \$$TMP_VARIANTS_LIST_POS_COL_IDX, \$$TMP_VARIANTS_LIST_REF_COL_IDX, \$$TMP_VARIANTS_LIST_ALT_COL_IDX }'"
cmd+=" $tmp_filtered_variants_list"
cmd+=" > $tmp_coors"
eval_cmd "$cmd"

info_msg
info_msg
info_msg ">>> extracting genes <<<"

tmp_genes="$working_dir/tmp_genes"
cmd="awk -F\$'\t'"
cmd+=" '{ printf \"%s\n\", \$$TMP_VARIANTS_LIST_GENE_REFGENE_COL_IDX }'"
cmd+=" $tmp_filtered_variants_list"
cmd+=" > $tmp_genes"
eval_cmd "$cmd"

info_msg
info_msg
info_msg ">>> merging genes, coordinates, snp_triplets, and zygosities together <<<"

tmp_genes_coors_snp_triplets_gtz="$working_dir/tmp_genes_coors_snp_triplets_gtz"
cmd="paste"
cmd+=" $tmp_genes"
cmd+=" $tmp_coors"
cmd+=" $tmp_snp_triplets"
cmd+=" $tmp_gtz"
cmd+=" > $tmp_genes_coors_snp_triplets_gtz"
eval_cmd "$cmd"

info_msg
info_msg
info_msg ">>> extracting transcription strands <<<"

tmp_strand_genes="$working_dir/tmp_strand_genes"
cmd="join"
cmd+=" -1 1"
cmd+=" -2 2"
cmd+=" -t $'\t'"
cmd+=" -o 2.1,1.1"
cmd+=" <( sort $tmp_genes )"
cmd+=" <( cut -f4,13 $ANNOVAR_HOME/humandb/hg19_refGene.txt | sort -k2,2 | uniq )"
cmd+=" | sort -k2,2"
cmd+=" | uniq -f1"
cmd+=" > $tmp_strand_genes"
eval_cmd "$cmd"

info_msg
info_msg
info_msg ">>> merging strands, genes, coordinates, snp_triplets, and zygosities together <<<"

echo -e "#Strand\tGene\tChr\tPos\tRef\tAlt\tsnp_triplet\t$samples_list" > $out_event_file

out_join_phrase="1.1,0,2.2,2.3,2.4,2.5,2.6"

for ((n_sample=1; n_sample<=$n_samples; n_sample++));
do
    out_join_phrase+=",2.$((n_sample+6))"
done

cmd="join"
cmd+=" -1 2"
cmd+=" -2 1"
cmd+=" -t $'\t'"
cmd+=" -o $out_join_phrase"
cmd+=" <( sort -k2,2 $tmp_strand_genes )"
cmd+=" <( sort -k1,1 $tmp_genes_coors_snp_triplets_gtz )"
cmd+=" >> $out_event_file"
eval_cmd "$cmd"

new_section_txt "F I N I S H <$script_name>"

