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
frequency_ratio="1"

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
-f {percent}        filtering frequency ratio (detault=1)
-e {file}           output events file (required)
-o {file}           output titv file (required)
-h                  this help
EOF
)

# parse option
while getopts ":i:g:e:r:f:o:h" OPTION; do
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
    f)
      frequency_ratio="$OPTARG"
      ;;
    e)
      out_event_file="$OPTARG"
      ;;
    o)
      out_titv_file="$OPTARG"
      ;;
    h)
      echo >&2 "$usage"
      ;;
    *)
      die "unrecognized option (-$OPTION) from executing: $0 $@"
      ;;
  esac
done

[ ! -z $out_event_file ] || die "output events file is required (-a)"
[ ! -z $out_titv_file ] || die "output titv count is required (-o)"
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
display_param "filtering frequency ratiio (-f)" "$frequency_ratio"
display_param "output event file (-e)" "$out_event_file"
display_param "output titv file (-o)" "$out_titv_file"
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
table_annovar_cmd+=" -protocol refGene,gnomad_exome"
table_annovar_cmd+=" -operation g,f"
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
TMP_VARIANTS_LIST_GNOMAD_EXOME_ALL_COL_IDX=6
TMP_VARIANTS_LIST_GNOMAD_EXOME_AFR_COL_IDX=7
TMP_VARIANTS_LIST_GNOMAD_EXOME_AMR_COL_IDX=8
TMP_VARIANTS_LIST_GNOMAD_EXOME_ASJ_COL_IDX=9
TMP_VARIANTS_LIST_GNOMAD_EXOME_EAS_COL_IDX=10
TMP_VARIANTS_LIST_GNOMAD_EXOME_FIN_COL_IDX=11
TMP_VARIANTS_LIST_GNOMAD_EXOME_NFE_COL_IDX=12
TMP_VARIANTS_LIST_GNOMAD_EXOME_OTH_COL_IDX=13
TMP_VARIANTS_LIST_FIRST_GT_COL_IDX=14

vcf_query_format="'"
vcf_query_format+="%CHROM"
vcf_query_format+="\t%POS"
vcf_query_format+="\t%REF"
vcf_query_format+="\t%ALT"
vcf_query_format+="\t%INFO/Gene_refGene"
vcf_query_format+="\t%INFO/gnomAD_exome_ALL"
vcf_query_format+="\t%INFO/gnomAD_exome_AFR"
vcf_query_format+="\t%INFO/gnomAD_exome_AMR"
vcf_query_format+="\t%INFO/gnomAD_exome_ASJ"
vcf_query_format+="\t%INFO/gnomAD_exome_EAS"
vcf_query_format+="\t%INFO/gnomAD_exome_FIN"
vcf_query_format+="\t%INFO/gnomAD_exome_NFE"
vcf_query_format+="\t%INFO/gnomAD_exome_OTH"
vcf_query_format+="[\t%$gt_format]"
vcf_query_format+="\n"
vcf_query_format+="'"

awk_filtering_condition="\$$TMP_VARIANTS_LIST_GNOMAD_EXOME_ALL_COL_IDX < $frequency_ratio"
awk_filtering_condition+=" && \$$TMP_VARIANTS_LIST_GNOMAD_EXOME_AFR_COL_IDX < $frequency_ratio"
awk_filtering_condition+=" && \$$TMP_VARIANTS_LIST_GNOMAD_EXOME_AMR_COL_IDX < $frequency_ratio"
awk_filtering_condition+=" && \$$TMP_VARIANTS_LIST_GNOMAD_EXOME_ASJ_COL_IDX < $frequency_ratio"
awk_filtering_condition+=" && \$$TMP_VARIANTS_LIST_GNOMAD_EXOME_EAS_COL_IDX < $frequency_ratio"
awk_filtering_condition+=" && \$$TMP_VARIANTS_LIST_GNOMAD_EXOME_FIN_COL_IDX < $frequency_ratio"
awk_filtering_condition+=" && \$$TMP_VARIANTS_LIST_GNOMAD_EXOME_NFE_COL_IDX < $frequency_ratio"
awk_filtering_condition+=" && \$$TMP_VARIANTS_LIST_GNOMAD_EXOME_OTH_COL_IDX < $frequency_ratio"

tmp_filtered_variants_list="$working_dir/tmp_filtered_variants_list"
cmd="vcf-query"
cmd+=" -f $vcf_query_format"
cmd+=" $tmp_ta_vcf.gz"
cmd+=" | awk '{ if ($awk_filtering_condition) print \$0 }'"
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
info_msg ">>> extracting tri_snpss <<<"

tmp_get_tri_snps_cmds="$working_dir/tmp_get_tri_snps_cmds"
cmd="awk -F\$'\t'"
cmd+=" '{ printf \"samtools faidx $reference %s:%s-%s | tail -1\n\", \$$TMP_VARIANTS_LIST_CHROM_COL_IDX, \$$TMP_VARIANTS_LIST_POS_COL_IDX-1, \$$TMP_VARIANTS_LIST_POS_COL_IDX+1 }'"
cmd+=" $tmp_filtered_variants_list"
cmd+=" > $tmp_get_tri_snps_cmds"
eval_cmd "$cmd"

tmp_tri_snpss="$working_dir/tmp_tri_snpss"
:>$tmp_tri_snpss
while read cmd; do
    eval "$cmd >> $tmp_tri_snpss"
done <  $tmp_get_tri_snps_cmds

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
info_msg ">>> merging genes, coordinates, tri_snpss, and zygosities together <<<"

tmp_genes_coors_tri_snpss_gtz="$working_dir/tmp_genes_coors_tri_snpss_gtz"
cmd="paste"
cmd+=" $tmp_genes"
cmd+=" $tmp_coors"
cmd+=" $tmp_tri_snpss"
cmd+=" $tmp_gtz"
cmd+=" > $tmp_genes_coors_tri_snpss_gtz"
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
info_msg ">>> merging strands, genes, coordinates, tri_snpss, and zygosities together <<<"

echo -e "#Strand\tGene\tChr\tPos\tRef\tAlt\ttri_snps\t$samples_list" > $out_event_file

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
cmd+=" <( sort -k1,1 $tmp_genes_coors_tri_snpss_gtz )"
cmd+=" >> $out_event_file"
eval_cmd "$cmd"

new_section_txt "counting Ti/Tv"

first_sample_idx=8

function count_subevents {
    events=$1

    events_count=""
    for sample_idx in $(seq 0 $((n_samples-1)))
    do
        sample_id="${array_samples_list[$sample_idx]}"
        # exclude normal tissue if the vcf data are from tumors
        if [[ $sample_id = *"NORMAL"* ]]
        then
            continue
        fi

        event_sample_idx=$(( sample_idx+first_sample_idx ))
        # count number of variants on the transcribed strand
        cmd="grep \"^ts\" $events"
        cmd+=" | cut -f$event_sample_idx"
        cmd+=" | grep -o 1"
        cmd+=" | wc -l"
        ts_count=$( eval $cmd  || true ) 

        # count number of variants on the untranscribed strand
        cmd="grep \"^uts\" $events"
        cmd+=" | cut -f$event_sample_idx"
        cmd+=" | grep -o 1"
        cmd+=" | wc -l"
        uts_count=$( eval $cmd  || true ) 
        events_count+="\t$(( ts_count+uts_count ))\t$ts_count\t$uts_count"
    done
    echo "$events_count"
}

function count_events {
    ref=$1
    alt=$2
    tri_snps=$3

    case "$ref" in
      C)
        rev_ref="G"
        ;;
      T)
        rev_ref="A"
        ;;
      *)
        die "unrecognized ref: $ref"
        ;;
    esac

    case "$alt" in
      A)
        rev_alt="T"
        ;;
      C)
        rev_alt="G"
        ;;
      G)
        rev_alt="C"
        ;;
      T)
        rev_alt="A"
        ;;
      *)
        die "unrecognized alt: $alt"
        ;;
    esac

    case "$tri_snps" in
      ACA)
        rev_cpl_tri_snps="TGT"
        ;;
      ACC)
        rev_cpl_tri_snps="GGT"
        ;;
      ACG)
        rev_cpl_tri_snps="CGT"
        ;;
      ACT)
        rev_cpl_tri_snps="AGT"
        ;;
      CCA)
        rev_cpl_tri_snps="TGG"
        ;;
      CCC)
        rev_cpl_tri_snps="GGG"
        ;;
      CCG)
        rev_cpl_tri_snps="CGG"
        ;;
      CCT)
        rev_cpl_tri_snps="AGG"
        ;;
      GCA)
        rev_cpl_tri_snps="TGC"
        ;;
      GCC)
        rev_cpl_tri_snps="GGC"
        ;;
      GCG)
        rev_cpl_tri_snps="CGC"
        ;;
      GCT)
        rev_cpl_tri_snps="AGC"
        ;;
      TCA)
        rev_cpl_tri_snps="TGA"
        ;;
      TCC)
        rev_cpl_tri_snps="GGA"
        ;;
      TCG)
        rev_cpl_tri_snps="CGA"
        ;;
      TCT)
        rev_cpl_tri_snps="AGA"
        ;;
      ATA)
        rev_cpl_tri_snps="TAT"
        ;;
      ATC)
        rev_cpl_tri_snps="GAT"
        ;;
      ATG)
        rev_cpl_tri_snps="CAT"
        ;;
      ATT)
        rev_cpl_tri_snps="AAT"
        ;;
      CTA)
        rev_cpl_tri_snps="TAG"
        ;;
      CTC)
        rev_cpl_tri_snps="GAG"
        ;;
      CTG)
        rev_cpl_tri_snps="CAG"
        ;;
      CTT)
        rev_cpl_tri_snps="AAG"
        ;;
      GTA)
        rev_cpl_tri_snps="TAC"
        ;;
      GTC)
        rev_cpl_tri_snps="GAC"
        ;;
      GTG)
        rev_cpl_tri_snps="CAC"
        ;;
      GTT)
        rev_cpl_tri_snps="AAC"
        ;;
      TTA)
        rev_cpl_tri_snps="TAA"
        ;;
      TTC)
        rev_cpl_tri_snps="GAA"
        ;;
      TTG)
        rev_cpl_tri_snps="CAA"
        ;;
      TTT)
        rev_cpl_tri_snps="AAA"
        ;;
      *)
        die "unrecognized tri_snps: $tri_snps"
        ;;
    esac

    tmp_events="$working_dir/tmp_events"

    ts_fwd_ref=$ref
    ts_fwd_alt=$alt
    ts_fwd_tri_snps=$tri_snps
    ts_rev_ref=$rev_ref
    ts_rev_alt=$rev_alt
    ts_rev_tri_snps=$rev_cpl_tri_snps

    cmd="awk -F\$'\t'"
    cmd+=" '{ if ((\$1 == \"+\" && \$5 == \"$ts_fwd_ref\" && \$6 == \"$ts_fwd_alt\" && \$7 == \"$ts_fwd_tri_snps\") || (\$1 == \"-\" && \$5 == \"$ts_rev_ref\" && \$6 == \"$ts_rev_alt\" && \$7 == \"$ts_rev_tri_snps\")) printf \"ts %s\n\",\$0 }'"
    cmd+=" $out_event_file"
    cmd+=" > $tmp_events"
    eval "$cmd"

    uts_fwd_ref=$rev_ref
    uts_fwd_alt=$rev_alt
    uts_fwd_tri_snps=$rev_cpl_tri_snps
    uts_rev_ref=$ref
    uts_rev_alt=$alt
    uts_rev_tri_snps=$tri_snps

    cmd="awk -F\$'\t'"
    cmd+=" '{ if ((\$1 == \"+\" && \$5 == \"$uts_fwd_ref\" && \$6 == \"$uts_fwd_alt\" && \$7 == \"$uts_fwd_tri_snps\") || (\$1 == \"-\" && \$5 == \"$uts_rev_ref\" && \$6 == \"$uts_rev_alt\" && \$7 == \"$uts_rev_tri_snps\")) printf \"uts %s\n\",\$0 }'"
    cmd+=" $out_event_file"
    cmd+=" >> $tmp_events"
    eval "$cmd"

    count_subevents $tmp_events
}

# write header
out_header="Substitution Type\tTrinucleotide\tSomatic Mutation Type"
array_samples_list=( $samples_list )  
for sample_idx in $(seq 0 $((n_samples-1)))
do
    sample_id="${array_samples_list[$sample_idx]}"
    # exclude normal tissue if the vcf data are from tumors
    if [[ $sample_id = *"NORMAL"* ]]
    then
        continue
    fi
    out_header+="\t$sample_id"
    out_header+="\t$sample_id"_transcribed
    out_header+="\t$sample_id"_untranscribed
done
echo -e "$out_header" > $out_titv_file

# counting events in alphabetical order of Somatic Mutation Type
for prime3 in A C G T
do
    for ref in C T
    do
        for alt in A C G T
        do
            if [ "$ref" != "$alt" ]
            then
                for prime5 in A C G T
                do
                    substition_type="$ref>$alt"
                    trinucleotide="$prime3$ref$prime5"
                    somatic_mutation_type="$prime3[$substition_type]$prime5"
                    tmp_count="$substition_type\t$trinucleotide\t$somatic_mutation_type"
                    # everything below will be kept in a line according to COSMIC signatures format
                    info_msg "counting somatic muation type: $somatic_mutation_type"
                    tmp_count+="$( count_events $ref $alt $trinucleotide )"
                    echo -e "$tmp_count" >> "$out_titv_file"
                done
            fi
        done
    done
done

info_msg "Done!! Ti/Tv count was exported to $out_titv_file"

new_section_txt "F I N I S H <$script_name>"

