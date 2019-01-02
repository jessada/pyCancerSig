#!/bin/bash 

set -e
set -u
set -o pipefail

if [[ -z ${PYCANCERSIG_SCRIPTS_DIR:-} ]]
then
    script_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
    source "$script_dir/../pyCancerSig_envs.sh"
fi

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-i {dir}          input directory with reference fasta chromosome by chromosome (required)
-o {file}         output list of microsatellite (required)
-h                this help
EOF
)

# parse option
while getopts ":i:o:h" OPTION; do
  case "$OPTION" in
    i)
      input_refs_dir="$OPTARG"
      ;;
    o)
      out_list="$OPTARG"
      ;;
    h)
      echo >&2 "$usage"
      ;;
    *)
      die "unrecognized option (-$OPTION) from executing: $0 $@"
      ;;
  esac
done

[ ! -z $input_refs_dir ] || die "input references directory is missing (-i)"
[ -d $input_refs_dir ] || die "$input_refs_dir is not found"
[ ! -z $out_list ] || die "output file name is missing (-o)"

time_stamp=$( date )

cd $PYCANCERSIG_SCRIPTS_DIR
revision_no=`git rev-list HEAD | wc -l`
revision_code=`git rev-parse HEAD`
cd - > /dev/null

working_dir=`mktemp -d`

## ****************************************  display configuration  ****************************************
## display required configuration
new_section_txt "S T A R T <$script_name>"
info_msg
info_msg "description"
info_msg "  Using msisensor to generate possilble STR loci chromosome by chromosome"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "References directory (-i)" "$input_refs_dir"
display_param "Output list of microsatellite (-o)" "$out_list"

# ****************************************  executing  ****************************************
 
new_section_txt "Scaning reference genome"

ms_list_temp="$working_dir/ms_list_temp"

echo -e "chromosome\tlocation\trepeat_unit_length\trepeat_unit_binary\trepeat_times\tleft_flank_binary\tright_flank_binary\trepeat_unit_bases\tleft_flank_bases\tright_flank_bases" > "$out_list"

for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
#    ls -altr "$input_refs_dir/chr$chrom.fa"
    cmd="msisensor scan"
    cmd+=" -d $input_refs_dir/chr$chrom.fa"
    cmd+=" -o $ms_list_temp"
    eval_cmd "$cmd"

    cmd="grep -v \"^chromosome\" $ms_list_temp"
    cmd+=" >> $out_list"
    eval_cmd "$cmd"
done


new_section_txt "F I N I S H <$script_name>"
