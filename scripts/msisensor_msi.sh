#!/bin/bash 

set -e
set -u
set -o pipefail

if [[ -z ${PYCANCERSIG_DIR:-} ]]
then
    script_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
    source "$script_dir/pyCancerSig_envs.sh"
fi

script_name=$(basename $0)
params="$@"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-d {file}         homopolymer and microsates file (required)
-n {file}         normal bam file (required)
-t {file}         tumor bam file (required)
-o {file_prefix}  output distribution file prefix (required)
-h                this help
EOF
)

# parse option
while getopts ":d:n:t:o:h" OPTION; do
  case "$OPTION" in
    d)
      ms_file="$OPTARG"
      ;;
    n)
      normal_bam_file="$OPTARG"
      ;;
    t)
      tumor_bam_file="$OPTARG"
      ;;
    o)
      out_prefix="$OPTARG"
      ;;
    h)
      echo >&2 "$usage"
      ;;
    *)
      die "unrecognized option (-$OPTION) from executing: $0 $@"
      ;;
  esac
done

[ ! -z $ms_file ] || die "input microsates is missing (-d)"
[ -f $ms_file ] || die "$ms_file is not found"
[ ! -z $normal_bam_file ] || die "input normal bam file is missing (-d)"
[ -f $normal_bam_file ] || die "$normal_bam_file is not found"
[ ! -z $tumor_bam_file ] || die "input tumor bam file is missing (-d)"
[ -f $tumor_bam_file ] || die "$tumor_bam_file is not found"
[ ! -z $out_prefix ] || die "output file prefix is missing (-o)"

time_stamp=$( date )

cd $PYCANCERSIG_DIR
revision_no=`git rev-list HEAD | wc -l`
revision_code=`git rev-parse HEAD`
cd - > /dev/null

working_dir=`mktemp -d`

## ****************************************  display configuration  ****************************************
## display required configuration
new_section_txt "S T A R T <$script_name>"
info_msg
info_msg "description"
info_msg "  An encapsulation of 'msisensor msi' for slurm"
info_msg
info_msg "version and script configuration"
display_param "revision no" "$revision_no"
display_param "revision code" "$revision_code"
display_param "parameters" "$params"
display_param "time stamp" "$time_stamp"
info_msg
info_msg "overall configuration"
display_param "Microsates list (-d)" "$ms_file"
display_param "Normal bam file (-n)" "$normal_bam_file"
display_param "Tumor bam file (-t)" "$tumor_bam_file"
display_param "Output file prefix (-o)" "$out_prefix"

# ****************************************  executing  ****************************************
 
new_section_txt "Running msisensor msi"

cmd="msisensor msi"
cmd+=" -d $ms_file"
cmd+=" -n $normal_bam_file"
cmd+=" -t $tumor_bam_file"
cmd+=" -o $out_prefix"
eval_cmd "$cmd"

new_section_txt "F I N I S H <$script_name>"
