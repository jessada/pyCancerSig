#!/bin/bash 

set -e
set -u
set -o pipefail

script_dir=""
source $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../pyCancerSig_envs.sh

echo $PYCANCERSIG_SCRIPTS_DIR
#module load Nextflow



