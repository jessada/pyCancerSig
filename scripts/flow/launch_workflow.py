#!/usr/bin/env python

import argparse
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--config',
                    metavar="FILE",
                    dest='raw_msisensor_report',
                    help='Nextflow workflow.config (-C) file',
#                    required=True,
                    )
parser.add_argument('--bam',
                    metavar="DIR",
                    dest='bam',
                    help='A directory containing sub-directories with bam files. Each sub-directory is for one sample, with both tumor and germline bam.',
                    required=True,
                    )
args = parser.parse_args()


path = os.path.dirname(sys.argv[0])
cmd = os.path.join(path, "launch_full_workflow.sh")
os.system(cmd)
