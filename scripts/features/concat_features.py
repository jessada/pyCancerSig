#!/usr/bin/env python
import argparse
import datetime
import sys
from parse_n_join_features import create_features_dict

VARIANT_TYPE = "variant type"
VARIANT_INFO = "variant info"

def __log_info(msg=""):
    msg = datetime.datetime.now().strftime("## [INFO] %Y-%m-%d %H:%M:%S - ") + msg
    print(msg, file=sys.stderr)

def main(features_files,
         output_filename,
         ):
    __log_info()
    __log_info("features files: " + str(features_files))
    __log_info("output file name: " + output_filename)

    # load features
    sample_profiles = {}
    for features_file in features_files:
        sample_profile = {}
        with open(features_file) as f_f:
            header = f_f.readline().strip().split("\t")
            sample_id = header[3]
            for line in f_f:
                content = line.strip().split("\t")
                feature_id = content[2]
                weight = content[3]
                sample_profile[feature_id] = weight
        sample_profiles[sample_id] = sample_profile
    

    features_dict = create_features_dict()
    sample_ids = sample_profiles.keys()
    print(sample_ids)
    __log_info()
    __log_info("Writing output profiles . . .")
    with open(output_filename, "w") as f_o:
        header = VARIANT_TYPE
        header += "\t" + VARIANT_INFO
        header += "\tfeature_id"
        header += "\t" + "\t".join(sample_ids)
        f_o.write(header+"\n")
        for feature_id in features_dict:
            content = features_dict[feature_id][VARIANT_TYPE]
            content += "\t" + features_dict[feature_id][VARIANT_INFO]
            content += "\t" + feature_id
            for sample_id in sample_ids:
                content += "\t" + sample_profiles[sample_id][feature_id]
            f_o.write(content+"\n")

    __log_info()
    __log_info("Done !!!")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(""" cancat features of all the input samples""")

    parser.add_argument('--features_files',
                        dest='features_files',
                        required=True,
                        )
    parser.add_argument('--output_filename',
                        dest='output_filename',
                        required=True,
                        )
    args = parser.parse_args()
    parsed_features_files=list(map(lambda x: x.strip(','), args.features_files.strip(']').strip('[').split()))
    main(features_files=parsed_features_files,
         output_filename=args.output_filename,
         )
