#!/usr/bin/env python

import argparse
import os
import fnmatch
import datetime
import sys
from os.path import join as join_path
from collections import OrderedDict

VARIANT_TYPE = "variant type"
VARIANT_INFO = "variant info"

SV_BND_TYPE = "Translocation"
SV_INV_TYPE = "Inversion"
SV_DEL_TYPE = "Deletion"
SV_DUP_TYPE = "Duplication"

MSI_TYPE = "MSI"

LOG10_2_3_INFO = "log10_2_3"
LOG10_3_4_INFO = "log10_3_4"
LOG10_4_5_INFO = "log10_4_5"
LOG10_5_6_INFO = "log10_5_6"
LOG10_6_7_INFO = "log10_6_7"
LOG10_7_8_INFO = "log10_7_8"
LOG10_8_9_INFO = "log10_8_9"
LOG10_9UP_INFO = "log10_9up"

def create_features_dict():
    # generate Ti/TV features
    features_dict = OrderedDict()
    for prime_5 in ["A", "C", "G", "T"]:
        for ref in ["C", "T"]:
            for alt in ["A", "C", "G", "T"]:
                for prime_3 in ["A", "C", "G", "T"]:
                    if ref != alt:
                        variant_type = ref + ">" + alt
                        variant_info = prime_5 + ref + prime_3
                        feature_id = prime_5 + "[" + variant_type + "]" + prime_3
                        features_dict[feature_id] = {VARIANT_TYPE: variant_type,
                                                     VARIANT_INFO: variant_info,
                                                     }
    # and manually generate SV
    features_dict["BND_log10_2_3"] = {VARIANT_TYPE: SV_BND_TYPE,
                                      VARIANT_INFO: LOG10_2_3_INFO,
                                      }
    features_dict["BND_log10_3_4"] = {VARIANT_TYPE: SV_BND_TYPE,
                                      VARIANT_INFO: LOG10_3_4_INFO,
                                      }
    features_dict["BND_log10_4_5"] = {VARIANT_TYPE: SV_BND_TYPE,
                                      VARIANT_INFO: LOG10_4_5_INFO,
                                      }
    features_dict["BND_log10_5_6"] = {VARIANT_TYPE: SV_BND_TYPE,
                                      VARIANT_INFO: LOG10_5_6_INFO,
                                      }
    features_dict["BND_log10_6_7"] = {VARIANT_TYPE: SV_BND_TYPE,
                                      VARIANT_INFO: LOG10_6_7_INFO,
                                      }
    features_dict["BND_log10_7_8"] = {VARIANT_TYPE: SV_BND_TYPE,
                                      VARIANT_INFO: LOG10_7_8_INFO,
                                      }
    features_dict["BND_log10_8_9"] = {VARIANT_TYPE: SV_BND_TYPE,
                                      VARIANT_INFO: LOG10_8_9_INFO,
                                      }
    features_dict["BND_log10_9up"] = {VARIANT_TYPE: SV_BND_TYPE,
                                      VARIANT_INFO: LOG10_9UP_INFO,
                                      }
    features_dict["INV_log10_2_3"] = {VARIANT_TYPE: SV_INV_TYPE,
                                      VARIANT_INFO: LOG10_2_3_INFO,
                                      }
    features_dict["INV_log10_3_4"] = {VARIANT_TYPE: SV_INV_TYPE,
                                      VARIANT_INFO: LOG10_3_4_INFO,
                                      }
    features_dict["INV_log10_4_5"] = {VARIANT_TYPE: SV_INV_TYPE,
                                      VARIANT_INFO: LOG10_4_5_INFO,
                                      }
    features_dict["INV_log10_5_6"] = {VARIANT_TYPE: SV_INV_TYPE,
                                      VARIANT_INFO: LOG10_5_6_INFO,
                                      }
    features_dict["INV_log10_6_7"] = {VARIANT_TYPE: SV_INV_TYPE,
                                      VARIANT_INFO: LOG10_6_7_INFO,
                                      }
    features_dict["INV_log10_7_8"] = {VARIANT_TYPE: SV_INV_TYPE,
                                      VARIANT_INFO: LOG10_7_8_INFO,
                                      }
    features_dict["INV_log10_8_9"] = {VARIANT_TYPE: SV_INV_TYPE,
                                      VARIANT_INFO: LOG10_8_9_INFO,
                                      }
    features_dict["INV_log10_9up"] = {VARIANT_TYPE: SV_INV_TYPE,
                                      VARIANT_INFO: LOG10_9UP_INFO,
                                      }
    features_dict["DEL_log10_2_3"] = {VARIANT_TYPE: SV_DEL_TYPE,
                                      VARIANT_INFO: LOG10_2_3_INFO,
                                      }
    features_dict["DEL_log10_3_4"] = {VARIANT_TYPE: SV_DEL_TYPE,
                                      VARIANT_INFO: LOG10_3_4_INFO,
                                      }
    features_dict["DEL_log10_4_5"] = {VARIANT_TYPE: SV_DEL_TYPE,
                                      VARIANT_INFO: LOG10_4_5_INFO,
                                      }
    features_dict["DEL_log10_5_6"] = {VARIANT_TYPE: SV_DEL_TYPE,
                                      VARIANT_INFO: LOG10_5_6_INFO,
                                      }
    features_dict["DEL_log10_6_7"] = {VARIANT_TYPE: SV_DEL_TYPE,
                                      VARIANT_INFO: LOG10_6_7_INFO,
                                      }
    features_dict["DEL_log10_7_8"] = {VARIANT_TYPE: SV_DEL_TYPE,
                                      VARIANT_INFO: LOG10_7_8_INFO,
                                      }
    features_dict["DEL_log10_8_9"] = {VARIANT_TYPE: SV_DEL_TYPE,
                                      VARIANT_INFO: LOG10_8_9_INFO,
                                      }
    features_dict["DEL_log10_9up"] = {VARIANT_TYPE: SV_DEL_TYPE,
                                      VARIANT_INFO: LOG10_9UP_INFO,
                                      }
    features_dict["DUP_log10_2_3"] = {VARIANT_TYPE: SV_DUP_TYPE,
                                      VARIANT_INFO: LOG10_2_3_INFO,
                                      }
    features_dict["DUP_log10_3_4"] = {VARIANT_TYPE: SV_DUP_TYPE,
                                      VARIANT_INFO: LOG10_3_4_INFO,
                                      }
    features_dict["DUP_log10_4_5"] = {VARIANT_TYPE: SV_DUP_TYPE,
                                      VARIANT_INFO: LOG10_4_5_INFO,
                                      }
    features_dict["DUP_log10_5_6"] = {VARIANT_TYPE: SV_DUP_TYPE,
                                      VARIANT_INFO: LOG10_5_6_INFO,
                                      }
    features_dict["DUP_log10_6_7"] = {VARIANT_TYPE: SV_DUP_TYPE,
                                      VARIANT_INFO: LOG10_6_7_INFO,
                                      }
    features_dict["DUP_log10_7_8"] = {VARIANT_TYPE: SV_DUP_TYPE,
                                      VARIANT_INFO: LOG10_7_8_INFO,
                                      }
    features_dict["DUP_log10_8_9"] = {VARIANT_TYPE: SV_DUP_TYPE,
                                      VARIANT_INFO: LOG10_8_9_INFO,
                                      }
    features_dict["DUP_log10_9up"] = {VARIANT_TYPE: SV_DUP_TYPE,
                                      VARIANT_INFO: LOG10_9UP_INFO,
                                      }

    # and manually generate MSI
    features_dict["A"] = {VARIANT_TYPE: MSI_TYPE,
                          VARIANT_INFO: "A",
                          }
    features_dict["C"] = {VARIANT_TYPE: MSI_TYPE,
                          VARIANT_INFO: "C",
                          }
    features_dict["AC"] = {VARIANT_TYPE: MSI_TYPE,
                           VARIANT_INFO: "AC",
                           }
    features_dict["AG"] = {VARIANT_TYPE: MSI_TYPE,
                           VARIANT_INFO: "AG",
                           }
    features_dict["AT"] = {VARIANT_TYPE: MSI_TYPE,
                           VARIANT_INFO: "AT",
                           }
    features_dict["CG"] = {VARIANT_TYPE: MSI_TYPE,
                           VARIANT_INFO: "CG",
                           }
    features_dict["AAC"] = {VARIANT_TYPE: MSI_TYPE,
                           VARIANT_INFO: "AAC",
                           }
    features_dict["AAG"] = {VARIANT_TYPE: MSI_TYPE,
                           VARIANT_INFO: "AAG",
                           }
    features_dict["AAT"] = {VARIANT_TYPE: MSI_TYPE,
                           VARIANT_INFO: "AAT",
                           }
    features_dict["ACC"] = {VARIANT_TYPE: MSI_TYPE,
                           VARIANT_INFO: "ACC",
                           }
    features_dict["ACG"] = {VARIANT_TYPE: MSI_TYPE,
                           VARIANT_INFO: "ACG",
                           }
    features_dict["ACT"] = {VARIANT_TYPE: MSI_TYPE,
                           VARIANT_INFO: "ACT",
                           }
    features_dict["AGC"] = {VARIANT_TYPE: MSI_TYPE,
                           VARIANT_INFO: "AGC",
                           }
    features_dict["AGG"] = {VARIANT_TYPE: MSI_TYPE,
                           VARIANT_INFO: "AGG",
                           }
    features_dict["ATC"] = {VARIANT_TYPE: MSI_TYPE,
                           VARIANT_INFO: "ATC",
                           }
    features_dict["CCG"] = {VARIANT_TYPE: MSI_TYPE,
                            VARIANT_INFO: "CCG",
                            }
    features_dict["Repeat_unit_length_4"] = {VARIANT_TYPE: MSI_TYPE,
                                             VARIANT_INFO: "Length_4",
                                             }
    features_dict["Repeat_unit_length_5"] = {VARIANT_TYPE: MSI_TYPE,
                                             VARIANT_INFO: "Length_5",
                                             }
    return features_dict

def load_titv_features(titv_file):
    __log_info(titv_file)
    sum_titv_count = 0
    features_count = {}
    with open(titv_file) as f_t:
        f_t.readline()
        for line in f_t:
            content = line.strip().split("\t")
            feature_id = content[2]
            titv_count = float(content[3])
            features_count[feature_id] = titv_count
            sum_titv_count += titv_count
    return features_count, sum_titv_count

def load_sv_features(sv_file):
    sum_sv_count = 0
    features_count = {}
    with open(sv_file) as f_v:
        header = f_v.readline().strip().split("\t")
        content = f_v.readline().strip().split("\t")
        for idx in range(1, len(header)):
            feature_id = header[idx]
            sv_count = float(content[idx])
            features_count[feature_id] = sv_count
            sum_sv_count += sv_count
    return features_count, sum_sv_count

def load_msi_features(msi_file):
    sum_msi_count = 0
    features_count = {}
    with open(msi_file) as f_m:
        header = f_m.readline().strip().split("\t")
        content = f_m.readline().strip().split("\t")
        for idx in range(1, len(header)):
            feature_id = header[idx]
            msi_count = float(content[idx])
            features_count[feature_id] = msi_count
            sum_msi_count += msi_count
    return features_count, sum_msi_count

def __log_info(msg=""):
    msg = datetime.datetime.now().strftime("## [INFO] %Y-%m-%d %H:%M:%S - ") + msg
    print(msg, file=sys.stderr)

def main(sample_id,
         msi_features,
         titv_features,
         sv_features,
         output_filename,
         ):
    __log_info()
    __log_info("Sample ID: " + sample_id)
    __log_info("MSI features: " + msi_features)
    __log_info("Ti/Tv features: " + titv_features)
    __log_info("SV features: " + sv_features)
    __log_info("output file name: " + output_filename)

    sample_profile = {}

    features_dict = create_features_dict()

    # load and parse Ti/Tv features
    features_count, sum_titv_count = load_titv_features(titv_features)
    for feature_id in features_count:
        titv_count = features_count[feature_id]
        sample_profile[feature_id] = (titv_count*7)/(sum_titv_count*10)

    # load and parse SV features
    features_count, sum_sv_count = load_sv_features(sv_features)
    for feature_id in features_count:
        sv_count = features_count[feature_id]
        sample_profile[feature_id] = (sv_count*2)/(sum_sv_count*10)

    # load and parse MSI features
    features_count, sum_msi_count = load_msi_features(msi_features)
    for feature_id in features_count:
        msi_count = features_count[feature_id]
        if sum_msi_count > 0.1:
            sample_profile[feature_id] = (msi_count*1)/(sum_msi_count*10)
        else:
            sample_profile[feature_id] = msi_count

    __log_info()
    __log_info("Writing output profiles . . .")
    with open(output_filename, "w") as f_o:
        header = VARIANT_TYPE
        header += "\t" + VARIANT_INFO
        header += "\tfeature_id"
        header += "\t" + sample_id
        f_o.write(header+"\n")
        for feature_id in features_dict:
            content = features_dict[feature_id][VARIANT_TYPE]
            content += "\t" + features_dict[feature_id][VARIANT_INFO]
            content += "\t" + feature_id
            content += "\t{:8.6f}".format(sample_profile[feature_id])
            f_o.write(content+"\n")

    __log_info()
    __log_info("Done !!!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser("""join features from FindSV, MSI and TiTv """)

    parser.add_argument('--sample_id',
                        dest='sample_id',
                        help='sample id to be written in the output file',
                        required=True,
                        )
    parser.add_argument('--msi_features',
                        dest='msi_features',
                        help='A file with MSI profiles',
                        required=True,
                        )
    parser.add_argument('--titv_features',
                        dest='titv_features',
                        help='A file with Ti/Tv profiles',
                        required=True,
                        )
    parser.add_argument('--sv_features',
                        dest='sv_features',
                        help='A file with SV profiles',
                        required=True,
                        )
    parser.add_argument('--output_filename',
                        dest='output_filename',
                        help='output file name',
                        required=True,
                        )
    args = parser.parse_args()
    main(sample_id=args.sample_id,
         msi_features=args.msi_features,
         titv_features=args.titv_features,
         sv_features=args.sv_features,
         output_filename=args.output_filename,
         )
