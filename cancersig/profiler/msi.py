import math
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict
from collections import OrderedDict

REPEAT_UNIT_LENGTH_4 = "Repeat_unit_length_4"
REPEAT_UNIT_LENGTH_5 = "Repeat_unit_length_5"

def __is_repeat_sub_unit(seq):
    seq_len = len(seq)
    for unit_len in range(1, math.floor(seq_len/2)+1):
        if seq_len % unit_len > 0:
            continue
        if seq[:unit_len] * int(seq_len/unit_len) == seq:
            return True
    return False

def __get_rotated_sequences(seq):
    tmp_seq = seq
    rotated_seqs = []
    for i in range(1, len(seq)):
        tmp_seq = tmp_seq[-1] + tmp_seq[:-1]
        rotated_seqs.append(tmp_seq)
    return rotated_seqs

def __get_base_combinations(unit_len):
    if unit_len == 0:
        return [""]
    shorter_base_combinations = __get_base_combinations(unit_len-1)
    result = []
    for base in ["A", "C", "G", "T"]:
        result += map(lambda x: base+x, shorter_base_combinations)
    return result

def __init_msi_features_template():
    max_repeat_unit_len = 3
    features_map = defaultdict(list)
    features_hash = {}
    for repeat_unit_len in range(1, max_repeat_unit_len+1):
        for feature in __get_base_combinations(repeat_unit_len):
            dna_seq = Seq(feature, generic_dna)
            reverse_complement = str(dna_seq.reverse_complement())
            if __is_repeat_sub_unit(feature):
                continue
            if reverse_complement in features_map:
                features_map[reverse_complement].append(feature)
                features_hash [feature] = reverse_complement
                continue
            rotated_seq_found = False
            for rotated_seq in __get_rotated_sequences(feature):
                if rotated_seq in features_map:
                    features_map[rotated_seq].append(feature)
                    features_hash [feature] = rotated_seq
                    rotated_seq_found = True
                    break
            if rotated_seq_found:
                continue
            for rotated_seq in __get_rotated_sequences(reverse_complement):
                if rotated_seq in features_map:
                    features_map[rotated_seq].append(feature)
                    features_hash [feature] = rotated_seq
                    rotated_seq_found = True
                    break
            if rotated_seq_found:
                continue
            features_map[feature].append(feature)
            features_hash[feature] = feature

    features_template = OrderedDict()
    for uniq_feature in features_map:
        features_template[uniq_feature] = 0
    features_template[REPEAT_UNIT_LENGTH_4] = 0
    features_template[REPEAT_UNIT_LENGTH_5] = 0
    return features_template

MSI_FEATURES_TEMPLATE = __init_msi_features_template()

#from os.path import basename
#import argparse
#import sys
#import glob
#import datetime
#import subprocess
#import copy
#from os.path import join as join_path
#from Bio.Seq import Seq
#
#MSI_H_CUTOFF = 3.5
#
#class MSIFeaturesExtractByBaseComposition():
#
#    def __init__(self):
#        self.__init_features_template()
#
#    def __info(self, msg=""):
#        msg = datetime.datetime.now().strftime("## [INFO] %Y-%m-%d %H:%M:%S - ") + str(msg)
#        print(msg, file=sys.stderr)
#
#    def __exec_sh(self, cmd, silent=False):
#        self.__info("executing: " + repr(cmd))
#        p = subprocess.Popen(cmd,
#                             shell=True,
#                             stdout=subprocess.PIPE,
#                             stderr=subprocess.PIPE,
#                             )
#        stdout_data, stderr_data = p.communicate()
#        return_code = p.returncode
#        if not silent:
#            print(stdout_data)
#            if return_code:
#                mylogger.throw("Error found during execute command '%s' with error code: %d, %s" % (cmd, return_code, stderr_data))
#            print(stderr_data, file=sys.stderr)
#        elif stderr_data:
#            print(stderr_data, file=sys.stderr)
#        return p, stdout_data
#
#    def __get_base_combinations(self, unit_len):
#        if unit_len == 0:
#            return [""]
#        shorter_base_combinations = self.__get_base_combinations(unit_len-1)
#        result = []
#        for base in ["A", "C", "G", "T"]:
#            result += map(lambda x: base+x, shorter_base_combinations)
#        return result
#
#    def __is_repeat_sub_unit(self, seq):
#        seq_len = len(seq)
#        for unit_len in range(1, math.floor(seq_len/2)+1):
#            if seq_len % unit_len > 0:
#                continue
#            if seq[:unit_len] * int(seq_len/unit_len) == seq:
#                return True
#        return False
#
#    def __get_rotated_sequences(self, seq):
#        tmp_seq = seq
#        rotated_seqs = []
#        for i in range(1, len(seq)):
#            tmp_seq = tmp_seq[-1] + tmp_seq[:-1]
#            rotated_seqs.append(tmp_seq)
#        return rotated_seqs
#
#    def __init_features_template(self):
#        max_repeat_unit_len = 3
#        features_map = defaultdict(list)
#        self.__features_hash = {}
#        for repeat_unit_len in range(1, max_repeat_unit_len+1):
#            for feature in self.__get_base_combinations(repeat_unit_len):
#                dna_seq = Seq(feature, generic_dna)
#                reverse_complement = str(dna_seq.reverse_complement())
#                if self.__is_repeat_sub_unit(feature):
#                    continue
#                if reverse_complement in features_map:
#                    features_map[reverse_complement].append(feature)
#                    self.__features_hash [feature] = reverse_complement
#                    continue
#                rotated_seq_found = False
#                for rotated_seq in self.__get_rotated_sequences(feature):
#                    if rotated_seq in features_map:
#                        features_map[rotated_seq].append(feature)
#                        self.__features_hash [feature] = rotated_seq
#                        rotated_seq_found = True
#                        break
#                if rotated_seq_found:
#                    continue
#                for rotated_seq in self.__get_rotated_sequences(reverse_complement):
#                    if rotated_seq in features_map:
#                        features_map[rotated_seq].append(feature)
#                        self.__features_hash [feature] = rotated_seq
#                        rotated_seq_found = True
#                        break
#                if rotated_seq_found:
#                    continue
#                features_map[feature].append(feature)
#                self.__features_hash[feature] = feature
#
#        self.__features_template = OrderedDict()
#        for uniq_feature in features_map:
#            self.__features_template[uniq_feature] = 0
#        self.__features_template[REPEAT_UNIT_LENGTH_4] = 0
#        self.__features_template[REPEAT_UNIT_LENGTH_5] = 0
#
##    def __extract_zeros_features(self, sample_id, output_file):
##        with open(output_file, "w") as f_o:
##            header = "profile_id"
##            header += "\t" + "\t".join(self.__features_template.keys())
##            f_o.write(header+"\n")
##            content = sample_id
##            content += "\t0.000001" * len(self.__features_template.keys())
##            f_o.write(content+"\n")
##
#    def __extract_features(self, raw_msi_somatic, sample_id, output_file, zero_features=False):
#        with open(output_file, "w") as f_o:
#            header = "variant type"
#            header += "\tvariant subgroup"
#            header += "\tfeature_id"
#            header += "\t" + sample_id
#            f_o.write(header+"\n")
#            features_count = copy.copy(self.__features_template)
#            if not zero_features:
#                cmd = "cut -f 5 " + raw_msi_somatic
#                cmd += " | sort"
#                cmd += " | uniq -c"
#                p, stdout_data = self.__exec_sh(cmd, silent=True)
#                unstable_loci_counts = list(map(lambda x: x.strip(), stdout_data.decode('utf-8').strip().split("\n")))
#                total_count = 0
#                for unstable_loci_count in unstable_loci_counts:
#                    _count, _repeat = unstable_loci_count.split()
#                    _count = float(_count)
#                    total_count += _count
#                    if len(_repeat) == 4:
#                        features_count[REPEAT_UNIT_LENGTH_4] += _count
#                        continue
#                    if len(_repeat) == 5:
#                        features_count[REPEAT_UNIT_LENGTH_5] += _count
#                        continue
#                    feature = self.__features_hash[_repeat]
#                    features_count[feature] += _count
#            for feature in features_count:
#                feature_count = features_count[feature]
#                content = "MSI"
#                if feature == "Repeat_unit_length_4":
#                    content += "\tLength_4"
#                elif feature == "Repeat_unit_length_5":
#                    content += "\tLength_5"
#                else:
#                    content += "\t" + feature
#                content += "\t" + feature
#                if zero_features:
#                    content += "\t0.00000001"
#                else:
#                    content += "\t" + "{:.8f}".format(feature_count/total_count+0.00000001)
#                f_o.write(content+"\n")
##            content += "\t" + "\t".join(map(lambda x: "{:8.6f}".format(x/total_count+0.000001), features_count.values()))
#
#    def extract_features(self,
#                         raw_msi_report,
#                         raw_msi_somatic,
#                         sample_id,
#                         output_file,
#                         ):
#        # MSI-H is one of the main factor to decide how the features look like
#        with open(raw_msi_report, "r") as f_r:
#            f_r.readline()
#            msi_score = float(f_r.readline().strip().split()[2])
#            if msi_score < MSI_H_CUTOFF:
#                self.__extract_features(raw_msi_somatic, sample_id, output_file, zero_features=True)
#            else:
#                self.__extract_features(raw_msi_somatic, sample_id, output_file)
