import copy
import math
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict
from collections import OrderedDict
from cancersig.utils import exec_sh
from cancersig.template import pyCancerSigBase

REPEAT_UNIT_LENGTH_4 = "Repeat_unit_length_4"
REPEAT_UNIT_LENGTH_5 = "Repeat_unit_length_5"

MSI_H_CUTOFF = 3.5

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
    return features_template, features_hash

MSI_FEATURES_TEMPLATE, MSI_FEATURES_HASH = __init_msi_features_template()

class MSIProfiler(pyCancerSigBase):

    def __init__(self, *args, **kwargs):
        super(MSIProfiler, self).__init__(*args, **kwargs)

    def __extract_features(self,
                           raw_msisensor_out_somatic,
                           sample_id,
                           output_file,
                           msi_positive,
                           ):
        with open(output_file, "w") as f_o:
            header = "variant type"
            header += "\tvariant subgroup"
            header += "\tfeature_id"
            header += "\t" + sample_id
            f_o.write(header+"\n")
            features_count = copy.copy(MSI_FEATURES_TEMPLATE)
            if msi_positive:
                cmd = "cut -f 5 " + raw_msisensor_out_somatic
                cmd += " | sort"
                cmd += " | uniq -c"
                p, stdout_data = exec_sh(cmd, silent=True)
                unstable_loci_counts = list(map(lambda x: x.strip(), stdout_data.decode('utf-8').strip().split("\n")))
                total_count = 0
                for unstable_loci_count in unstable_loci_counts:
                    _count, _repeat = unstable_loci_count.split()
                    _count = float(_count)
                    total_count += _count
                    if len(_repeat) == 4:
                        features_count[REPEAT_UNIT_LENGTH_4] += _count
                        continue
                    if len(_repeat) == 5:
                        features_count[REPEAT_UNIT_LENGTH_5] += _count
                        continue
                    feature = MSI_FEATURES_HASH[_repeat]
                    features_count[feature] += _count
            for feature in features_count:
                feature_count = features_count[feature]
                content = "MSI"
                if feature == "Repeat_unit_length_4":
                    content += "\tLength_4"
                elif feature == "Repeat_unit_length_5":
                    content += "\tLength_5"
                else:
                    content += "\t" + feature
                content += "\t" + feature
                if msi_positive:
                    content += "\t" + "{:.8f}".format(feature_count/total_count+0.00000001)
                else:
                    content += "\t" + "{:.8f}".format(0.00000001)
                f_o.write(content+"\n")

    def extract_features(self,
                         raw_msisensor_out,
                         raw_msisensor_out_somatic,
                         sample_id,
                         output_file,
                         ):
        # MSI-H is one of the main factor to decide how the features look like
        with open(raw_msisensor_out, "r") as f_r:
            f_r.readline()
            msi_score = float(f_r.readline().strip().split()[2])
            msi_positive = msi_score > MSI_H_CUTOFF
            self.__extract_features(raw_msisensor_out_somatic,
                                    sample_id,
                                    output_file,
                                    msi_positive=msi_positive,
                                    )
