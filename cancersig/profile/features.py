import math
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from collections import OrderedDict
from collections import defaultdict

PROFILE_TYPE_SNV = "SNV"
PROFILE_TYPE_SV = "SV"
PROFILE_TYPE_MSI = "MSI"
PROFILE_WEIGHTS = OrderedDict()
PROFILE_WEIGHTS[PROFILE_TYPE_SNV] = 7
PROFILE_WEIGHTS[PROFILE_TYPE_SV] = 2
PROFILE_WEIGHTS[PROFILE_TYPE_MSI] = 1
PROFILE_TYPES = list(PROFILE_WEIGHTS.keys())

VARIANT_TYPE = "variant type"
VARIANT_SUBGROUP = "variant subgroup"
FEATURE_ID = "feature id"
FEATURE_QUANTITY = "feature quantity"

VAR_TYPE_NA = "NA"
SMALL_QUANTITY = 0.00000001

def __init_snv_features_template():
    features_hash = {}
    for ref in ["A", "C", "G", "T"]:
        features_hash[ref] = {}
        for alt in ["A", "C", "G", "T"]:
            if ref != alt:
                features_hash[ref][alt] = {}
    features_template = OrderedDict()
    for prime_5 in ["A", "C", "G", "T"]:
        for ref in ["C", "T"]:
            for alt in ["A", "C", "G", "T"]:
                for prime_3 in ["A", "C", "G", "T"]:
                    if ref != alt:
                        # forward strand
                        var_type = ref + ">" + alt
                        var_subgroup= prime_5 + ref + prime_3
                        feature_id = prime_5 + "[" + var_type + "]" + prime_3
                        features_template[feature_id] = {}
                        features_template[feature_id][VARIANT_TYPE] = var_type
                        features_template[feature_id][VARIANT_SUBGROUP] = var_subgroup
                        features_template[feature_id][FEATURE_QUANTITY] = 0
                        features_hash[ref][alt][var_subgroup] = feature_id

                        # reverse strand
                        if ref == "C":
                            rev_ref = "G"
                        elif ref == "T":
                            rev_ref = "A"
                        if alt == "A":
                            rev_alt = "T"
                        elif alt == "C":
                            rev_alt = "G"
                        elif alt == "G":
                            rev_alt = "C"
                        elif alt == "T":
                            rev_alt = "A"
                        features_hash[rev_ref][rev_alt][str(Seq(var_subgroup).reverse_complement())] = feature_id
    return features_template, features_hash
        
SNV_FEATURES_TEMPLATE, SNV_FEATURES_HASH = __init_snv_features_template()

# define SV constant
SV_EVENT_TYPE_LIST = []
SV_EVENT_TYPE_BND = "BND"
SV_VAR_TYPE_BND = "Translocation"
SV_EVENT_TYPE_LIST.append(SV_EVENT_TYPE_BND)
SV_EVENT_TYPE_DEL = "DEL"
SV_VAR_TYPE_DEL = "Deletion"
SV_EVENT_TYPE_LIST.append(SV_EVENT_TYPE_DEL)
SV_EVENT_TYPE_DUP = "DUP"
SV_VAR_TYPE_DUP = "Duplication"
SV_EVENT_TYPE_LIST.append(SV_EVENT_TYPE_DUP)
SV_EVENT_TYPE_INV = "INV"
SV_VAR_TYPE_INV = "Inversion"
SV_EVENT_TYPE_LIST.append(SV_EVENT_TYPE_INV)
SV_VAR_SUBGROUP_WHOLE_CHROM = "whole_chrom"
SV_VAR_SUBGROUP_INTER_CHROM = "inter_chrom"
SV_LEN_LOG10LIST = []
SV_LEN_LOG10_2_3 = "log10_2_3"
SV_LEN_LOG10LIST.append(SV_LEN_LOG10_2_3)
SV_LEN_LOG10_3_4 = "log10_3_4"
SV_LEN_LOG10LIST.append(SV_LEN_LOG10_3_4)
SV_LEN_LOG10_4_5 = "log10_4_5"
SV_LEN_LOG10LIST.append(SV_LEN_LOG10_4_5)
SV_LEN_LOG10_5_6 = "log10_5_6"
SV_LEN_LOG10LIST.append(SV_LEN_LOG10_5_6)
SV_LEN_LOG10_6_7 = "log10_6_7"
SV_LEN_LOG10LIST.append(SV_LEN_LOG10_6_7)
SV_LEN_LOG10_7_8 = "log10_7_8"
SV_LEN_LOG10LIST.append(SV_LEN_LOG10_7_8)
SV_LEN_LOG10_8_9 = "log10_8_9"
SV_LEN_LOG10LIST.append(SV_LEN_LOG10_8_9)
SV_LEN_LOG10_9up = "log10_9up"
SV_LEN_LOG10LIST.append(SV_LEN_LOG10_9up)

def __init_sv_features_template():
    features_template = OrderedDict()
    features_hash = {}
    for event_type in SV_EVENT_TYPE_LIST:
        var_type = VAR_TYPE_NA
        features_hash[event_type] = {}
        if event_type == SV_EVENT_TYPE_BND:
            var_type = SV_VAR_TYPE_BND
        if event_type == SV_EVENT_TYPE_DEL:
            var_type = SV_VAR_TYPE_DEL
        if event_type == SV_EVENT_TYPE_DUP:
            var_type = SV_VAR_TYPE_DUP
        if event_type == SV_EVENT_TYPE_INV:
            var_type = SV_VAR_TYPE_INV
        for len_log10 in SV_LEN_LOG10LIST:
            if (event_type == SV_EVENT_TYPE_BND) and (len_log10 == SV_LEN_LOG10_9up):
                var_subgroup = SV_VAR_SUBGROUP_INTER_CHROM
            elif (event_type == SV_EVENT_TYPE_DEL) and (len_log10 == SV_LEN_LOG10_9up):
                var_subgroup = SV_VAR_SUBGROUP_WHOLE_CHROM
            elif (event_type == SV_EVENT_TYPE_DUP) and (len_log10 == SV_LEN_LOG10_9up):
                var_subgroup = SV_VAR_SUBGROUP_WHOLE_CHROM
            elif (event_type == SV_EVENT_TYPE_INV) and (len_log10 == SV_LEN_LOG10_9up):
                var_subgroup = SV_VAR_SUBGROUP_INTER_CHROM
            else:
                var_subgroup = len_log10
            feature_id = event_type + "_" + var_subgroup
            features_template[feature_id] = {}
            features_template[feature_id][VARIANT_TYPE] = var_type
            features_template[feature_id][VARIANT_SUBGROUP] = var_subgroup
            features_template[feature_id][FEATURE_QUANTITY] = 0
            features_hash[event_type][len_log10] = feature_id
    return features_template, features_hash
        
SV_FEATURES_TEMPLATE, SV_FEATURES_HASH = __init_sv_features_template()

# define MSI constant
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
    features_map[REPEAT_UNIT_LENGTH_4] = []
    features_map[REPEAT_UNIT_LENGTH_5] = []

    features_template = OrderedDict()
    for uniq_feature in features_map:
        if uniq_feature == REPEAT_UNIT_LENGTH_4:
            var_subgroup = "Length_4"
        elif uniq_feature == REPEAT_UNIT_LENGTH_5:
            var_subgroup = "Length_5"
        else:
            var_subgroup = uniq_feature
        features_template[uniq_feature] = {}
        features_template[uniq_feature][VARIANT_TYPE] = "MSI"
        features_template[uniq_feature][VARIANT_SUBGROUP] = var_subgroup
        features_template[uniq_feature][FEATURE_QUANTITY] = 0
    return features_template, features_hash

MSI_FEATURES_TEMPLATE, MSI_FEATURES_HASH = __init_msi_features_template()
