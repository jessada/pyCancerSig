from collections import OrderedDict

VARIANT_TYPE = "variant type"
VARIANT_SUBGROUP = "variant subgroup"
FEATURE_ID = "feature id"
FEATURE_QUANTITY = "feature quantity"

VAR_TYPE_NA = "NA"
SMALL_QUANTITY = 0.00000001

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
            features_template[feature_id][FEATURE_QUANTITY] = 0.0
            features_hash[event_type][len_log10] = feature_id
    return features_template, features_hash
        
SV_FEATURES_TEMPLATE, SV_FEATURES_HASH = __init_sv_features_template()
