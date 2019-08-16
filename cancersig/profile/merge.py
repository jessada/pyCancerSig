import copy
import sys
import os
import fnmatch
from collections import defaultdict
from os.path import join as join_path
from cancersig.template import pyCancerSigBase
from cancersig.profile.features import PROFILE_TYPES
from cancersig.profile.features import PROFILE_TYPE_SNV
from cancersig.profile.features import PROFILE_TYPE_SV
from cancersig.profile.features import PROFILE_TYPE_MSI
from cancersig.profile.features import PROFILE_WEIGHTS
from cancersig.profile.features import VARIANT_TYPE
from cancersig.profile.features import VARIANT_SUBGROUP
from cancersig.profile.features import FEATURE_ID
from cancersig.profile.features import SNV_FEATURES_TEMPLATE
from cancersig.profile.features import SV_FEATURES_TEMPLATE
from cancersig.profile.features import MSI_FEATURES_TEMPLATE

class ProfileMerger(pyCancerSigBase):

    def __init__(self, *args, **kwargs):
        super(ProfileMerger, self).__init__(*args, **kwargs)

    def __load_profile(self, profile_file):
        profile_type = PROFILE_TYPE_SNV
        sum_quantity = 0
        quantity_dict = {}
        with open(profile_file) as f_p:
            profile_id = f_p.readline().strip().split("\t")[3]
            for line in f_p:
                content = line.strip().split("\t")
                feature_id = content[2]
                quantity = float(content[3])
                quantity_dict[feature_id] = quantity
                sum_quantity += quantity
                if feature_id == "INV_log10_6_7":
                    profile_type = PROFILE_TYPE_SV
                if feature_id == "Repeat_unit_length_4":
                    profile_type = PROFILE_TYPE_MSI
        return quantity_dict, sum_quantity, profile_type, profile_id

    def __validate_merged_profile(self,
                                  sample_profiles,
                                  input_profile_types,
            
                                  ):
        samples_list = sample_profiles.keys()
        samples_to_be_merged = []
        for sample_id in samples_list:
            exclude = False
            for profile_type in input_profile_types:
                if profile_type not in sample_profiles[sample_id]:
                    self.warning("******************************   W A R N I N G   ******************************")
                    self.warning("   There is no profile type: " + profile_type + " for sample: " + sample_id)
                    self.warning("   The sample will be excluded from the merged file")
                    self.warning("   You can complete the missing profile type and rerun cancersig profile merge again")
                    exclude = True
                    break
            if exclude:
                continue
            samples_to_be_merged.append(sample_id)
        return samples_to_be_merged

    def __write_output_features(self,
                                f_o,
                                sample_profiles,
                                samples_list,
                                profile_type,
                                ):
        if profile_type == PROFILE_TYPE_SNV:
            features_dict = SNV_FEATURES_TEMPLATE
        if profile_type == PROFILE_TYPE_SV:
            features_dict = SV_FEATURES_TEMPLATE
        if profile_type == PROFILE_TYPE_MSI:
            features_dict = MSI_FEATURES_TEMPLATE
        for feature_id in features_dict:
            content = features_dict[feature_id][VARIANT_TYPE]
            content += "\t" + features_dict[feature_id][VARIANT_SUBGROUP]
            content += "\t" + feature_id
            for sample_id in samples_list:
                if profile_type in sample_profiles[sample_id]:
                    content += "\t{:.8f}".format(sample_profiles[sample_id][profile_type][feature_id])
                else:
                    content += "\t{:.8f}".format(0.00000001)
            f_o.write(content+"\n")
    
    def __merge(self,
                input_dirs,
                output_file,
                input_profile_types,
                ):
        total_weight = 0
        for profile_type in input_profile_types:
            total_weight += PROFILE_WEIGHTS[profile_type]
        sample_profiles = defaultdict(dict)
        for input_dir in input_dirs:
            self.info("Scanning: " + input_dir)
            self.info()
            for file_name in os.listdir(input_dir):
                if fnmatch.fnmatch(file_name, "*profile.txt"):
                    profile_file = join_path(input_dir, file_name)
                elif fnmatch.fnmatch(file_name, "*feature.txt"):
                    profile_file = join_path(input_dir, file_name)
                else:
                    continue
                self.info(">> Loading profile: " + profile_file)
                quantity_dict, sum_quantity, profile_type, profile_id = self.__load_profile(profile_file)
                weight = 0
                if profile_type == PROFILE_TYPE_SNV:
                    weight = PROFILE_WEIGHTS[PROFILE_TYPE_SNV]/total_weight
                    sample_profiles[profile_id][profile_type] = copy.deepcopy(SNV_FEATURES_TEMPLATE)
                    expected_feature = list(SNV_FEATURES_TEMPLATE.keys())
                if profile_type == PROFILE_TYPE_SV:
                    weight = PROFILE_WEIGHTS[PROFILE_TYPE_SV]/total_weight
                    sample_profiles[profile_id][profile_type] = copy.deepcopy(SV_FEATURES_TEMPLATE)
                    expected_feature = list(SV_FEATURES_TEMPLATE.keys())
                if profile_type == PROFILE_TYPE_MSI:
                    weight = PROFILE_WEIGHTS[PROFILE_TYPE_MSI]/total_weight
                    sample_profiles[profile_id][profile_type] = copy.deepcopy(MSI_FEATURES_TEMPLATE)
                    expected_feature = list(MSI_FEATURES_TEMPLATE.keys())
                for feature_id in quantity_dict:
                    if feature_id not in expected_feature:
                        raise ValueError("Unknown feature: " + feature_id + " in file: " + profile_file)
                    if sum_quantity < 0.0001:
                        sample_profiles[profile_id][profile_type][feature_id] = 0.00000001
                    else:
                        sample_profiles[profile_id][profile_type][feature_id] = ((quantity_dict[feature_id]*weight)/sum_quantity) + 0.00000001
            self.info()
        self.info("Validating and merging all input profiles")
        samples_to_be_merged = self.__validate_merged_profile(sample_profiles,
                                                              input_profile_types)
        samples_to_be_merged.sort()
        with open(output_file, "w") as f_o:
            header = VARIANT_TYPE
            header += "\t" + VARIANT_SUBGROUP
            header += "\t" + FEATURE_ID
            header += "\t" + "\t".join(samples_to_be_merged)
            f_o.write(header+"\n")
            for profile_type in PROFILE_TYPES:
                if profile_type in input_profile_types:
                    self.__write_output_features(f_o,
                                                 sample_profiles,
                                                 samples_to_be_merged,
                                                 profile_type)
        self.info("DONE!! Profiles have been merged and written to: " + output_file)

    def merge(self,
              input_dirs,
              output_file,
              profile_types=PROFILE_TYPES,
              ):
        self.__merge(input_dirs,
                     output_file,
                     profile_types,
                     )
