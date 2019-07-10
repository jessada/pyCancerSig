import copy
import numpy
from cancersig.template import pyCancerSigBase
from cancersig.utils import readVCF
from cancersig.profile.features import SV_FEATURES_TEMPLATE
from cancersig.profile.features import SV_FEATURES_HASH
from cancersig.profile.features import VARIANT_TYPE
from cancersig.profile.features import VARIANT_SUBGROUP
from cancersig.profile.features import FEATURE_ID
from cancersig.profile.features import FEATURE_QUANTITY
from cancersig.profile.features import SV_LEN_LOG10_2_3
from cancersig.profile.features import SV_LEN_LOG10_3_4
from cancersig.profile.features import SV_LEN_LOG10_4_5
from cancersig.profile.features import SV_LEN_LOG10_5_6
from cancersig.profile.features import SV_LEN_LOG10_6_7
from cancersig.profile.features import SV_LEN_LOG10_7_8
from cancersig.profile.features import SV_LEN_LOG10_8_9
from cancersig.profile.features import SV_LEN_LOG10_9up
from cancersig.profile.features import SMALL_QUANTITY

class SVProfiler(pyCancerSigBase):

    def __init__(self, *args, **kwargs):
        super(SVProfiler, self).__init__(*args, **kwargs)

    def __profile(self,
                  input_vcf_file,
                  output_file,
                  sample_id=None,
                  ):
        features_count = copy.deepcopy(SV_FEATURES_TEMPLATE)
        with open(input_vcf_file) as f_i:
            for line in f_i:
                if line[0:2] == "##":
                    continue
                if (line[0:2] == "#C"):
                    variation = line.strip().split("\t")
                    if (len(variation) > 8) and (sample_id is None):
                        sample_id = variation[9]
                    continue
                chrA, posA, chrB, posB,event_type,INFO,format = readVCF.readVCFLine(line)
        
                if chrB not in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]:
                    continue
        
                # set default signal strength to be very high to keep result sfrom CNVnator
                signalPE=1000000
        
                if "PE" in format:
                    signalPE =int(format["PE"][0])
                elif  "DV" in format:
                    signalPE = int(format["DV"][0])
        
                if signalPE < 10:
                    continue
        
                length = float("inf")
                if chrA == chrB:
                    length = int(posB)-int(posA)
        
                len_log10 = numpy.log10(length)
                if length == float("inf"):
                    feature_id = SV_FEATURES_HASH[event_type][SV_LEN_LOG10_9up]
                elif len_log10 > 9:
                    feature_id = SV_FEATURES_HASH[event_type][SV_LEN_LOG10_9up]
                elif len_log10 > 8:
                    feature_id = SV_FEATURES_HASH[event_type][SV_LEN_LOG10_8_9]
                elif len_log10 > 7:
                    feature_id = SV_FEATURES_HASH[event_type][SV_LEN_LOG10_7_8]
                elif len_log10 > 6:
                    feature_id = SV_FEATURES_HASH[event_type][SV_LEN_LOG10_6_7]
                elif len_log10 > 5:
                    feature_id = SV_FEATURES_HASH[event_type][SV_LEN_LOG10_5_6]
                elif len_log10 > 4:
                    feature_id = SV_FEATURES_HASH[event_type][SV_LEN_LOG10_4_5]
                elif len_log10 > 3:
                    feature_id = SV_FEATURES_HASH[event_type][SV_LEN_LOG10_3_4]
                elif len_log10 > 2:
                    feature_id = SV_FEATURES_HASH[event_type][SV_LEN_LOG10_2_3]
                else:
                    self.info("variant at " + str(chrA) + ":" + str(posA) + " is too small to be considered")
                    # too small to be considered as any features
                    continue
                features_count[feature_id][FEATURE_QUANTITY] += 1
       
        # count total event
        total_event = 0
        for feature_id in features_count:
            total_event += features_count[feature_id][FEATURE_QUANTITY]
        
        with open(output_file, "w") as f_o:
            header = VARIANT_TYPE
            header += "\t" + VARIANT_SUBGROUP
            header += "\t" + FEATURE_ID
            header += "\t" + sample_id
            f_o.write(header+"\n")
            for feature_id in features_count:
                f_o.write("{:s}\t{:s}\t{:s}\t{:d}\n".format(features_count[feature_id][VARIANT_TYPE],
                                                            features_count[feature_id][VARIANT_SUBGROUP],
                                                            feature_id,
                                                            features_count[feature_id][FEATURE_QUANTITY],
                                                            ))

    def profile(self,
                input_vcf_file,
                output_file,
                sample_id=None,
                ):
        self.__profile(input_vcf_file,
                       output_file,
                       sample_id,
                       )
