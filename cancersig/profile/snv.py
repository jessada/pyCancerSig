import copy
from os.path import join as join_path
from cancersig.config import CANCERSIG_SCRIPTS_DIR
from cancersig.utils import exec_sh
from cancersig.template import pyCancerSigBase
from cancersig.profile.features import SNV_FEATURES_TEMPLATE
from cancersig.profile.features import SNV_FEATURES_HASH
from cancersig.profile.features import VARIANT_TYPE
from cancersig.profile.features import VARIANT_SUBGROUP
from cancersig.profile.features import FEATURE_ID
from cancersig.profile.features import FEATURE_QUANTITY

SCRIPT_COUNT_SNV_EVENTS = join_path(CANCERSIG_SCRIPTS_DIR, "identify_snv_events.sh")

class SNVProfiler(pyCancerSigBase):

    def __init__(self, *args, **kwargs):
        super(SNVProfiler, self).__init__(*args, **kwargs)

    def __count_snv_events(self,
                           input_vcf_file,
                           ref_genome_file,
                           output_event_file,
                           gt_format="GTR",
                           sample_id=None,
                           ):
        self.info("counting SNV events and saving them to " + output_event_file)
        cmd = SCRIPT_COUNT_SNV_EVENTS
        cmd +=" -i " + input_vcf_file
        cmd +=" -r " + ref_genome_file
        cmd +=" -o " + output_event_file
        cmd +=" -g " + gt_format
        p, stdout_data = exec_sh(cmd, silent=True)

    def __event_to_profile(self,
                           input_event_file,
                           output_file,
                           sample_id,
                           ):
        features_count = copy.deepcopy(SNV_FEATURES_TEMPLATE)
        with open(input_event_file) as f_i:
            for line in f_i:
                if line[0] == "#":
                    if sample_id is None:
                        sample_id = line.strip().split("\t")[7]
                    continue
                dummy, dummy, dummy, dummy, ref, alt, triplet, dummy = line.strip().split("\t")
                feature_id = SNV_FEATURES_HASH[ref][alt][triplet]
                features_count[feature_id][FEATURE_QUANTITY] += 1
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
                ref_genome_file,
                output_file,
                gt_format="GTR",
                sample_id=None,
                ):
        event_file = output_file + ".event"
        self.__count_snv_events(input_vcf_file,
                                  ref_genome_file,
                                  event_file,
                                  gt_format,
                                  sample_id,
                                  )
        self.__event_to_profile(event_file, output_file, sample_id)

