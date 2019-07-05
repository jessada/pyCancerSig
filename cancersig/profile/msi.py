import copy
from cancersig.utils import exec_sh
from cancersig.template import pyCancerSigBase
from cancersig.profile.features import MSI_FEATURES_TEMPLATE
from cancersig.profile.features import MSI_FEATURES_HASH
from cancersig.profile.features import VARIANT_TYPE
from cancersig.profile.features import VARIANT_SUBGROUP
from cancersig.profile.features import FEATURE_ID
from cancersig.profile.features import FEATURE_QUANTITY

MSI_H_CUTOFF = 3.5

class MSIProfiler(pyCancerSigBase):

    def __init__(self, *args, **kwargs):
        super(MSIProfiler, self).__init__(*args, **kwargs)

    def __profile(self,
                  raw_msisensor_out_somatic,
                  sample_id,
                  output_file,
                  msi_positive,
                  ):
        features_count = copy.deepcopy(MSI_FEATURES_TEMPLATE)
        if msi_positive:
            cmd = "cut -f 5 " + raw_msisensor_out_somatic
            cmd += " | sort"
            cmd += " | uniq -c"
            p, stdout_data = exec_sh(cmd, silent=True)
            unstable_loci_counts = list(map(lambda x: x.strip(), stdout_data.decode('utf-8').strip().split("\n")))
            total_count = 0
            for unstable_loci_count in unstable_loci_counts:
                _count, _repeat = unstable_loci_count.split()
                _count = int(_count)
                total_count += _count
                if len(_repeat) == 4:
                    features_count[REPEAT_UNIT_LENGTH_4][FEATURE_QUANTITY] += _count
                    continue
                if len(_repeat) == 5:
                    features_count[REPEAT_UNIT_LENGTH_5][FEATURE_QUANTITY] += _count
                    continue
                feature = MSI_FEATURES_HASH[_repeat]
                features_count[feature][FEATURE_QUANTITY] += _count
        with open(output_file, "w") as f_o:
            header = VARIANT_TYPE
            header += "\t" + VARIANT_SUBGROUP
            header += "\t" + FEATURE_ID
            header += "\t" + sample_id
            f_o.write(header+"\n")
            for feature_id in features_count:
                if msi_positive:
                    feature_quantity = features_count[feature_id][FEATURE_QUANTITY]
                else:
                    feature_quantity = 0
                f_o.write("{:s}\t{:s}\t{:s}\t{:d}\n".format(features_count[feature_id][VARIANT_TYPE],
                                                            features_count[feature_id][VARIANT_SUBGROUP],
                                                            feature_id,
                                                            feature_quantity,
                                                            ))

    def profile(self,
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
            self.__profile(raw_msisensor_out_somatic,
                           sample_id,
                           output_file,
                           msi_positive=msi_positive,
                           )
