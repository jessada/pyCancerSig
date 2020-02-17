import copy
import sys
import re
from os.path import join as join_path
from tempfile import mkdtemp
from pyfaidx import Faidx
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

    def profile(self,
                input_vcf_file,
                ref_genome_file,
                output_file,
                raw_gt_format="GTR",
                sample_id=None,
                ):
        # unzip decompose and  clean vcf file
        clean_vcf_file = join_path(mkdtemp(), "clean.vcf")
        if input_vcf_file.endswith(".gz"):
            cmd = "gunzip -c " + input_vcf_file
        else:
            cmd = "cat " + input_vcf_file
        cmd += " | vt decompose -s -"
        cmd += " | grep -Pv \"\t\*\t\""
        cmd += " | grep -v \"\\x3b\""
        cmd += " | grep -v \"^M\""
        cmd += " > " + clean_vcf_file
        p, stdout_data = exec_sh(cmd, silent=True)

        # parse the clean vcf file for the required fields
        vcf_query_format = "'"
        vcf_query_format += "%CHROM"
        vcf_query_format += "\t%POS"
        vcf_query_format += "\t%REF"
        vcf_query_format += "\t%ALT"
        vcf_query_format += "[\t%SAMPLE=%" + raw_gt_format + "]"
        vcf_query_format += "\n"
        vcf_query_format += "'"
        cmd = "vcf-query"
        cmd += " -f " + vcf_query_format
        cmd += " " + clean_vcf_file
        p, stdout_data = exec_sh(cmd, silent=True)

        # get list of smaple id and prepare data structure
        first_variant_record = stdout_data.decode('utf-8').split("\n")[0]
        variant_items = first_variant_record.strip().split("\t")
        samples_features = {}
        for sample_idx in range(4, len(variant_items)):
            gt_data = variant_items[sample_idx] 
            m = re.match(r"(?P<sample_id>.*)=(?P<raw_gt>.*)", gt_data)
            sample_id = m.group("sample_id")
            samples_features[sample_id] = copy.deepcopy(SNV_FEATURES_TEMPLATE)

        # iterate over all vcf record and count variants for each sample
        fa = Faidx(ref_genome_file)
        for variant_record in stdout_data.decode('utf-8').split("\n"):
            variant_items = variant_record.strip().split("\t")
            if len(variant_items) < 4:
                continue
            chrom = variant_items[0]
            pos = variant_items[1]
            ref = variant_items[2]
            alt = variant_items[3]
            if len(ref) > 1:
                continue
            if ref == "-":
                continue
            if len(alt) > 1:
                continue
            if alt == "-":
                continue
            triplet = fa.fetch(chrom, int(pos)-1, int(pos)+1).seq
            feature_id = SNV_FEATURES_HASH[ref][alt][triplet]
            # iterate over all samples in the record
            for sample_idx in range(4, len(variant_items)):
                gt_data = variant_items[sample_idx] 
                m = re.match(r"(?P<sample_id>.*)=(?P<raw_gt>.*)", gt_data)
                raw_gt = m.group("raw_gt")
                if raw_gt == "0/0":
                    continue
                samples_features[m.group("sample_id")][feature_id][FEATURE_QUANTITY] += 1
        fa.close()

        # write output feature file
        with open(output_file, "w") as f_o:
            header = VARIANT_TYPE
            header += "\t" + VARIANT_SUBGROUP
            header += "\t" + FEATURE_ID
            for sample_id in samples_features:
                header += "\t" + sample_id
            f_o.write(header+"\n")
            for feature_id in SNV_FEATURES_TEMPLATE:
                feature_info =  "{:s}\t{:s}\t{:s}".format(SNV_FEATURES_TEMPLATE[feature_id][VARIANT_TYPE],
                                                          SNV_FEATURES_TEMPLATE[feature_id][VARIANT_SUBGROUP],
                                                          feature_id,
                                                          )
                for sample_id in samples_features:
                    feature_info += "\t" + str(samples_features[sample_id][feature_id][FEATURE_QUANTITY])
                f_o.write(feature_info + "\n")

        self.info()
        self.info("Done!! The output file is at " + output_file)
