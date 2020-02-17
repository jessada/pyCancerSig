import unittest
import filecmp
from os.path import join as join_path
from cancersig.config import ENABLE_CPU_INTENSIVE_UNITTEST
from cancersig.template import Tester
from cancersig.profile.snv import SNVProfiler

class TestSNVProfiler(Tester):

    def __init__(self, methodName):
        super(TestSNVProfiler, self).__init__(methodName=methodName,
                                             test_module_name=__name__,
                                             )

    def setUp(self):
        self.__snv_profiler = SNVProfiler()
        self.__snv_profiler.debug_mode = True

    def test_profile_1(self):
        """ test extracting snv features on standard argument (sample id from vcf header) """

        self.init_test(self.current_func_name)
        input_vcf_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        ref_genome_file = join_path(self.data_dir,
                                    "ref_genome.fa")
        output_feature_file = join_path(self.working_dir,
                                        self.current_func_name+".txt")
        output_event_file = output_feature_file + ".event"
        self.__snv_profiler.profile(input_vcf_file,
                                    ref_genome_file,
                                    output_feature_file,
                                    )
        exp_output_feature_file = join_path(self.data_dir,
                                            "exp_snv_feature.txt")
        self.assertTrue(filecmp.cmp(output_feature_file, exp_output_feature_file), 'Malfunction in SNVProfiler.profile()')

    def test_profile_2(self):
        """ test simple vcf, a mockup one with as litle information as possible """

        self.init_test(self.current_func_name)
        input_vcf_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        ref_genome_file = join_path(self.data_dir,
                                    "ref_genome.fa")
        output_feature_file = join_path(self.working_dir,
                                        self.current_func_name+".txt")
        output_event_file = output_feature_file + ".event"
        self.__snv_profiler.profile(input_vcf_file,
                                    ref_genome_file,
                                    output_feature_file,
                                    )
        exp_output_feature_file = join_path(self.data_dir,
                                            "exp_snv_feature.txt")
        self.assertTrue(filecmp.cmp(output_feature_file, exp_output_feature_file), 'Malfunction in SNVProfiler.profile()')

    def test_profile_3(self):
        """ test a vcf with other non-standard chromosomes """

        self.init_test(self.current_func_name)
        input_vcf_file = join_path(self.data_dir,
                                   "input.vcf.gz")
        ref_genome_file = join_path(self.data_dir,
                                    "ref_genome.fa")
        output_feature_file = join_path(self.working_dir,
                                        self.current_func_name+".txt")
        output_event_file = output_feature_file + ".event"
        self.__snv_profiler.profile(input_vcf_file,
                                    ref_genome_file,
                                    output_feature_file,
                                    )
        exp_output_feature_file = join_path(self.data_dir,
                                            "exp_snv_feature.txt")
        self.assertTrue(filecmp.cmp(output_feature_file, exp_output_feature_file), 'Malfunction in SNVProfiler.profile()')

    def test_profile_4(self):
        """ test a vcf without gz suffix """

        self.init_test(self.current_func_name)
        input_vcf_file = join_path(self.data_dir,
                                   "input.vcf")
        ref_genome_file = join_path(self.data_dir,
                                    "ref_genome.fa")
        output_feature_file = join_path(self.working_dir,
                                        self.current_func_name+".txt")
        output_event_file = output_feature_file + ".event"
        self.__snv_profiler.profile(input_vcf_file,
                                    ref_genome_file,
                                    output_feature_file,
                                    )
        exp_output_feature_file = join_path(self.data_dir,
                                            "exp_snv_feature.txt")
        self.assertTrue(filecmp.cmp(output_feature_file, exp_output_feature_file), 'Malfunction in SNVProfiler.profile()')

    def test_profile_5(self):
        """ test a vcf with multiple samples """

        self.individual_debug = True
        self.init_test(self.current_func_name)
        input_vcf_file = join_path(self.data_dir,
                                   "input.vcf")
        ref_genome_file = join_path(self.data_dir,
                                    "ref_genome.fa")
        output_feature_file = join_path(self.working_dir,
                                        self.current_func_name+".txt")
        output_event_file = output_feature_file + ".event"
        self.__snv_profiler.profile(input_vcf_file,
                                    ref_genome_file,
                                    output_feature_file,
                                    )
        exp_output_feature_file = join_path(self.data_dir,
                                            "exp_snv_feature.txt")
        self.assertTrue(filecmp.cmp(output_feature_file, exp_output_feature_file), 'Malfunction in SNVProfiler.profile()')
