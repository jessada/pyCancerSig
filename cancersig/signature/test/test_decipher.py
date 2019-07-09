#import filecmp
from os.path import join as join_path
from cancersig.template import Tester
from cancersig.signature.decipher import CancerSigNMFController

class TestCancerSigController(Tester):

    def __init__(self, methodName):
        super(TestCancerSigController, self).__init__(methodName=methodName,
                                                      test_module_name=__name__,
                                                      )

    def setUp(self):
        self.__cancersig_nmf_controller = CancerSigNMFController(debug_mode=True, verbose=True)

    def test_decipher_1(self):
        """ test decipher complete profile (all profile type) with 2 signatures """

        self.init_test(self.current_func_name)
        mutation_profiles = join_path(self.data_dir,
                                      "mutation_profiles.txt")
        output_prefix = join_path(self.working_dir,
                                  self.current_func_name)
        self.__cancersig_nmf_controller.decipher(mutation_profiles,
                                                 output_prefix,
                                                 min_signatures=2,
                                                 max_signatures=2,
                                                 )

#        output_feature_file = join_path(self.working_dir,
#                                        self.current_func_name+".txt")
#        output_event_file = output_feature_file + ".event"
#        self.__snv_profiler.profile(input_vcf_file,
#                                    ref_genome_file,
#                                    output_feature_file,
#                                    sample_id="example_sample_id_from_argument",
#                                    )
#        exp_output_feature_file = join_path(self.data_dir,
#                                             "exp_output_feature_file")
#        exp_output_event_file = join_path(self.data_dir,
#                                          "exp_output_event_file")
#        self.assertTrue(filecmp.cmp(output_event_file, exp_output_event_file), 'Malfunction in SNVProfiler.profile()')
#        self.assertTrue(filecmp.cmp(output_feature_file, exp_output_feature_file), 'Malfunction in SNVProfiler.profile()')

#    def test_decipher_2(self):
#        """ test decipher partial profile (SV only) with 2 signatures """
#
#        self.init_test(self.current_func_name)
#        mutation_profiles = join_path(self.data_dir,
#                                      "mutation_profiles.txt")
#        output_prefix = join_path(self.working_dir,
#                                  self.current_func_name)
#        self.__cancersig_nmf_controller.decipher(mutation_profiles,
#                                                 output_prefix,
#                                                 min_signatures=2,
#                                                 max_signatures=2,
#                                                 )
