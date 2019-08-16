import unittest
from os.path import join as join_path
from cancersig.config import ENABLE_CPU_INTENSIVE_UNITTEST
from cancersig.template import Tester
from cancersig.signature.decipher import CancerSigNMFController

class TestCancerSigController(Tester):

    def __init__(self, methodName):
        super(TestCancerSigController, self).__init__(methodName=methodName,
                                                      test_module_name=__name__,
                                                      )

    def setUp(self):
        self.__cancersig_nmf_controller = CancerSigNMFController(debug_mode=True, verbose=False)

    @unittest.skipUnless(ENABLE_CPU_INTENSIVE_UNITTEST, "This test was disabled due to its computational burden, you can enaable it in cancersig.config")
    def test_decipher_1(self):
        """ test decipher complete profile (all profile type) with 2 signatures and standard model parameters"""

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

    @unittest.skipUnless(ENABLE_CPU_INTENSIVE_UNITTEST, "This test was disabled due to its computational burden, you can enaable it in cancersig.config")
    def test_decipher_2(self):
        """ test decipher partial profile (SV only) with 2 signatures """

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

    def test_decipher_3(self):
        """ test with parameters for just fast going through the model """

        self.init_test(self.current_func_name)
        mutation_profiles = join_path(self.data_dir,
                                      "mutation_profiles.txt")
        output_prefix = join_path(self.working_dir,
                                  self.current_func_name)
        cancersig_nmf_controller = CancerSigNMFController(debug_mode=True,
                                                          verbose=True,
                                                          resampling_size=100,
                                                          convergence_cutoff=0.01,
                                                          weak_signal_cutoff=0.00,
                                                          max_model_solutions=5,
                                                          )
        cancersig_nmf_controller.decipher(mutation_profiles,
                                          output_prefix,
                                          min_signatures=2,
                                          max_signatures=2,
                                          save_signature_pdfs=False,
                                          )
