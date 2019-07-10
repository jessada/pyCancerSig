import unittest
import filecmp
from os.path import join as join_path
from cancersig.config import ENABLE_CPU_INTENSIVE_UNITTEST
from cancersig.template import Tester
from cancersig.signature.visualize import Visualizer

class TestVisualizer(Tester):

    def __init__(self, methodName):
        super(TestVisualizer, self).__init__(methodName=methodName,
                                             test_module_name=__name__,
                                             )

    def setUp(self):
        self.__vz = Visualizer(debug_mode=True, verbose=True)

    @unittest.skipUnless(ENABLE_CPU_INTENSIVE_UNITTEST, "This test was disabled due to its computational burden, you can enaable it in cancersig.config")
    def test_visualize_1(self):
        """ test standard """

        self.init_test(self.current_func_name)
        mutation_profiles = join_path(self.data_dir,
                                      "mutation_profiles.txt")
        signatures_probabilities = join_path(self.data_dir,
                                      "signatures_probabilities.txt")
        output_dir = self.working_dir
        self.__vz.visualize(mutation_profiles,
                            signatures_probabilities,
                            output_dir,
                            )
        output_weights_file = self.__vz.weights_file
        exp_weights = join_path(self.data_dir,
                                "exp_weights")
        self.assertTrue(filecmp.cmp(output_weights_file, exp_weights), 'Malfunction in Visualizer.visualize()')
