import filecmp
from os.path import join as join_path
from cancersig.template import Tester
from cancersig.profile.merge import ProfileMerger

class TestProfileMerger(Tester):

    def __init__(self, methodName):
        super(TestProfileMerger, self).__init__(methodName=methodName,
                                                test_module_name=__name__,
                                                )

    def setUp(self):
        self.__profile_merger = ProfileMerger(verbose=False)
        self.__profile_merger.debug_mode = True

    def test_merge_1(self):
        """ test complete profile merge (default) """

        self.init_test(self.current_func_name)

        input_dirs = []
        input_dirs.append(join_path(join_path(self.data_dir,
                                              "SNV"),
                                    "cancer1"))
        input_dirs.append(join_path(join_path(self.data_dir,
                                              "SNV"),
                                    "cancer2"))
        input_dirs.append(join_path(join_path(self.data_dir,
                                              "SV"),
                                    "cancer1"))
        input_dirs.append(join_path(join_path(self.data_dir,
                                              "SV"),
                                    "cancer2"))
        input_dirs.append(join_path(join_path(self.data_dir,
                                              "MSI"),
                                    "cancer1"))
        input_dirs.append(join_path(join_path(self.data_dir,
                                              "MSI"),
                                    "cancer2"))
        output_file = join_path(self.working_dir,
                                self.current_func_name+".txt")
        self.__profile_merger.merge(input_dirs, output_file)
        exp_output_file = join_path(self.data_dir,
                                    "exp_output_file")
        self.assertTrue(filecmp.cmp(output_file, exp_output_file), 'Malfunction in ProfileMerger.merge()')

    def test_merge_2(self):
        """ test merge only one profile type """

        self.init_test(self.current_func_name)

        input_dirs = []
        input_dirs.append(join_path(join_path(self.data_dir,
                                              "SNV"),
                                    "cancer1"))
        input_dirs.append(join_path(join_path(self.data_dir,
                                              "SNV"),
                                    "cancer2"))
        input_dirs.append(join_path(join_path(self.data_dir,
                                              "SV"),
                                    "cancer1"))
        input_dirs.append(join_path(join_path(self.data_dir,
                                              "SV"),
                                    "cancer2"))
        output_file = join_path(self.working_dir,
                                self.current_func_name+".txt")
        self.__profile_merger.merge(input_dirs, output_file, ["SV"])
        exp_output_file = join_path(self.data_dir,
                                    "exp_output_file")
        self.assertTrue(filecmp.cmp(output_file, exp_output_file), 'Malfunction in ProfileMerger.merge()')
