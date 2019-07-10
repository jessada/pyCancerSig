from cancersig.template import pyCancerSigBase
import pandas
import numpy as np
import re
import math
import os
import matplotlib.pyplot as plt
from os.path import join as join_path
from os.path import basename
from os.path import splitext
from scipy.optimize import minimize_scalar
from collections import defaultdict
from collections import OrderedDict
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors as mcolors
from cancersig.signature.figure import plot_distribution
from cancersig.signature.figure import get_profile_dict
from cancersig.profile.features import VARIANT_TYPE
from cancersig.profile.features import VARIANT_SUBGROUP
from cancersig.profile.features import FEATURE_ID

WEIGHT_CUTOFF = 0.06

# define cancer processes file header
SIGNATURES_PROB_HEADER_REGEX_SYNTAX = "(" + VARIANT_TYPE + ")"
SIGNATURES_PROB_HEADER_REGEX_SYNTAX += "|(" + VARIANT_SUBGROUP + ")"
SIGNATURES_PROB_HEADER_REGEX_SYNTAX += "|(" + FEATURE_ID + ")"
SIGNATURES_PROB_HEADER_REGEX_SYNTAX += "|(Unnamed)"

class Visualizer(pyCancerSigBase):
    def __init__(self, *args, **kwargs):
        super(Visualizer, self).__init__(*args, **kwargs)


    @property
    def output_root_dir(self):
        return self.__output_root_dir

    @property
    def figures_dir(self):
        return join_path(self.output_root_dir, "figures")

    @property
    def weights_dir(self):
        if self.__weights_dir is None:
            self.__weights_dir = join_path(self.output_root_dir, "weights")
            if not os.path.exists(self.__weights_dir):
                os.makedirs(self.__weights_dir)
        return self.__weights_dir

    @property
    def weights_file(self):
        return join_path(self.weights_dir, "normalized_weights.txt")

    def __load_sample_profiles(self, file_path):
        """
        Load sample profiles, which are in the same format as the deciphered
        cancer processes. The first three columns are fix, "variant type",
        "variant subgroup", and "feature id". The rest are for the samples to be
        reconstructed. Each column represents one sample.
        """
        df = pandas.read_csv(file_path, sep='\t', engine='python')
        samples_id = df.columns.values[3:]
        samples_context_ratio = np.array(df.loc[:, filter(lambda x: not re.search(SIGNATURES_PROB_HEADER_REGEX_SYNTAX, x), df.columns.values)])
        samples_context_ratio = samples_context_ratio/sum(samples_context_ratio)
        # parse substition for each context for each sample
        n_rows = df.shape[0]
        samples_vars_dict = {}
        for sample_id in samples_id:
            sample_vars_dict = defaultdict(dict)
            for row_idx in range(n_rows):
                row_df = df.iloc[row_idx]
                variant_type = row_df[VARIANT_TYPE]
                variant_info = row_df[VARIANT_SUBGROUP]
                context_count = row_df[sample_id]
                sample_vars_dict[variant_type][variant_info] = float(context_count)
            samples_vars_dict[sample_id] = sample_vars_dict
        return samples_id, samples_context_ratio, samples_vars_dict

    def __load_signatures_probabilities(self, file_path):
        """ Load siganture probabilities"""
        df = pandas.read_csv('{}'.format(file_path), sep='\t', engine='python')
        self.__signature_names = df.columns.values[3:]
        self.__mutation_prob = np.array(df.loc[:, filter(lambda x: not re.search(SIGNATURES_PROB_HEADER_REGEX_SYNTAX, x), df.columns.values)])
        self.__features = df[FEATURE_ID]
        self.__feature_infos = df.iloc[:, 0:3]

    def __calculate_ignorable_signatures(self, sample_context_ratio):
        """
        Calculates which signatures can be ignored because they contain a peak for a context that is clearly not
        seen in the tumor data.
        :return: List of dicts of signatures to ignore along with contextual information for why.
        """
        signatures_to_ignore = []
        for i, signature_name in enumerate(self.__signature_names):
            context_fractions = self.__mutation_prob[:, i]
            for j, cf in enumerate(context_fractions):
                if cf > 0.2:
                    somatic_mutation_type = self.__features[j]
                    if sample_context_ratio[j] < 0.01:
                        signatures_to_ignore.append({'name': signature_name,
                                                     'index': i,
                                                     'outlier_context': somatic_mutation_type,
                                                     'context_fraction': cf})
                    break
        return signatures_to_ignore

    def __seed_weights(self, sample_context_ratio, ignorable_indices=None):
        """
        Find which of the signatures best approximates the tumor signature, and seed the weights such that that
        signature is assigned weight 1 and all other signatures are assigned weight zero. These are the seed weights
        upon which the algorithm will build as it tries to further reduce sum of squared error.
        :param tumor: normalized array of shape (1, 96) where each entry is a mutation context fraction
        :param signatures: array of shape (96, num_signatures) where each row represents a mutation context and each
        column is a signature
        :return: normalized array of shape (num_signatures, 1) representing weight of each signature
        """
        if ignorable_indices is None:
            ignorable_indices = []

        num_sigs = len(self.__mutation_prob[0])
        ss_errors = np.empty(num_sigs, )
        ss_errors.fill(math.inf)
        for i in range(num_sigs):
            if i not in ignorable_indices:
                tmp_weights = np.zeros((num_sigs,))
                tmp_weights[i] = 1
                error = self.__get_error(sample_context_ratio, tmp_weights)
                ss_errors[i] = error
        # Seed index that minimizes sum of squared error metric
        seed_index = np.argmin(ss_errors, axis=0)
        final_weights = np.zeros(num_sigs)
        final_weights[seed_index] = 1
        return final_weights

    def __log_normalized_weights(self, weights):
        """A standard way to print normalized weights given a vector of potentially not yet normalized weights"""
        normalized_weights = weights / sum(weights)
        for i, weight in enumerate(normalized_weights):
            if weight != 0:
                self.info("\t\tSignature {}: {}".format(str(i+1), weight))

    def __update_weights(self, sample_context_ratio, weights, ignorable_signature_indices=None):
        """
        Given a set of initial weights, update the weights array with new values that shrink the sum of squares
        error metric.
        :param tumor: normalized array of shape (1, 96) where each entry is a mutation context fraction for the tumor
        :param signatures: signatures: array of shape (96, num_signatures) where each row represents a mutation context
        and each column is a signature
        :param w: array of shape (num_signatures, 1) representing weight of each signature
        :param signatures_limit: How many of the total signatures to consider when assigning weights
        :param ignorable_signature_indices: an array of indices into the signatures array indicating which to ignore
        :return: a new weights array, w.
        """
        if ignorable_signature_indices is None:
            ignorable_signature_indices = []

        # The number of signatures already being used in the current linear combination of signatures
        num_sigs_present = len([weight for weight in weights if weight != 0])

        # The total number of signatures to choose from
        num_sigs = np.shape(self.__mutation_prob)[1]

        # The current sum of squares error given the present weights assigned for each signature
        error_old = self.__get_error(sample_context_ratio, weights)

        # Which weight indices to allow changes for; if we haven't reached the limit all weights are fair game
        changeable_indices = range(num_sigs)
        changeable_indices = [i for i in changeable_indices if i not in ignorable_signature_indices]

        # zero square matrix of num signatures dimensions
        v = np.zeros((num_sigs, num_sigs))

        # 1 * num signatures vector with values preset to infinity
        new_squared_errors = np.empty(num_sigs, )
        new_squared_errors.fill(math.inf)

        # Only consider adjusting the weights which are allowed to change
        for i in changeable_indices:
            # Find the weight x for the ith signature that minimizes the sum of squared error
            def to_minimize(x):
                # Initialize a temporary zero vector of length number of signatures
                tmp = np.zeros((1, num_sigs))
                tmp[0, i] = x
                return self.__get_error(sample_context_ratio,
                                   weights + tmp[0,])

            error_minimizer = minimize_scalar(to_minimize, bounds=(-weights[i], 1),
                                              method="bounded").x
            v[i, i] = error_minimizer
            w_new = weights + v[i]
            new_squared_errors[i] = self.__get_error(sample_context_ratio,
                                                w_new)

        # Find which signature can be added to the weights vector to best reduce the error
        min_new_squared_error = min(new_squared_errors)

        # Update that signature within the weights vector with the new value that best reduces the overall error
        if min_new_squared_error < error_old:
            index_of_min = np.argmin(new_squared_errors, axis=0)
            weights[index_of_min] = weights[index_of_min] + v[index_of_min, index_of_min]

        return weights

    def __calculate_weights(self, sample_context_ratio):
        ignorable_signatures = self.__calculate_ignorable_signatures(sample_context_ratio)
        self.info("")
        self.info('Signatures ignored because of outlying contexts:')
        if len(ignorable_signatures) == 0:
            self.info("    None")
        else:
            for s in ignorable_signatures:
                self.info('>>>  {} ignored because of outlying context {} with fraction {}'
                                   .format(s.get('name'),
                                           s.get('outlier_context'),
                                           s.get('context_fraction')))
    
        ignorable_indices = [ig['index'] for ig in ignorable_signatures]
        weights = self.__seed_weights(sample_context_ratio,
                                      ignorable_indices=ignorable_indices)

        self.info("")
        self.info("Initial seed weights assigned:")
        self.__log_normalized_weights(weights)

        # perform actual weights calculation
        self.info("")
        iteration = 0
        error_diff = math.inf
        error_threshold = 1e-3
        self.info("Calculating weights")
        while error_diff > error_threshold:
            iteration = iteration + 1
            error_pre = self.__get_error(sample_context_ratio,
                                         weights)
            if error_pre == 0:
                break
            self.info("Iter {}:".format(iteration))
            self.info(">>>  Pre error: {}".format(error_pre))
            weights = self.__update_weights(sample_context_ratio,
                                            weights,
                                            ignorable_signature_indices=ignorable_indices)
            error_post = self.__get_error(sample_context_ratio, weights)
            self.info(">>>  Post error: {}".format(error_post))
            error_diff = (error_pre - error_post) / error_pre
            self.info(">>>  New normalized weights: ")
            self.__log_normalized_weights(weights)

        normalized_weights = weights/sum(weights)

        # Filter out any weights less than cut off
        np.place(normalized_weights, normalized_weights < WEIGHT_CUTOFF, 0)
        return normalized_weights

    def __get_reconstructed_tumor_profile(self, weights):
        """Reconstruct a tumor profile given a set of signatures and a vector of signature weights"""
        w_norm = weights/sum(weights)
        return w_norm.dot(np.transpose(self.__mutation_prob))

    def __get_error(self, sample_context_ratio, weights):
        """
        Calculate the SSE between the true tumor signature and the calculated linear combination of different signatures
        :param tumor: normalized array of shape (1, 96) where each entry is a mutation context fraction for the tumor
        :param signatures: array of shape (96, num_signatures) where each row represents a mutation context and each
        column is a signature
        :param w: array of shape (num_signatures, 1) representing weight of each signature
        :return: sum of squares error between reconstructed tumor context fractions and actual tumor profile
        """
        sample_context_ratio = sample_context_ratio/sum(sample_context_ratio)
        reconstructed_tumor_profile = self.__get_reconstructed_tumor_profile(weights)
        error = sample_context_ratio - reconstructed_tumor_profile
        squared_error_sum = np.sum(error.dot(np.transpose(error)))
        return squared_error_sum

    def __plot_pie_chart(self, weights, sample_id):
        all_colors = {self.__signature_names[i]: list(mcolors.CSS4_COLORS.values())[30:][i] for i in range(len(self.__signature_names))}
        # Below code looks strange because I want to keep the second argument for the explanation (similar to COSMIC)
        signature_weights_text = zip(self.__signature_names,
                                     self.__signature_names,
                                     weights)
        # Data to plot
        non_zero_weights = []
        non_zero_labels = []
        colors = []
        for process, process_explanation, weight in sorted(signature_weights_text, key=lambda t: t[2], reverse=True):
            if weight != 0:
                label = '{}, {}%'.format(process, round(weight*100, 2))
                color = all_colors[process]
                non_zero_labels.append(label)
                non_zero_weights.append(weight)
                colors.append(color)

        # Fill in the missing piece of the pie (due to removed below-threshold signatures)
        if 1 - sum(non_zero_weights) > 1e-3:
            difference = 1 - sum(non_zero_weights)
            percentage = round(difference*100, 2)
            non_zero_labels.append('Other Signatures Below Significance Threshold, {}%'.format(percentage))
            non_zero_weights.append(difference)
            colors.append('black')

        # Plot
        fig = plt.figure(figsize=(20, 7))
        fig_cols = 2
        mutation_annotation_fig_col = None
        sample_phenotype_fig_col = None
        if self.__mutation_annotation_df is not None:
            mutation_annotation_fig_col = fig_cols
            fig_cols += 1
        if self.__sample_phenotype_df is not None:
            sample_phenotype_fig_col = fig_cols
            fig_cols += 1

        # Add pie chart
        ax_pie = plt.subplot2grid((3, fig_cols), (0, 0), rowspan=2, colspan=2)
        title = 'Cancer processes Weights'
        title = '{} for {}'.format(title, sample_id)
        ax_pie.axis("equal")
        pie = ax_pie.pie(non_zero_weights, startangle=90, colors=colors)
        ax_pie.set_title(title)

        # Add pie-chart legend
        ax_legend = plt.subplot2grid((3, fig_cols), (2, 0), colspan=2)
        ax_legend.axis("off")
        ax_legend.legend(pie[0], non_zero_labels, loc="center")

        # If there is information about the mutations
        if mutation_annotation_fig_col is not None:
            # Add a subplot to describe sample information
            ax_info = plt.subplot2grid((3, fig_cols), (0, mutation_annotation_fig_col), rowspan=3)
            ax_info.axis([0, 10, 0, 10])
            ax_info.axis("off")
            ax_info.set_title("Mutations")
            info_txt = "cBioPortal"
            if sample_id in self.__mutation_annotation_df.index:
                sample_mutations = self.__mutation_annotation_df.loc[sample_id]
                sample_mutations = sample_mutations[sample_mutations.notnull()]
                if len(sample_mutations) == 0:
                    info_txt += "\n  - None"
                else:
                    info_txt += "\n" + "\n".join(map(lambda x: "   - "+x[0]+" : "+x[1], zip(sample_mutations.index, sample_mutations.values)))
            else:
                info_txt += "\n  - None"
            ax_info.text(0, 9, info_txt, fontsize=12, verticalalignment='top')

        # If phenotype data available
        if sample_phenotype_fig_col is not None:
            # Add a subplot to describe sample information
            ax_info = plt.subplot2grid((3, fig_cols), (0, sample_phenotype_fig_col), rowspan=3)
            ax_info.axis([0, 10, 0, 10])
            ax_info.axis("off")
            ax_info.set_title("Sample phenotype")
            info_txt = "Phenotype"
            if sample_id in self.__sample_phenotype_df.index:
                sample_phenotype = self.__sample_phenotype_df.loc[sample_id]
                info_txt += "\n" + "\n".join(map(lambda x: "   - "+x[0]+" : "+str(x[1]), zip(sample_phenotype.index, sample_phenotype.values)))
            else:
                info_txt += "\n  - None"
            ax_info.text(0, 9, info_txt, fontsize=12, verticalalignment='top')

        fig.canvas.set_window_title(title)

    def __save_weights_figures(self,
                               sample_id,
                               sample_vars_dict,
                               normalized_weights,
                               ):
        sample_txt = sample_id
        figure_file = join_path(self.figures_dir, sample_id+"_cancersig_profile.pdf")
        pdf_file = PdfPages(figure_file)
        # Pie chart of signature components
        self.__plot_pie_chart(normalized_weights, sample_txt)
        pdf_file.savefig(bbox_inches='tight')
        # Sample profile (Ti/Tv distribution)
        y_max = self.__plot_sample_profile(sample_txt, sample_vars_dict)
        pdf_file.savefig(bbox_inches='tight')
        # Plot sample profile after reconstructed
        reconstructed_tumor_profile = self.__plot_reconstructed_profile(sample_id,
                                                                  normalized_weights,
                                                                  y_max=y_max)
        pdf_file.savefig(bbox_inches='tight')
        # Plot the difference between the reconstructed and the normal profile
        self.__plot_difference(sample_id,
                               reconstructed_tumor_profile,
                               sample_vars_dict,
                               y_max=y_max)
        pdf_file.savefig(bbox_inches='tight')
        self.info("Saving figure to: " + figure_file)
        pdf_file.close()
        plt.clf()
        plt.close('all')

    def __plot_sample_profile(self, sample_id, sample_vars_dict):
        """
        Plot the context profile for the original tumor sample given.
        Return the maximum y value on the y-axis in order to generate
        an appropriately scaled plot.
        """
        title = 'Tumor Profile for ' + sample_id
        return plot_distribution(title, profile_dict=sample_vars_dict)

    def __plot_reconstructed_profile(self,
                                     sample_id,
                                     weights,
                                     y_max=1):
        """Given a set of weights for each signature plot the reconstructed tumor profile"""
        reconstructed_tumor_profile = self.__get_reconstructed_tumor_profile(weights)
        title = 'Reconstructed Tumor Profile for ' + sample_id
        plot_distribution(title,
                          feature_infos=self.__feature_infos,
                          feature_fractions=reconstructed_tumor_profile,
                          y_max=y_max,
                          )
        return reconstructed_tumor_profile

    def __plot_difference(self,
                          sample_id,
                          reconstructed_tumor_profile,
                          sample_vars_dict,
                          y_max):
        """Plot the difference between the original and reconstructed tumor profile."""
        reconstructed_vars_dict = get_profile_dict(self.__feature_infos, reconstructed_tumor_profile)

        differnce_vars_dict = defaultdict(dict)
        for variant_type, variant_infos in sample_vars_dict.items():
            for variant_info in variant_infos:
                original_fraction = sample_vars_dict[variant_type][variant_info]
                reconstructed_fraction = reconstructed_vars_dict[variant_type][variant_info]
                difference_fraction = original_fraction - reconstructed_fraction
                differnce_vars_dict[variant_type][variant_info] = difference_fraction
        title = 'Difference Between Original and Reconstructed Tumor Profile for ' + sample_id
        return plot_distribution(title, profile_dict=differnce_vars_dict, y_max=y_max)

    def visualize(self,
                  mutation_profiles,
                  signatures_probabilities,
                  output_dir,
                  mutation_annotation=None,
                  sample_phenotype=None,
                  ):
        # prepare data and output directories
        self.__output_root_dir = output_dir
        self.__weights_dir = None
        self.__load_signatures_probabilities(signatures_probabilities)
        if mutation_annotation is not None:
            self.__mutation_annotation_df = pandas.read_csv(mutation_annotation, sep='\t', engine='python').set_index(CASE_SUBMITTER_ID_COL_NAME)
        else:
            self.__mutation_annotation_df = None
        if sample_phenotype is not None:
            self.__sample_phenotype_df = pandas.read_csv(sample_phenotype, sep='\t', engine='python').set_index(CASE_SUBMITTER_ID_COL_NAME)
        else:
            self.__sample_phenotype_df = None
        self.__samples_normalized_weights = OrderedDict()
        if not os.path.exists(self.figures_dir):
            os.makedirs(self.figures_dir)
        # perform the actual work, finding underlying cancer processes for each sample
        self.__run_reconstruct_profiles(mutation_profiles)
        self.__write_weights(self.__samples_normalized_weights)

    def __run_reconstruct_profiles(self, mutation_profiles):
        input_file = basename(mutation_profiles)
        file_prefix, file_extention = splitext(input_file)
        self.info("Loading input file: " + mutation_profiles)
        samples_id, samples_context_ratio, samples_vars_dict = self.__load_sample_profiles(mutation_profiles)

        # Remove signatures from possibilities if they have a "strong" peak for a context that
        # is not seen in the tumor sample
        num_samples = len(samples_id)
        for sample_idx in range(num_samples):
            sample_id = samples_id[sample_idx]
            sample_context_ratio = samples_context_ratio[:, sample_idx]
            sample_vars_dict = samples_vars_dict[sample_id]

            self.info("")
            msg = " Processing sample: " + sample_id + " "
            self.info(msg.center(140, "*"))
            normalized_weights = self.__calculate_weights(sample_context_ratio)
            self.__samples_normalized_weights[sample_id] = normalized_weights
            self.__save_weights_figures(sample_id,
                                        sample_vars_dict,
                                        normalized_weights,
                                        )

    def __write_weights(self, samples_weights):
        """
        Save weights of signatures that give minimum square error.
        The file can be directly used for generating heatmap.
        """
        samples_normalized_weights = self.__samples_normalized_weights
        with open(self.weights_file, "w") as f_w:
            # write output header
            header = "\t".join(["Signature id"]+list(samples_normalized_weights.keys()))
            f_w.write(header+"\n")
            for sig_idx in range(len(self.__signature_names)):
                content = self.__signature_names[sig_idx]
                content += "\t" + "\t".join(map(lambda x: "{:14.12f}".format(samples_normalized_weights[x][sig_idx]),
                                                samples_normalized_weights,
                                                ))
                f_w.write(content+"\n")
            content = "Other signatures"
            content += "\t" + "\t".join(map(lambda x: "{:14.12f}".format(1-sum(samples_normalized_weights[x])),
                                            samples_normalized_weights,
                                            ))
            f_w.write(content+"\n")
        self.info("")
        self.info("Done!! Saving cancer processes weight to: " + self.weights_file)
