import pandas
import numpy as np
import matplotlib.pyplot as plt
import string
from sklearn.decomposition import NMF
from sklearn.metrics.pairwise import cosine_similarity
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PatchCollection
from cancersig.template import pyCancerSigBase
from cancersig.utils.disp import disp_subparam
from cancersig.utils.disp import disp_header
from cancersig.profile.features import VARIANT_TYPE
from cancersig.profile.features import VARIANT_SUBGROUP
from cancersig.profile.features import FEATURE_ID
from cancersig.signature.figure import plot_distribution

DEFAULT_RESAMPLING_SIZE = 100000
DEFAULT_CONVERGENCE_CUTOFF = 0.01

SIGNATURE_INDEX_CODES = string.ascii_uppercase

DEFAULT_MIN_SIGNATURES = 2
DEFAULT_MAX_SIGNATURES = 15

DEFAULT_WEAK_SIGNAL_CUTOFF = 0.00
DEFAULT_MAX_MODEL_SOLUTIONS = 200

class DecipheredProcesses(pyCancerSigBase):
    """
    Containing all necessary information regarding the deciphered 
    cancer signature processes
    """
    def __init__(self,
                 n_mutation_types,
                 n_signatures,
                 max_model_solutions,
                 genome_ids,
                 feature_infos=None,
                 *args,
                 **kwargs
                 ):
        self.__n_mutation_types = n_mutation_types
        self.__n_signatures = n_signatures
        self.__genome_ids = genome_ids
        self.__feature_infos = feature_infos
        self.__iter_count = 0
        self.__centroids = None
        self.__exposures = None
        self.__avg_proc_sim = None
        super(DecipheredProcesses, self).__init__(*args, **kwargs)

        # allocate maximum space for the calculation but only a fraction of them
        # are used
        self.__est_p_all = np.zeros((max_model_solutions, n_mutation_types, n_signatures))
        self.__est_e_all = np.zeros((max_model_solutions, n_signatures, len(genome_ids)))

    @property
    def avg_proc_sim(self):
        if self.__avg_proc_sim is not None:
            return self.__avg_proc_sim
        # calculate cosine similarity for each signature cluster
        sum_proc_sim = np.zeros(self.__n_signatures)
        iter_count = self.__iter_count
        for n_sig in range(self.__n_signatures):
            sum_proc_sim[n_sig] = (sum(cosine_similarity(self.__est_p_all[:iter_count,:,n_sig], [self.__centroids[:,n_sig]]))/iter_count)[0]
        self.__avg_proc_sim = sum(sum_proc_sim)/self.__n_signatures
        return self.__avg_proc_sim

    @property
    def processes(self):
        return self.__centroids

    @property
    def exposures(self):
        return self.__exposures

    def add_stochastic_iteration(self, P, E):
        """
        Add another set of P and E to the object and calculate the difference
        between the new centriods and the old ones, aka. improvement, using Forbenius norm
        """
        self.__est_p_all[self.__iter_count] = P
        self.__est_e_all[self.__iter_count] = E
        self.__iter_count += 1
        if self.__centroids is None:
            # for adding P and E for the first time
            # just do the initialization
            self.__centroids = P
            self.__exposures = E
            return 1000
        new_centroids = ((self.__centroids * (self.__iter_count-1)) + P) / self.__iter_count
        self.__exposures = ((self.__exposures * (self.__iter_count-1)) + E) / self.__iter_count
        improvement = np.linalg.norm(new_centroids - self.__centroids)
        self.__centroids = new_centroids
        return improvement

    def export_processes_to_txt(self, file_name):
        feature_infos = self.__feature_infos
        with open(file_name, "w") as f_p:
            # write output header
            header = "\t".join(feature_infos.columns.values)
            header += "\t" + "\t".join(map(lambda x: "Signature "+SIGNATURE_INDEX_CODES[x], range(self.__n_signatures)))
            f_p.write(header+"\n")
            # write content
            for mut_idx in range(self.__n_mutation_types):
                content = "\t".join(list(feature_infos.iloc[mut_idx, :]))
                content += "\t" + "\t".join(map(lambda x: "{:14.12f}".format(x), list(self.processes[mut_idx])))
                f_p.write(content+"\n")
        self.info()
        msg = "Deciphered processes with " + str(self.__n_signatures) + " signatures"
        msg += " have been exported to " + file_name
        self.info(msg)

    def export_exposures_to_txt(self, file_name):
        exposures = self.exposures
        with open(file_name, "w") as f_e:
            # write output header
            header = "\t".join(["signature id"]+list(self.__genome_ids))
            f_e.write(header+"\n")
            n_sig, n_samples = exposures.shape
            for sig_idx in range(n_sig):
                content = "signature " + SIGNATURE_INDEX_CODES[sig_idx]
                content += "\t" + "\t".join(map(lambda x: "{:14.12f}".format(x), exposures[sig_idx]))
                f_e.write(content+"\n")
        self.info()
        msg = "Deciphered exposures"
        msg += " have been exported to " + file_name
        self.info(msg)

    def export_processes_to_pdf(self, file_name):
        pdf_file = PdfPages(file_name)
        for sig_idx in range(self.__n_signatures):
            title = "Signature " + SIGNATURE_INDEX_CODES[sig_idx]
            plot_distribution(title=title,
                              feature_infos=self.__feature_infos,
                              feature_fractions=self.processes[:, sig_idx],
                              )
            pdf_file.savefig(bbox_inches='tight')

        pdf_file.close()
        plt.clf()
        plt.close('all')
        self.info()
        msg = "Deciphered processes with " + str(self.__n_signatures) + " signatures"
        msg += " have been exported to " + file_name
        self.info(msg)


class CancerSigNMF(pyCancerSigBase):
    """
    Deciphering Signatures of Mutational Processes as described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/pdf/main.pdf
    """
    def __init__(self,
                 resampling_size=DEFAULT_RESAMPLING_SIZE,
                 convergence_cutoff=DEFAULT_CONVERGENCE_CUTOFF,
                 weak_signal_cutoff=DEFAULT_WEAK_SIGNAL_CUTOFF,
                 max_model_solutions=DEFAULT_MAX_MODEL_SOLUTIONS,
                 *args,
                 **kwargs
                 ):
        self.__resampling_size = resampling_size
        self.__convergence_cutoff = convergence_cutoff
        self.__weak_signaal_cutoff = weak_signal_cutoff
        self.__max_model_solutions = max_model_solutions
        super(CancerSigNMF, self).__init__(*args, **kwargs)

    @property
    def n_mutation_types(self):
        return self.__n_mutation_types

    @property
    def n_signatures(self):
        return self.__n_signatures

    @property
    def resampling_size(self):
        return self.__resampling_size

    @property
    def convergence_cutoff(self):
        return self.__convergence_cutoff

    @property
    def weak_signal_cutoff(self):
        return self.__weak_signaal_cutoff

    @property
    def max_model_solutions(self):
        return self.__max_model_solutions

    @property
    def genomes_original_df(self):
        return self.__genomes_original_df

    def __input_error_n_exit(self, variable_name):
        self.throw("ERROR: Please specify an input file containing the variables: " + variable_name)
        exit(1)

    def __bootstrap_genomes(self, genomes):
        mutation_types = genomes.index.get_values()
        bootstrapped_genomes = genomes.copy(deep=True)
        for genome in genomes.columns.values:
            # Resampling the probabllity of each mutation in each genome
            #
            # The probability to observe each mutation in the sampling is
            # the correponding one from the original genome
            sampling = np.random.choice(mutation_types,
                                        size=self.resampling_size,
                                        replace=True,
                                        p=genomes[genome],
                                        )
            mutations, mutations_count = np.unique(sampling, return_counts=True)
            for mutation_type in mutation_types:
                count = mutations_count[np.where(mutations==mutation_type)]
                if count.size > 0:
                    count = float(count)
                else:
                    count = 0.0
                bootstrapped_genomes.at[mutation_type, genome] = count/self.resampling_size
        return bootstrapped_genomes

    def __remove_weak_signals(self, genomes_df):
        """
        Removing mutation type with weak signals
        """
        # sum for each mutation type
        total_mutations = np.sum(genomes_df, axis=1)
        # then divide the sum by the sum of them self (the sum of this result will be one)
        total_mutations = total_mutations/np.sum(total_mutations)
        self.__weak_signal_mutations = total_mutations[total_mutations < self.weak_signal_cutoff].index
        new_genomes = genomes_df.drop(self.__weak_signal_mutations)
        scaled_genomes = new_genomes/np.sum(new_genomes)
        return scaled_genomes

    def __extract_signatures(self, genomes):
        """
        Resampling genomes and using non-negative matrix factorization to decipher
        P and E matrix of each resampling

        return DecipheredProcesses object
        """
        deproc = DecipheredProcesses(genomes.shape[0],
                                     self.n_signatures,
                                     self.max_model_solutions,
                                     genomes.columns.values,
                                     feature_infos=self.__mutation_profiles.iloc[:, 0:3],
                                     )
        # iterate until the NMF solution doesn't get better or until it reach the max tries
        max_iter = 1
        for i in range(1, self.max_model_solutions+1):
            bootstrapped_genomes = self.__bootstrap_genomes(genomes)
            nmf_model = NMF(n_components=self.n_signatures, tol=1e-09)
            est_p = nmf_model.fit_transform(bootstrapped_genomes)
            est_e = nmf_model.components_
            improvement = deproc.add_stochastic_iteration(est_p/est_p.sum(axis=0)[None,:],
                                                          est_e/est_e.sum(axis=0)[None,:])
            if improvement > self.convergence_cutoff:
                # adding another iteration is still not yet converged
                # so double the maximum iteration
                max_iter = (i+1) * 2
            if i >= max_iter:
                break
        self.info("Deciphering was done at iteration " + str(i))
        return deproc

    def decipher_mutational_processes(self,
                                      mutation_profiles,
                                      n_signatures,
                                      ):
        """
        This function will load a mutation matrix (G) from mutation_profiles
        and attempt to decipher the signatures matrix (P) and the corresponding
        exposures matrix (E), that G = PxE.

        n_signatures is the number of signatures expected to be deciphered, aka.
        number of rows in matrix P
        """
        status_txt = " Computing " + str(n_signatures) + " optimal cancer signature processes "
        self.info(status_txt.center(120, "-"))
        self.__n_signatures = n_signatures
        # load samples profile
        self.__mutation_profiles = pandas.read_csv(mutation_profiles,
                                                   sep='\t',
                                                   engine='python')
        genomes_original_df = self.__mutation_profiles.drop([VARIANT_TYPE, VARIANT_SUBGROUP], axis=1)
        if FEATURE_ID not in self.__mutation_profiles:
            self.__input_error_n_exit(FEATURE_ID)
        if VARIANT_TYPE not in self.__mutation_profiles:
            self.__input_error_n_exit(VARIANT_TYPE)
        if VARIANT_SUBGROUP not in self.__mutation_profiles:
            self.__input_error_n_exit(VARIANT_SUBGROUP)
        genomes_original_df = genomes_original_df.set_index(FEATURE_ID)
        self.__genomes_original_df = genomes_original_df.div(genomes_original_df.sum(axis=0))
        self.__n_mutation_types, self.__n_genomes  = self.__genomes_original_df.shape
        strong_genomes_df = self.__remove_weak_signals(genomes_original_df)
        self.info("Removing mutation types that together account for less than " + str(self.__weak_signaal_cutoff*100) + "% of mutations in all genomes")
        self.info("-- After removing --")
        self.info("Total mutation types left: " + str(strong_genomes_df.shape[0]))
        self.info("Total number of samples: " + str(self.__n_genomes))
        return self.__extract_signatures(strong_genomes_df)

class CancerSigNMFController(CancerSigNMF):

    def __init__(self, *args, **kwargs):
        super(CancerSigNMFController, self).__init__(*args, **kwargs)

    def __save_performace_figure(self,
                                 reproducibility_rates,
                                 reconstruction_errors,
                                 output_prefix):
        SIG_REP_CO_LNAME = "Signature reproducibility"
        RECON_ERROR_COL_NAME = "Frobenius reconstruction error"
        df = pandas.DataFrame({"Number_of_processes": range(2, len(reproducibility_rates)+2),
                               SIG_REP_CO_LNAME: reproducibility_rates,
                               RECON_ERROR_COL_NAME: reconstruction_errors,
                               })
        df = df.set_index("Number_of_processes")
        # find optimal performance
        optimal_value = 0
        optimal_idx = 0
        for idx, row in df.iterrows():
            reproducibility_rate = row[SIG_REP_CO_LNAME]
            reconstruction_error = row[RECON_ERROR_COL_NAME]
            performance = reproducibility_rate - reconstruction_error
            if performance > optimal_value:
                optimal_value = performance
                optimal_idx = idx
        # plot model performance
        plt.style.use('classic')
        fig, ax = plt.subplots()
        ax.plot(df.index.values, df[SIG_REP_CO_LNAME], '-r', marker="o", label=SIG_REP_CO_LNAME)
        ax.plot(df.index.values, df[RECON_ERROR_COL_NAME], '-b', marker="s", label=RECON_ERROR_COL_NAME)
        ax.set_xlim([min(df.index.values)-0.5, max(df.index.values)+0.5])
        min_ylim = min(min(df[SIG_REP_CO_LNAME]), min(df[RECON_ERROR_COL_NAME]))
        max_ylim = max(max(df[SIG_REP_CO_LNAME]), max(df[RECON_ERROR_COL_NAME]))
        diff_ylim = max_ylim - min_ylim
        ylim_gap = diff_ylim * 0.05
        min_ylim = min_ylim - ylim_gap
        max_ylim = max_ylim + ylim_gap
        ax.set_ylim([min_ylim, max_ylim])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        plt.xticks(np.arange(2, max(df.index.values)+1, step=1))
        ax.set_xlabel("Number of mutational signatures")
        optimal_reproducibility_rate = df.loc[optimal_idx][SIG_REP_CO_LNAME]
        optimal_reconstruction_error = df.loc[optimal_idx][RECON_ERROR_COL_NAME]
        points = [[optimal_idx-0.2, optimal_reconstruction_error-ylim_gap*0.5],
                  [optimal_idx-0.2, optimal_reproducibility_rate+ylim_gap*0.5],
                  [optimal_idx+0.2, optimal_reproducibility_rate+ylim_gap*0.5],
                  [optimal_idx+0.2, optimal_reconstruction_error-ylim_gap*0.5]]
        patch = plt.Polygon(points, facecolor="lightcyan", linestyle=":")
        ax.add_patch(patch)
        ax.annotate("optimal solution",
                    xy=(optimal_idx+0.3, optimal_reproducibility_rate+ylim_gap*0.5),
                    xytext=(optimal_idx+0.5, optimal_reproducibility_rate+ylim_gap*2.5),
                    va="center",
                    ha="center",
                    arrowprops=dict(facecolor="black", shrink=0.05))
        ax.legend(loc="center right")
        output_figure_file = output_prefix + "_model_performance.pdf"
        fig.savefig(output_figure_file)
        output_performance_txt = output_prefix + "_model_performance.txt"
        df.to_csv(output_performance_txt, sep="\t")
        self.info()
        self.info("Model performance figure has been saved to: " + output_figure_file)
        self.info("Model performance has been exported to: " + output_performance_txt)

    def decipher(self,
                 mutation_profiles,
                 output_prefix,
                 min_signatures=DEFAULT_MIN_SIGNATURES,
                 max_signatures=DEFAULT_MAX_SIGNATURES,
                 save_signature_pdfs=True,
                 ):
        disp_header("Model parameters")
        disp_subparam("Resampling size", self.resampling_size)
        disp_subparam("Convergence cutoff", self.convergence_cutoff)
        disp_subparam("Maximum model solutions", self.max_model_solutions)
        self.info()

        reproducibility_rates = []
        reconstruction_errors = []
        for n_sig in range(min_signatures, max_signatures+1):
            deproc = self.decipher_mutational_processes(mutation_profiles,
                                                        n_sig)
            self.info("Signatures processes have been deciphered with reproducibility rate: " + str(deproc.avg_proc_sim))
            reproducibility_rates.append(deproc.avg_proc_sim)
            reconstructed_genomes = np.matmul(deproc.processes, deproc.exposures)
            reconstruction_error = np.linalg.norm(self.genomes_original_df - reconstructed_genomes)
            reconstruction_errors.append(reconstruction_error)
            self.info("Reconstruction error: " + str(reconstruction_error))
            out_txt_deciphered_prociphered_processes = output_prefix + "_" + str(n_sig) + "_processes.txt"
            deproc.export_processes_to_txt(out_txt_deciphered_prociphered_processes)
            out_txt_deciphered_prociphered_exposures = output_prefix + "_" + str(n_sig) + "_exposures.txt"
            deproc.export_exposures_to_txt(out_txt_deciphered_prociphered_exposures)
            if save_signature_pdfs:
                out_pdf_deciphered_prociphered_processes = output_prefix + "_" + str(n_sig) + "_processes.pdf"
                deproc.export_processes_to_pdf(out_pdf_deciphered_prociphered_processes)
        self.__save_performace_figure(reproducibility_rates,
                                      reconstruction_errors,
                                      output_prefix,
                                      )
        self.info()
        self.info("DONE!!")
