import pandas
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
from sklearn.metrics.pairwise import cosine_similarity
from matplotlib.backends.backend_pdf import PdfPages
from cancersig.template import pyCancerSigBase
from cancersig.profile.features import VARIANT_TYPE
from cancersig.profile.features import VARIANT_SUBGROUP
from cancersig.profile.features import FEATURE_ID
from cancersig.signature.figure import plot_distribution

RESAMPLING_SIZE = 100000
CONVERGENCE_CUTOFF = 0.01

ALPHABETS = "ABCDEFGHIJKLMNOPQRSTU"

DEFAULT_MIN_SIGNATURES = 2
DEFAULT_MAX_SIGNATURES = 15

DEFAULT_WEAK_SIGNAL_CUTOFF = 0.00
DEFAULT_N_ITERATIONS = 20

class DecipheredProcesses(pyCancerSigBase):
    """
    Containing all necessary information regarding the deciphered 
    cancer signature processes
    """
    def __init__(self,
                 n_mutation_types,
                 n_signatures,
                 n_iterations,
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
        self.__est_p_all = np.zeros((n_iterations, n_mutation_types, n_signatures))
        self.__est_e_all = np.zeros((n_iterations, n_signatures, len(genome_ids)))

    @property
    def avg_proc_sim(self):
        if self.__avg_proc_sim is not None:
            return self.__avg_proc_sim
        # calculate cosine similarity for each signature cluster
        sum_proc_sim = np.zeros(self.__n_signatures)
        iter_count = self.__iter_count
        for n_sig in range(self.__n_signatures):
            sum_proc_sim[n_sig] = (sum(cosine_similarity(self.__est_p_all[:iter_count,:,n_sig], [self.processes[:,n_sig]]))/iter_count)[0]
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
            header += "\t" + "\t".join(map(lambda x: "Signature "+ALPHABETS[x], range(self.__n_signatures)))
            f_p.write(header+"\n")
            # write content
            for mut_idx in range(self.__n_mutation_types):
                content = "\t".join(list(feature_infos.iloc[mut_idx, :]))
                content += "\t" + "\t".join(map(lambda x: "{:14.12f}".format(x), list(self.processes[mut_idx])))
                f_p.write(content+"\n")
        self.info()
        msg = "Deciphered processes with " + str(self.__n_signatures) + " signatures"
        msg += " have been export to " + file_name
        self.info(msg)

    def export_exposures_to_txt(self, file_name):
        exposures = self.exposures
        with open(file_name, "w") as f_e:
            # write output header
            header = "\t".join(["signature id"]+list(self.__genome_ids))
            f_e.write(header+"\n")
            n_sig, n_samples = exposures.shape
            for sig_idx in range(n_sig):
                content = "signature " + ALPHABETS[sig_idx]
                content += "\t" + "\t".join(map(lambda x: "{:14.12f}".format(x), exposures[sig_idx]))
                f_e.write(content+"\n")
        self.info()
        msg = "Deciphered exposures"
        msg += " have been export to " + file_name
        self.info(msg)

    def export_processes_to_pdf(self, file_name):
        pdf_file = PdfPages(file_name)
        for sig_idx in range(self.__n_signatures):
            title = "Signature " + ALPHABETS[sig_idx]
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
        msg += " have been export to " + file_name
        self.info(msg)


class CancerSigNMF(pyCancerSigBase):
    """
    Deciphering Signatures of Mutational Processes as described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/pdf/main.pdf
    """

    def __init__(self,
                 n_iterations=DEFAULT_N_ITERATIONS,
                 signal_cutoff=DEFAULT_WEAK_SIGNAL_CUTOFF,
                 *args,
                 **kwargs
                 ):
        self.__n_iterations = n_iterations
        self.__signaal_cutoff = signal_cutoff
        super(CancerSigNMF, self).__init__(*args, **kwargs)

    @property
    def n_mutation_types(self):
        return self.__n_mutation_types

    @property
    def n_signatures(self):
        return self.__n_signatures

    @property
    def n_iterations(self):
        return self.__n_iterations

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
                                        size=RESAMPLING_SIZE,
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
                bootstrapped_genomes.at[mutation_type, genome] = count/RESAMPLING_SIZE
        return bootstrapped_genomes

    def __remove_weak_signals(self, genomes_df):
        """ removing mutation type with weak signals """
        total_mutations = np.sum(genomes_df, axis=1)
        total_mutations = total_mutations/np.sum(total_mutations)
        new_genomes = genomes_df.drop(total_mutations[total_mutations < self.__signaal_cutoff].index)
        scaled_genomes = new_genomes/np.sum(new_genomes)
        return scaled_genomes

    def __extract_signatures(self, genomes):
        """
        Resampling genomes and using non-negative matrix factorization to decipher
        P and E matrix of each resampling

        return DecipheredProcesses object
        """
        deproc = DecipheredProcesses(self.n_mutation_types,
                                     self.n_signatures,
                                     self.n_iterations,
                                     genomes.columns.values,
                                     feature_infos=self.__mutation_profiles.iloc[:, 0:3],
                                     )
        max_iter = self.n_iterations
        self.info("Deciphering signatures processes at convergence cut-off at " + str(CONVERGENCE_CUTOFF))
        for i in range(self.n_iterations):
            bootstrapped_genomes = self.__bootstrap_genomes(genomes)
            nmf_model = NMF(n_components=self.n_signatures, tol=1e-09, max_iter=self.n_iterations)
            est_p = nmf_model.fit_transform(bootstrapped_genomes)
            est_e = nmf_model.components_
            improvement = deproc.add_stochastic_iteration(est_p/est_p.sum(axis=0)[None,:],
                                                          est_e/est_e.sum(axis=0)[None,:])
            if improvement > CONVERGENCE_CUTOFF:
                # adding another iteration is still not yet converged
                # so double the maximum iteration
                max_iter = i * 2
            if i > max_iter:
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
        self.info("--------------------------------------------------------------------------------")
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
        self.info("Extracting " + str(n_signatures) + " mutational signatures for " + str(self.n_iterations) + " iterations")
        genomes_original_df = genomes_original_df.set_index(FEATURE_ID)
        self.__genomes_original_df = genomes_original_df
        strong_genomes_df = self.__remove_weak_signals(genomes_original_df)
        self.__n_mutation_types, self.__n_genomes  = strong_genomes_df.shape
        self.info("Removing mutation types that together account for less than " + str(DEFAULT_WEAK_SIGNAL_CUTOFF*100) + "% of mutations in all genomes")
        self.info("-- After removing --")
        self.info("Total mutation types left: " + str(self.__n_mutation_types))
        self.info("Total number of samples: " + str(self.__n_genomes))
        return self.__extract_signatures(strong_genomes_df)

class CancerSigNMFController(CancerSigNMF):

    def __init__(self, *args, **kwargs):
        super(CancerSigNMFController, self).__init__(*args, **kwargs)

    def decipher(self,
                 mutation_profiles,
                 output_prefix,
                 min_signatures=DEFAULT_MIN_SIGNATURES,
                 max_signatures=DEFAULT_MAX_SIGNATURES,
                 ):
        for n_sig in range(min_signatures, max_signatures+1):
            deproc = self.decipher_mutational_processes(mutation_profiles,
                                                        n_sig)
            self.info("Signatures processes have been deciphered with reproducibility rate: " + str(deproc.avg_proc_sim))
            reconstructed_genomes = np.matmul(deproc.processes, deproc.exposures)
            reconstruction_error = np.linalg.norm(self.genomes_original_df - reconstructed_genomes)
            self.info("Reconstruction error: " + str(reconstruction_error))
            out_txt_deciphered_prociphered_processes = output_prefix + "_" + str(n_sig) + "_processes.txt"
            deproc.export_processes_to_txt(out_txt_deciphered_prociphered_processes)
            out_txt_deciphered_prociphered_exposures = output_prefix + "_" + str(n_sig) + "_exposures.txt"
            deproc.export_exposures_to_txt(out_txt_deciphered_prociphered_exposures)
            out_pdf_deciphered_prociphered_processes = output_prefix + "_" + str(n_sig) + "_processes.pdf"
            deproc.export_processes_to_pdf(out_pdf_deciphered_prociphered_processes)
        self.info("DONE!!")
