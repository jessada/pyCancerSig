#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import numpy as np
import pandas
import itertools
from collections import defaultdict
from collections import OrderedDict
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import ticker
from matplotlib.font_manager import FontProperties
from sklearn.decomposition import NMF
from sklearn.metrics.pairwise import cosine_similarity
from sigfig import plot_distribution

FEATURES_COL_NAME = "feature_id"
VARIANT_INFO_COL_NAME = "variant info"
VARIANT_TYPE_COL_NAME = "variant type"
RESAMPLING_SIZE = 100000
CONVERGENCE_CUTOFF = 0.01

DEFAULT_WEAK_SIGNAL_CUTOFF = 0.00
DEFAULT_N_ITERATIONS = 20

REPRODUCIBITILITY_RATE = "reproducibility rate"
RECONSTRUCTION_ERROR = "reconstruction error"

class DecipheredProcesses:
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
                 ):
        self.__n_mutation_types = n_mutation_types
        self.__n_signatures = n_signatures
        self.__genome_ids = genome_ids
        self.__feature_infos = feature_infos
        self.__iter_count = 0
        self.__centroids = None
        self.__exposures = None
        self.__avg_proc_sim = None

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
    def centroids(self):
        return self.__centroids

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
            header += "\t" + "\t".join(map(lambda x: "Signature "+str(x+1), range(self.__n_signatures)))
            f_p.write(header+"\n")
            # write content
            for mut_idx in range(self.__n_mutation_types):
                content = "\t".join(list(feature_infos.iloc[mut_idx, :]))
                content += "\t" + "\t".join(map(lambda x: "{:14.12f}".format(x), list(self.processes[mut_idx])))
                f_p.write(content+"\n")
        print("", file=sys.stderr)
        msg = "Deciphered processes with " + str(self.__n_signatures) + " signatures"
        msg += " have been export to " + file_name
        print(">> "+msg, file=sys.stderr)

    def export_exposures_to_txt(self, file_name):
        exposures = self.exposures
        with open(file_name, "w") as f_e:
            # write output header
            header = "\t".join(["signature id"]+list(self.__genome_ids))
            f_e.write(header+"\n")
            n_sig, n_samples = exposures.shape
            for sig_idx in range(n_sig):
                content = "signature " + str(sig_idx+1)
                content += "\t" + "\t".join(map(lambda x: "{:14.12f}".format(x), exposures[sig_idx]))
                f_e.write(content+"\n")
        print("", file=sys.stderr)
        msg = "Deciphered exposures"
        msg += " have been export to " + file_name
        print(">> "+msg, file=sys.stderr)

    def export_processes_to_pdf(self, file_name):
        pdf_file = PdfPages(file_name)
        for sig_idx in range(self.__n_signatures):
            title = "Signature #" + str(sig_idx+1)
            plot_distribution(title=title,
                              feature_infos=self.__feature_infos,
                              feature_fractions=self.processes[:, sig_idx],
                              )
            pdf_file.savefig(bbox_inches='tight')

        pdf_file.close()
        plt.clf()
        plt.close('all')
        print("", file=sys.stderr)
        msg = "Deciphered processes with " + str(self.__n_signatures) + " signatures"
        msg += " have been export to " + file_name
        print(">> "+msg, file=sys.stderr)


class CancerSig:
    """
    Deciphering Signatures of Mutational Processes as described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/pdf/main.pdf
    """
    def __init__(self,
                 n_iterations=DEFAULT_N_ITERATIONS,
                 signal_cutoff=DEFAULT_WEAK_SIGNAL_CUTOFF,
                 ):
        self.__n_iterations = n_iterations
        self.__signaal_cutoff = signal_cutoff

    @property
    def n_mutation_types(self):
        return self.__n_mutation_types

    @property
    def n_signatures(self):
        return self.__n_signatures

    @property
    def n_genomes(self):
        return self.__n_genomes

    @property
    def n_iterations(self):
        return self.__n_iterations

    @property
    def genomes_original_df(self):
        return self.__genomes_original_df

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
                if count:
                    count = float(count)
                else:
                    count = 0.0
                bootstrapped_genomes.at[mutation_type, genome] = count/RESAMPLING_SIZE
        return bootstrapped_genomes

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
        print(">> Deciphering signatures processes at convergence cut-off at " + str(CONVERGENCE_CUTOFF), file=sys.stderr)
        for i in range(self.n_iterations):
            bootstrapped_genomes = self.__bootstrap_genomes(genomes)
            nmf_model = NMF(n_components=self.n_signatures, tol=1e-09, max_iter=10000)
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
        print(">> Deciphering was done at iteration " + str(i), file=sys.stderr)
        return deproc

    def __remove_weak_signals(self, genomes_df):
        """ removing mutation type with weak signals """
        total_mutations = np.sum(genomes_df, axis=1)
        total_mutations = total_mutations/np.sum(total_mutations)
        new_genomes = genomes_df.drop(total_mutations[total_mutations < self.__signaal_cutoff].index)
        scaled_genomes = new_genomes/np.sum(new_genomes)
        return scaled_genomes

    def __input_error_n_exit(self, variable_name):
            print("", file=sys.stderr)
            print("ERROR: Please specify an input file containing the variables: " + variable_name,
                  file=sys.stderr)
            print("", file=sys.stderr)
            exit(1)

    def decipher_mutational_processes(self,
                                      mutation_catalog_file,
                                      n_signatures,
                                      ):
        """
        This function will load a mutation matrix (G) from mutation_catalog_file
        and attempt to decipher the signatures matrix (P) and the corresponding
        exposures matrix (E), that G = PxE.

        n_signatures is the number of signatures expected to be deciphered, aka.
        number of rows in matrix P
        """
        print("-------------------------------------------------", file=sys.stderr)
        print(">> Deciphering signature processes from : " + mutation_catalog_file, file=sys.stderr)
        self.__n_signatures = n_signatures
        # load samples profile
        self.__mutation_profiles = pandas.read_csv(mutation_catalog_file,
                                                   sep='\t',
                                                   engine='python')
        genomes_original_df = self.__mutation_profiles.drop([VARIANT_INFO_COL_NAME, VARIANT_TYPE_COL_NAME], axis=1)
        if FEATURES_COL_NAME not in self.__mutation_profiles:
            self.__input_error_n_exit(FEATURES_COL_NAME)
        if VARIANT_INFO_COL_NAME not in self.__mutation_profiles:
            self.__input_error_n_exit(VARIANT_INFO_COL_NAME)
        if VARIANT_TYPE_COL_NAME not in self.__mutation_profiles:
            self.__input_error_n_exit(VARIANT_TYPE_COL_NAME)

        msg = ">> Extracting " + str(n_signatures) + " mutational signatures"
        msg += " for " + str(self.n_iterations) + " iterations"
        print("", file=sys.stderr)
        print(msg, file=sys.stderr)
        genomes_original_df = genomes_original_df.set_index(FEATURES_COL_NAME)
        self.__genomes_original_df = genomes_original_df
        strong_genomes_df = self.__remove_weak_signals(genomes_original_df)
        self.__n_mutation_types, self.__n_genomes  = strong_genomes_df.shape
        print(">> Removing mutation types that together account for less than " + str(DEFAULT_WEAK_SIGNAL_CUTOFF*100) + "% of mutations in all genomes",  file=sys.stderr)
        print(">> -- After removing --",  file=sys.stderr)
        print(">> Total mutation types left: " + str(self.__n_mutation_types), file=sys.stderr)
        print(">> Total number of samples: " + str(self.__n_genomes), file=sys.stderr)
        return self.__extract_signatures(strong_genomes_df)

if __name__ == '__main__':
    parser = argparse.ArgumentParser("""Evaluating Non-negative matrix factorization""")

    parser.add_argument('--mutation_profiles',
                        dest='mutation_profiles',
                        help='input mutation calalog to be deciphered',
                        required=True,
                        )
    parser.add_argument('--min_signatures',
                        dest='min_signatures',
                        type=int,
                        help='minimum number of signatures to be deciphered',
                        default=2,
                        )
    parser.add_argument('--max_signatures',
                        dest='max_signatures',
                        type=int,
                        help='maximum number of signatures to be deciphered',
                        default=15,
                        )
    parser.add_argument('--out_prefix',
                        dest='out_prefix',
                        help='output prefix',
                        required=True,
                        )
    args = parser.parse_args()
    print("", file=sys.stderr)
    print(">> Input mutation profiles: " + args.mutation_profiles, file=sys.stderr)
    print(">> Minimum number of signatures to be deciphered: " + str(args.min_signatures), file=sys.stderr)
    print(">> Maximum number of signatures to be deciphered: " + str(args.max_signatures), file=sys.stderr)
    print(">> Output file prefix: " + args.out_prefix, file=sys.stderr)
    print("", file=sys.stderr)
    stats = {}
    for n_sig in range(args.min_signatures, args.max_signatures+1):
        stat = {}
        cancersig = CancerSig(n_iterations=1000)
        deproc = cancersig.decipher_mutational_processes(mutation_catalog_file=args.mutation_profiles,
                                                         n_signatures=n_sig)
        reproducibility_rate = deproc.avg_proc_sim
        print(">> Signatures processes have been deciphered with reproducibility rate: " + str(reproducibility_rate), file=sys.stderr)
        stat[REPRODUCIBITILITY_RATE] = reproducibility_rate
        reconstructed_genomes = np.matmul(deproc.processes, deproc.exposures)
        reconstruction_error = np.linalg.norm(cancersig.genomes_original_df - reconstructed_genomes)
        print(">> Reconstruction error: " + str(reconstruction_error), file=sys.stderr)
        stat[RECONSTRUCTION_ERROR] = reconstruction_error
        out_txt_deciphered_prociphered_processes = args.out_prefix
        out_txt_deciphered_prociphered_processes += "_" + str(n_sig) + "_processes.txt"
        deproc.export_processes_to_txt(out_txt_deciphered_prociphered_processes)
        out_txt_deciphered_prociphered_exposures = args.out_prefix
        out_txt_deciphered_prociphered_exposures += "_" + str(n_sig) + "_exposures.txt"
        deproc.export_exposures_to_txt(out_txt_deciphered_prociphered_exposures)
        out_pdf_deciphered_prociphered_processes = args.out_prefix
        out_pdf_deciphered_prociphered_processes += "_" + str(n_sig) + "_processes.pdf"
        deproc.export_processes_to_pdf(out_pdf_deciphered_prociphered_processes)
        stats[n_sig] = stat
    
    print("-------------------------------------------------", file=sys.stderr)
    stat_file = args.out_prefix + ".stat" 
    with open(stat_file, "w") as f_s:
        header = "n_signatures\treproducibility_rate\treconstruction_error"
        f_s.write(header + "\n")
        for n_sig in range(args.min_signatures, args.max_signatures+1):
            content = str(n_sig)
            content += "\t" + str(stats[n_sig][REPRODUCIBITILITY_RATE])
            content += "\t" + str(stats[n_sig][RECONSTRUCTION_ERROR])
            f_s.write(content + "\n")
    msg = "Model statistics have been export to " + stat_file
    print("", file=sys.stderr)
    print(">> "+msg, file=sys.stderr)
    print("", file=sys.stderr)
    print("-------------------------------------------------", file=sys.stderr)

    print("", file=sys.stderr)
    print(">> DONE!!", file=sys.stderr)
    print("", file=sys.stderr)
