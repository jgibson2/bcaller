"""
Call variants
"""
import multiprocessing
import functools
import argparse
import vcf
import numpy as np
import itertools
import scipy.special
from scipy.integrate import quad
from scipy.stats.distributions import beta, binom
import h5py

def generate_prior_SNP_frequency_distributions(prior_data, window_size, avg_coverage=20):
    prior_coords = np.uint64(prior_data[:, 0])
    prior_alpha = (prior_data[:, 1] * avg_coverage)
    prior_beta = avg_coverage - prior_alpha
    return prior_coords, prior_alpha, prior_beta

def generate_prior_mutation_frequency_distribution(prior_data_file):
    bases = {'A', 'T', 'C', 'G'}
    mut_freq_dict = {}
    permutations = itertools.permutations(bases, 2)
    for b1, b2 in permutations:
        if b1 != b2:
            freq = prior_data_file['/freqs/{}/{}'.format(b1, b2)].value
            mut_freq_dict[(b1, b2)] = freq
    # scaling of frequencies
    for b1 in bases:
        total = np.sum([mut_freq_dict[(b1, b2)] for b2 in bases if b2 != b1])
        for b2 in bases:
            if b1 != b2:
                mut_freq_dict[(b1, b2)] /= total
    return mut_freq_dict

def p_theta_given_X_B1_B2_beta_parameters(prior_coords, prior_alpha, prior_beta, window_size, mutation_frequencies, X, B1, B2):
    window = int(X // window_size)
    if window < prior_coords.shape[0]:
        assert(prior_coords[window] <= X < prior_coords[window] + window_size)
        return ((prior_alpha[window] * mutation_frequencies[(B1, B2)]) + 0.5, prior_beta[window] + (prior_alpha[window] * (1.0 - mutation_frequencies[(B1, B2)]))+ 0.5)
    return (0.5, 0.5)

def p_theta_given_D_X_B1_B2_beta_parameters(B2_in_D, total_D, prior_coords, prior_alpha, prior_beta, window_size, mutation_frequencies, X, B1, B2):
    prior_alpha, prior_beta = p_theta_given_X_B1_B2_beta_parameters(prior_coords, prior_alpha, prior_beta, window_size, mutation_frequencies, X, B1, B2)
    return (prior_alpha + B2_in_D, prior_beta + total_D - B2_in_D)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '---in-file', type=str, help='Input HDF5 file', action='store', required=True)
    parser.add_argument('-o', '---out-file', type=str, help='Output BED File', action='store', required=True)
    parser.add_argument('-v', '--vcf', type=str, help='Input VCF file', action='store', required=True)
    parser.add_argument('--hdf5-contigs', type=str, nargs='+', help='Ordered list of contigs to use prior from in input HDF5 file', action='store', required=False)
    parser.add_argument('--vcf-contigs', type=str, nargs='+', help='Ordered list of contigs to predict on in input VCF file', action='store', required=False)
    parser.add_argument('--coverage', type=int, default=20, help='Average coverage', action='store')
    parser.add_argument('--min-probability', type=float, default=0.0, help='Minimum probability to output', action='store')

    args = parser.parse_args()

    input = h5py.File(args.in_file, 'r')
    if not args.hdf5_contigs:
        args.hdf5_contigs = list(sorted([k for k in input.keys() if k != 'freqs']))
    if not args.vcf_contigs:
        args.vcf_contigs = args.hdf5_contigs

    assert(len(args.vcf_contigs) == len(args.hdf5_contigs))

    reader = vcf.Reader(filename=args.vcf)

    mutation_type_prior = generate_prior_mutation_frequency_distribution(input)

    outfile = open(args.out_file, 'w')

    for hdf5_contig, vcf_contig in zip(args.hdf5_contigs, args.vcf_contigs):
        records = reader.fetch(vcf_contig)
        prior_data = input[hdf5_contig][:]
        window_size = input.attrs['/{}/window_size'.format(hdf5_contig)]

        prior_coords, prior_alpha, prior_beta = generate_prior_SNP_frequency_distributions(prior_data, window_size, args.coverage)

        for rec in records:
            ref = rec.REF
            alts = rec.ALT
            total_depth = rec.INFO['DP']
            base_depths = rec.INFO['AD']
            alt_and_depths = zip(alts, base_depths[1:])
            for alt, dp in alt_and_depths:
                if str(alt) not in {'A', 'T', 'C', 'G'} or str(ref) not in {'A', 'T', 'C', 'G'}:
                    continue
                # used 0-based coords for windows
                a, B = p_theta_given_D_X_B1_B2_beta_parameters(dp, total_depth, prior_coords, prior_alpha, prior_beta, window_size, mutation_type_prior, rec.POS - 1, str(ref), str(alt))
                prob = 1.0E-8
                if a + B - 2 > 0:
                    prob = (a - 1) / (a + B - 2)
                if prob >= args.min_probability:
                    outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(vcf_contig, rec.POS, rec.POS + 1, ref, alt, prob))
    print('DONE')