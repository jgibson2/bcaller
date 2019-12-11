"""
Ref: https://github.com/cornellius-gp/gpytorch/blob/master/examples/04_Scalable_GP_Regression_1D/KISSGP_Regression_1D.ipynb
"""

import h5py
import argparse
import matplotlib
import numpy as np
import scipy.special
from matplotlib import pyplot as plt
matplotlib.use('TkAgg')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '---in-file', type=str, help='Input HDF5 file', action='store', required=True)
    parser.add_argument('-c', '--contigs', type=str, nargs='+', help='Contigs to train on', action='store', required=False)

    args = parser.parse_args()

    input = h5py.File(args.in_file, 'r')
    if not args.contigs:
        args.contigs = list(input.keys())

    models = {}
    for contig in args.contigs:
        train_data = input[contig][:]
        # int limits...
        train_x = train_data[:, 0]
        # train_y = scipy.special.logit(np.array(train_data[:, 1]))
        train_y = np.log(train_data[:, 1])

        plt.clf()
        plt.plot(train_x, train_y)
        plt.title('Mutation probabilities over genetic coordinate for {} with window size {}'.format(contig, input.attrs['/{}/window_size'.format(contig)]))
        plt.xlabel('Genomic coordinate')
        plt.ylabel('log(Probability)')
        plt.show()