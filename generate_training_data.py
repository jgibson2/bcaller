"""
Generates training data from the 1000 Genomes VCF data
"""
import multiprocessing
import subprocess
import argparse
import vcf
import numpy as np

def process_vcf(fname):
    subprocess.check_call(["bcaller/bin/vcf_to_hdf5", fname, fname + "_data.hdf5"], shell=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '---vcf-files', type=str, nargs='+', help='Input VCF files', action='store', required=True)
    parser.add_argument('--merge-only', action='store_true', required=False)
    parser.add_argument('-o', '--out-file', type=str, help='Output file', action='store', required=True)
    parser.add_argument('-w', '--window-size', type=int, help='Window size', action='store', default=300, required=False)

    args = parser.parse_args()

    count = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=count)
    if not args.merge_only:
        pool.map(process_vcf, args.vcf_files)

    for vcf_file in args.vcf_files:
        reader = vcf.Reader(filename=vcf_file)
        contigs = list(reader.contigs.keys())
        for contig in contigs:
            subprocess.check_call(["bcaller/bin/hdf5_to_window_frequencies", vcf_file + "_data.hdf5", args.out_file, contig, str(args.window_size)], shell=False)


