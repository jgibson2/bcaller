import sys
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd

plt.rcParams["figure.figsize"] = (20,6)

infile = sys.argv[1]
outfile = sys.argv[2]

gtruthfile = None
if len(sys.argv) == 4:
    gtruthfile = sys.argv[3]

data = pd.read_csv(infile, sep='\t', header=None, names=('chr', 'start', 'end', 'ref', 'alt', 'prob'))

gtruth = None
if gtruthfile is not None:
    gtruth = pd.read_csv(gtruthfile, sep='\t', header=None, names=('pos', 'ref', 'alt'))


x = data.start.values
y = data.prob.values

colors = ['blue' for pt in x]
if gtruth is not None:
    q = set(gtruth.pos.values)
    for i,pt in enumerate(x):
        if pt in q:
            colors[i] = 'red'

plt.scatter(x, y, c=colors)
plt.title('Probabilities of SNPs over chromosome 1')
plt.ylabel('Probability of SNP')
plt.xlabel('Genomic coordinate')
if gtruth is not None:
    legend_elements = [Line2D([0], [0], marker='o', color='blue', label='Sequencing Errors', markerfacecolor='blue', markersize=5), Line2D([0], [0], marker='o', color='red', label='True SNPs',markerfacecolor='red', markersize=5)]
    plt.legend(handles=legend_elements)
plt.savefig(outfile)
