import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from tqdm.notebook import tqdm
import math
from scipy.stats import hypergeom
from scipy.stats import fisher_exact

neutral_gene_ratios = {}
for gl, changes in llp.items():
    #print(changes)
    gene = changes['gid']
    if gene not in neutral_gene_ratios:
        neutral_gene_ratios[gene] = [0,0]
    for b in changes['syn']:
        if len(b) > 0:
            mutation = (changes['ref'], b)
            mutation_probability = norm_mttypes[mutation]
            neutral_gene_ratios[gene][0] += mutation_probability
    for b in changes['non']:
        if len(b) > 0:
            mutation = (changes['ref'], b)
            mutation_probability = norm_mttypes[mutation]
            neutral_gene_ratios[gene][1] += mutation_probability
neutral_gene_ratios


gnr = {}
for g,ps in neutral_gene_ratios.items():
    #print(g, ps[1]/ps[0])
    gnr[g] = ps[1]/ps[0]
gnr

mtvc = midf[midf.GID == 'Gene_8234_9572'].Type.value_counts()
mtvc

neutral_gene_ratios['FBgn0013680']
(mtvc['missense_variant']/mtvc['synonymous_variant'])
(mtvc['missense_variant']/mtvc['synonymous_variant'])/gnr['FBgn0013680']

(sum(midf[(midf.GID == 'Gene_8234_9572') & (midf.Type == 'missense_variant')].Pi)/sum(midf[(midf.GID == 'Gene_8234_9572') & (midf.Type == 'synonymous_variant')].Pi)) / gnr['FBgn0013680']
