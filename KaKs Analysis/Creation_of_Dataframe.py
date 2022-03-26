import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from tqdm.notebook import tqdm
import math
from scipy.stats import hypergeom
from scipy.stats import fisher_exact

#path to tsv created usinng circle seq pipeline

midf = pd.read_csv("../mbaylis/mitochondion_pileups_/fixed_pileups/fixed_all_mutations.tsv",sep='\t').drop("Unnamed: 0",axis=1)

midf['Missense'] = midf.Type.apply(lambda x:('missense' in x))
midf['Synonymous'] = midf.Type.apply(lambda x:('synonymous' in x))
midf.Missense.value_counts()

mttypes = {}
for b1 in 'ACGT':
    for b2 in 'ACGT':
        if b1 != b2:
            mttypes[(b1,b2)] = 0
for i,d in midf.iterrows():
    mttypes[(d.Ref, d.Alt)] += 1
mttypes 

norm_mttypes = {}
for k,v in mttypes.items():
    norm_mttypes[k] = v/midf.shape[0]
norm_mttypes

def parse_lookupplus(path):
    lookd = {}
    with open(path) as inf:
        for entry in inf:
            spent = entry.strip().split('\t')
            #store these results as a dictionary of dictionaries, outer key chro-loc, inner keys each of the bases, gene id, strand
            #when accessing with mutdf, grab the mutation's spot, check the reference, then check the results for the alternative
            subd = {}
            try:
                subd['ref'] = spent[2]
                subd['syn'] = spent[3].split(',')
                subd['non'] = spent[4].split(',')
                subd['sgain'] = spent[5].split(',')
                subd['sloss'] = spent[6].split(',')
                subd['gid'] = spent[7]
                subd['strand'] = spent[8]
                lookd[(spent[0], int(spent[1]))] = subd
            except:
                print("Can't parse entry")
                print(entry)
                continue
    return lookd
llp = parse_lookupplus("../mbaylis/mitochondion_pileups_/fixed_pileups/v2.lookupplus")
