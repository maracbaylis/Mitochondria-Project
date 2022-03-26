import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from tqdm.notebook import tqdm
import math
from scipy.stats import hypergeom
from scipy.stats import fisher_exact

adultmidf = midf[midf.Stage == "a"]
adultmidf

amttypes = {}
for b1 in 'ACGT':
    for b2 in 'ACGT':
        if b1 != b2:
            amttypes[(b1,b2)] = 0
for i,d in adultmidf.iterrows():
    amttypes[(d.Ref, d.Alt)] += 1
amttypes

norm_mttypes = {}
for k,v in mttypes.items():
    norm_mttypes[k] = v/adultmidf.shape[0]
norm_mttypes
amtvc = adultmidf[adultmidf.GID == 'Gene_8234_9572'].Type.value_counts()
amtvc
(amtvc['missense_variant']/amtvc['synonymous_variant'])/gnr['FBgn0013680']
(sum(adultmidf[(adultmidf.GID == 'Gene_8234_9572') & (adultmidf.Type == 'missense_variant')].Pi)/sum(adultmidf[(adultmidf.GID == 'Gene_8234_9572') & (adultmidf.Type == 'synonymous_variant')].Pi)) / gnr['FBgn0013680']

larvamidf = midf[midf.Stage == "l"]
larvamidf
lmttypes = {}
for b1 in 'ACGT':
    for b2 in 'ACGT':
        if b1 != b2:
            lmttypes[(b1,b2)] = 0
for i,d in larvamidf.iterrows():
    lmttypes[(d.Ref, d.Alt)] += 1
lmttypes
lnorm_mttypes = {}
for k,v in mttypes.items():
    lnorm_mttypes[k] = v/larvamidf.shape[0]
lnorm_mttypes
lmtvc = larvamidf[larvamidf.GID == 'Gene_8234_9572'].Type.value_counts()
lmtvc
(lmtvc['missense_variant']/lmtvc['synonymous_variant'])/gnr['FBgn0013680'])

(sum(larvamidf[(larvamidf.GID == 'Gene_8234_9572') & (larvamidf.Type == 'missense_variant')].Pi)/sum(larvamidf[(larvamidf.GID == 'Gene_8234_9572') & (larvamidf.Type == 'synonymous_variant')].Pi)) / gnr['FBgn0013680']


pupamidf = midf[midf.Stage == "p"]
pupamidf
pmttypes = {}
for b1 in 'ACGT':
    for b2 in 'ACGT':
        if b1 != b2:
            pmttypes[(b1,b2)] = 0
for i,d in pupamidf.iterrows():
    pmttypes[(d.Ref, d.Alt)] += 1
pmttypes
pnorm_mttypes = {}
for k,v in mttypes.items():
    pnorm_mttypes[k] = v/pupamidf.shape[0]
pnorm_mttypes

pmtvc = pupamidf[pupamidf.GID == 'Gene_8234_9572'].Type.value_counts()
(pmtvc['missense_variant']/pmtvc['synonymous_variant'])/gnr['FBgn0013680']

(sum(pupamidf[(pupamidf.GID == 'Gene_8234_9572') & (pupamidf.Type == 'missense_variant')].Pi)/sum(pupamidf[(pupamidf.GID == 'Gene_8234_9572') & (pupamidf.Type == 'synonymous_variant')].Pi)) / gnr['FBgn0013680']
