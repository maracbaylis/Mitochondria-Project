import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob 
import sys
import math


pred = {}
colnames = ["ID", "BasedPosition", "RefrenceBase", "NumOfReads", "ReadBases",'BaseQualities','Source']
for c in colnames:
    pred[c] = []
path = "../mbaylis/mitochondion_pileups_/mitochondion_pileups_variants/" 
all_files = glob.glob(path + "/*.pileup")
print(all_files)
for file in all_files:
    with open(file) as inf:
        for i,entry in enumerate(inf.readlines()):
            if i >= 0:
                row = entry.strip().split("\t")
                pred['Source'].append(file)
                '''
                for value in pred["SomaticMutation"]:
                    pred['ReadBases'].
                        print("False")
                else
                    print("True")
                '''
                for i, value in enumerate(row):
                    name = colnames[i]
                    pred[name].append(value)
mitodf = pd.DataFrame(pred)

useless_wrds = ["../mbaylis/mitochondion_pileups_/mitochondion_pileups_variants/", "mitochondion pileups", "mitochondion", "pileup", "_", "variants."]
for wrd in useless_wrds:
    mitodf['Source'] = mitodf['Source'].str.replace(wrd, ' ')

is_somatic = []
for bases in mitodf.ReadBases:
    bases = [b for b in bases if b != 'N']
    dotcount = bases.count(".")
    if dotcount > .75 * len(bases):
        is_somatic.append(True)
    else:
        is_somatic.append(False)
len(is_somatic)
mitodf['IsSomatic'] = is_somatic
mitodf[mitodf.IsSomatic]
#ADD STAGE DATA 
staged = {'f':"adult", 'l':'larva', 'p': 'pupa'}
stages = []
for source in mitodf.Source:
    #print(source)
    letter = source.strip()[1]
    #print(letter)
    if letter in staged:
        stages.append(staged[letter])
    else:
        stages.append("NA")
stages
mitodf["Stage"] = stages

mitodf["Source"] = mitodf.Source.apply(lambda x:x.strip())
mitodf[~mitodf.isin([np.nan, np.inf, -np.inf]).any(1)]
mitodf = mitodf.replace("NA", np.nan)
mitodf.dropna()


mitodf
Somaticmitodf = mitodf[mitodf.IsSomatic]
Somaticmitodf["JustMut"] = Somaticmitodf['ReadBases'].str.replace('[^\w\s]','')

Somaticmitodf = Somaticmitodf[~Somaticmitodf['JustMut'].str.contains("[N]").fillna(False)]

Somaticmitodf["NumberOfJustMut"] = Somaticmitodf['JustMut'].map(lambda calc: len(calc))
Somaticmitodf["SampleFreq"] = Somaticmitodf[('NumberOfJustMut')].astype(int)/Somaticmitodf[("NumOfReads")].astype(int)

Somaticmitodf
