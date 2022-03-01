import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob 
import sys
import math

LowDownsampledSomaticbardf = {k:[] for k in ['Stage','Ind','NormCount']}
for stage, stdf in LowDownsampledSomaticMitodf[LowDownsampledSomaticMitodf.IsSomatic].groupby("Stage"):
    indvc = stdf.Source.value_counts()
    for ind in indvc.index:
        LowDownsampledSomaticbardf['Stage'].append(stage)
        LowDownsampledSomaticbardf['Ind'].append(ind)
        #needs to be changed
        #need to divide Norm Count by number of source stages 
        LowDownsampledSomaticbardf['NormCount'].append(indvc[ind] / np.sqrt(New_Mean_DepthValues[ind])) 
LowDownsampledSomaticbardf = pd.DataFrame(LowDownsampledSomaticbardf)
#LowDownsampledSomaticbardf.sort_values(by=['Stage'], inplace=True, ascending=False)

LowDownsampledSomaticbardf = LowDownsampledSomaticbardf.reindex([2,3,4,5,6,7,8,0,1])
LowDownsampledSomaticbardf

sns.set_theme(style="whitegrid")

ax = sns.barplot(x="Stage", y = "NormCount", data = LowDownsampledSomaticbardf[LowDownsampledSomaticbardf.Ind != "sf1"],  ci = None, palette="Blues")
ax.set_xlabel("Life Stage")
ax.set_ylabel("Normalized Mutation Count")
ax.set_title('Somatic Mutation Based on Lifestage')

