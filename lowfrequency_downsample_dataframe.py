import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob 
import sys
import math

DownsampledSomaticMitodf["DownJustMut"] = DownsampledSomaticMitodf['DownReadBases'].str.replace('[^\w\s]','')
DownsampledSomaticMitodf.dropna(subset=['DownJustMut'],inplace=True)
DownsampledSomaticMitodf["DownNumberOfJustMut"] = DownsampledSomaticMitodf['DownJustMut'].map(lambda calc: len(calc))
DownsampledSomaticMitodf["DownSampleFreq"] = DownsampledSomaticMitodf[('DownNumberOfJustMut')].astype(int)/DownsampledSomaticMitodf[("DownNumOfReads")].astype(int)

DownsampledSomaticMitodf
'''
rf2dfdf_mean = rf2df["DownNumOfReads"].mean()
print("rf2df_mean", rf2dfdf_mean)
'''
New_Mean_DepthValues = {'rf2':330.3195266272189,
                        'rf1': 296.8031914893617,
                        'rp1': 460.42148760330576,
                        'rp2': 515.7870370370371,
                        'sl2': 471.96610169491527,
                        'rl2': 562.6060606060606,
                        'rl1': 592.4787234042553,
                        'sl1': 598.752688172043,
                        'sp1': 807.1594202898551}
np.quantile(DownsampledSomaticMitodf.SampleFreq, [0.05,.25,.50,.75,.95])
LowDownsampledSomaticMitodf = DownsampledSomaticMitodf[DownsampledSomaticMitodf.SampleFreq < 0.00067069]
LowDownsampledSomaticMitodf


