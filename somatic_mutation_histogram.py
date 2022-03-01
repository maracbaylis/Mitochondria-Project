import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob 
import sys
import math

sns.set(font_scale = .8)
  
y,x=np.histogram((DownsampledSomaticMitodf.SampleFreq),bins=10,density=True)    
x = [round(10**(i),4) for i in x]
ax=sns.barplot(x=x[1:],y=y,color="seagreen")

ax.set_title("Binned Log Circle Frequency Spectrum")
ax.set_ylabel("Density")
ax.set_xlabel("Maximum Sample Frequency")

plt.savefig("_cfs.png", figsize = (10,28),dpi=800)
