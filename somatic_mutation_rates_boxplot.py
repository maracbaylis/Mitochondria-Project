import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob 
import sys
import math

mutdf = DownsampledSomaticMitodf

#SSN = Source or maybe ID
#REF = RefrenceBase
#ALT = JustMut

#mutdf['Type'] = mutdf.Ref + ">" + mutdf.Alt
#print some basic statistics.
ssnvc = mutdf.Source.value_counts() #SINGLE SAMPLE IDENTIFIER 
for ssn in ssnvc.index:
    print("Individual {} has {} total mutations.".format(ssn, ssnvc[ssn]))
    mtvc = mutdf[mutdf.Source == ssn].Source.value_counts()
    for t in mtvc.index:
        print("with {} {} mutations".format(mtvc[t],t))

print("Generating graph of mutation rates.")
#need to create a diconary of basecount values  from dataframe given
#key = base 
#value = number of times you see that base 

alist = mutdf.RefrenceBase.tolist()
print(alist)
def CountFrequency(my_list):
    bcd = {}
    for item in my_list:
        if (item in bcd):
            bcd[item] += 1
        else:
            bcd[item] = 1 
    for key, value in bcd.items():
        print ((key, value))
    return bcd
        
#if __name__ == "__main__":
   # my_list == alist    
    #CountFrequency(my_list)
    
bcd = CountFrequency(alist)
pairs = []
for a in 'ACGT':
    for b in 'ACGT':
        if a != b:
            pairs.append((a,b))
            
cdf = {k:[] for k in ['ind', 'Mutation Type', 'Rate (Detections per Site per Depth)']}
for ssn in mutdf.Source.value_counts().index:
    for ref, alt in pairs:
        subdf = mutdf[(mutdf.SampleFreq < .25) & (mutdf.Source == ssn) & (mutdf.RefrenceBase == ref) & (mutdf.JustMut == alt)]
        
        #for i,v in subdf.iterrows():
            #print(f"{i=},{v=}")
            #print(f"{v.NumOfReads=}")
            
        count = sum([(math.floor((int(v.NumOfReads) * v.SampleFreq))) for i,v in subdf.iterrows()])
        
        cdf['ind'].append(ssn)
        cdf['Mutation Type'].append(ref + '>' + alt)

        if ssn in bcd.keys():
            c = bcd[ssn][ref]
        else:
            c = bcd[ref]
        
        cdf['Rate (Detections per Site per Depth)'].append(count/c)
    
        
cdf = pd.DataFrame(cdf)
ax = sns.boxplot(x = 'Mutation Type', y = 'Rate (Detections per Site per Depth)', data = cdf, color = "seagreen") #pretty color (color = ?)
ax.set_yticklabels(['{:.2e}'.format(v) for v in list(ax.get_yticks())])
ax.set_xlabel("Mutation Type")
ax.set_ylabel("Rate (Detections per Site per Depth)")
ax.set_title('Somatic Mutation Rates')
plt.savefig("thegraph.png",dpi=800)
