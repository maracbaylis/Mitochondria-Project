import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob 
import sys
import math


#Downsample the Somatic Mutation Dataframe to normalize for depth
#Determine which depth value to downsample to 
total = Somaticmitodf.loc[Somaticmitodf['Source'] == 'rf1', 'NumOfReads'].astype(int).sum()
print("rf1:", total)
total = Somaticmitodf.loc[Somaticmitodf['Source'] == 'rf2', 'NumOfReads'].astype(int).sum()
print("rf2:",total)
total = Somaticmitodf.loc[Somaticmitodf['Source'] == 'rp1', 'NumOfReads'].astype(int).sum()
print("rp1:", total)
total = Somaticmitodf.loc[Somaticmitodf['Source'] == 'rp2', 'NumOfReads'].astype(int).sum()
print("rp2:", total)
total = Somaticmitodf.loc[Somaticmitodf['Source'] == 'sl2', 'NumOfReads'].astype(int).sum()
print("sl2:", total)
total = Somaticmitodf.loc[Somaticmitodf['Source'] == 'rl2', 'NumOfReads'].astype(int).sum()
print("rl2:", total)
total = Somaticmitodf.loc[Somaticmitodf['Source'] == 'rl1', 'NumOfReads'].astype(int).sum()
print("rl1:", total)
total = Somaticmitodf.loc[Somaticmitodf['Source'] == 'sl1', 'NumOfReads'].astype(int).sum()
print("sl1:", total)
total = Somaticmitodf.loc[Somaticmitodf['Source'] == 'sp1', 'NumOfReads'].astype(int).sum()
print("sp1:", total)

#Remove individuals from the Somaticmitodf dataframe that are dramatic outliers (very low depth)
index_names = Somaticmitodf[Somaticmitodf['Source'] == 'sf1' ].index
Somaticmitodf.drop(index_names, inplace = True)
Somaticmitodf
index_names = Somaticmitodf[Somaticmitodf['Source'] == 'sp2' ].index
Somaticmitodf.drop(index_names, inplace = True)
Somaticmitodf


#Downsample the Somatic Mutation Dataframe to normalize for depth 
def summingPreviousLengths(old_array):
    """
    takes a list, makes it into a list where each term is the sum of all previous terms
    """
    new_array = [old_array[0]]
    
    for i in old_array[1:]:
        new_array.append(new_array[-1] + i)
    return new_array

test_array = [2, 2, 2, 2]

my_cool_new_list = summingPreviousLengths(test_array)
print(my_cool_new_list)

rf2df = Somaticmitodf[Somaticmitodf.Source == "rf2"]
rf2df["NumOfReads"] = pd.to_numeric(rf2df["NumOfReads"])

def downsample(pileup_dataframe = rf2df, target=55649):
    """
    input is a pileupdataframe with columns named NumOfReads and ReadBases
    """
    sum_of_lens = sum(pileup_dataframe.NumOfReads)
    rand_indices = sorted(np.random.choice(range(0, sum_of_lens), size=(sum_of_lens - target), replace=False))[::-1]
    sum_of_read_lengths = [0]+summingPreviousLengths(list(pileup_dataframe.NumOfReads))

    new_column_values = [c for c in pileup_dataframe["ReadBases"]]
    for index in rand_indices:
        for i in range(len(sum_of_read_lengths)):
            if index < sum_of_read_lengths[i]:
                row_to_drop = i
                break
                
        term_to_drop = index - sum_of_read_lengths[i-1]
        
        rowstring = new_column_values[row_to_drop-1]
        new_rowstring_list = list(rowstring)
        new_rowstring_list[term_to_drop] = "_"
        new_rowstring = "".join(new_rowstring_list)
        new_column_values[row_to_drop-1] = new_rowstring
    pileup_dataframe["DownReadBases"] = new_column_values
    pileup_dataframe["DownReadBases"] = pileup_dataframe["DownReadBases"].apply(lambda x:x.replace("_",""))
    pileup_dataframe["DownNumOfReads"] = pileup_dataframe["DownReadBases"].apply(lambda x:len(x))
    return pileup_dataframe

data_set = rf2df
downsample(data_set, target = 55649)


rf1df = Somaticmitodf[Somaticmitodf.Source == "rf1"]
def downsample(pileup_dataframe = rf1df, target=55649):
    """
    input is a pileupdataframe with columns named NumOfReads and ReadBases
    """
    sum_of_lens = (sum(pileup_dataframe.NumOfReads.astype(int)))
    rand_indices = sorted(np.random.choice(range(0, sum_of_lens), size=(sum_of_lens - target), replace=False))[::-1]
    sum_of_read_lengths = [0]+summingPreviousLengths(list(pileup_dataframe.NumOfReads.astype(int)))
 
    new_column_values = [c for c in pileup_dataframe["ReadBases"]]
    for index in rand_indices:
        for i in range(len(sum_of_read_lengths)):
            if index < sum_of_read_lengths[i]:
                row_to_drop = i
                break
                
        term_to_drop = index - sum_of_read_lengths[i-1]
        
        rowstring = new_column_values[row_to_drop-1]
        new_rowstring_list = list(rowstring)
        new_rowstring_list[term_to_drop] = "_"
        new_rowstring = "".join(new_rowstring_list)
        new_column_values[row_to_drop-1] = new_rowstring
    pileup_dataframe["DownReadBases"] = new_column_values
    pileup_dataframe["DownReadBases"] = pileup_dataframe["DownReadBases"].apply(lambda x:x.replace("_",""))
    pileup_dataframe["DownNumOfReads"] = pileup_dataframe["DownReadBases"].apply(lambda x:len(x))
    return pileup_dataframe

data_set = rf1df
downsample(data_set, target = 55649)

rp1df = Somaticmitodf[Somaticmitodf.Source == "rp1"]
def downsample(pileup_dataframe = rp1df, target=55649):
    """
    input is a pileupdataframe with columns named NumOfReads and ReadBases
    """
    sum_of_lens = (sum(pileup_dataframe.NumOfReads.astype(int)))
    rand_indices = sorted(np.random.choice(range(0, sum_of_lens), size=(sum_of_lens - target), replace=False))[::-1]
    sum_of_read_lengths = [0]+summingPreviousLengths(list(pileup_dataframe.NumOfReads.astype(int)))
 
    new_column_values = [c for c in pileup_dataframe["ReadBases"]]
    for index in rand_indices:
        for i in range(len(sum_of_read_lengths)):
            if index < sum_of_read_lengths[i]:
                row_to_drop = i
                break
                
        term_to_drop = index - sum_of_read_lengths[i-1]
        
        rowstring = new_column_values[row_to_drop-1]
        new_rowstring_list = list(rowstring)
        new_rowstring_list[term_to_drop] = "_"
        new_rowstring = "".join(new_rowstring_list)
        new_column_values[row_to_drop-1] = new_rowstring
    pileup_dataframe["DownReadBases"] = new_column_values
    pileup_dataframe["DownReadBases"] = pileup_dataframe["DownReadBases"].apply(lambda x:x.replace("_",""))
    pileup_dataframe["DownNumOfReads"] = pileup_dataframe["DownReadBases"].apply(lambda x:len(x))
    return pileup_dataframe

data_set = rp1df
downsample(data_set, target = 55649)

rp2df = Somaticmitodf[Somaticmitodf.Source == "rp2"]
def downsample(pileup_dataframe = rp2df, target=55649):
    """
    input is a pileupdataframe with columns named NumOfReads and ReadBases
    """
    sum_of_lens = (sum(pileup_dataframe.NumOfReads.astype(int)))
    rand_indices = sorted(np.random.choice(range(0, sum_of_lens), size=(sum_of_lens - target), replace=False))[::-1]
    sum_of_read_lengths = [0]+summingPreviousLengths(list(pileup_dataframe.NumOfReads.astype(int)))
 
    new_column_values = [c for c in pileup_dataframe["ReadBases"]]
    for index in rand_indices:
        for i in range(len(sum_of_read_lengths)):
            if index < sum_of_read_lengths[i]:
                row_to_drop = i
                break
                
        term_to_drop = index - sum_of_read_lengths[i-1]
        
        rowstring = new_column_values[row_to_drop-1]
        new_rowstring_list = list(rowstring)
        new_rowstring_list[term_to_drop] = "_"
        new_rowstring = "".join(new_rowstring_list)
        new_column_values[row_to_drop-1] = new_rowstring
    pileup_dataframe["DownReadBases"] = new_column_values
    pileup_dataframe["DownReadBases"] = pileup_dataframe["DownReadBases"].apply(lambda x:x.replace("_",""))
    pileup_dataframe["DownNumOfReads"] = pileup_dataframe["DownReadBases"].apply(lambda x:len(x))
    return pileup_dataframe

data_set = rp2df
downsample(data_set, target = 55649)

sl2df = Somaticmitodf[Somaticmitodf.Source == "sl2"]
def downsample(pileup_dataframe = sl2df, target=55649):
    """
    input is a pileupdataframe with columns named NumOfReads and ReadBases
    """
    sum_of_lens = (sum(pileup_dataframe.NumOfReads.astype(int)))
    rand_indices = sorted(np.random.choice(range(0, sum_of_lens), size=(sum_of_lens - target), replace=False))[::-1]
    sum_of_read_lengths = [0]+summingPreviousLengths(list(pileup_dataframe.NumOfReads.astype(int)))
 
    new_column_values = [c for c in pileup_dataframe["ReadBases"]]
    for index in rand_indices:
        for i in range(len(sum_of_read_lengths)):
            if index < sum_of_read_lengths[i]:
                row_to_drop = i
                break
                
        term_to_drop = index - sum_of_read_lengths[i-1]
        
        rowstring = new_column_values[row_to_drop-1]
        new_rowstring_list = list(rowstring)
        new_rowstring_list[term_to_drop] = "_"
        new_rowstring = "".join(new_rowstring_list)
        new_column_values[row_to_drop-1] = new_rowstring
    pileup_dataframe["DownReadBases"] = new_column_values
    pileup_dataframe["DownReadBases"] = pileup_dataframe["DownReadBases"].apply(lambda x:x.replace("_",""))
    pileup_dataframe["DownNumOfReads"] = pileup_dataframe["DownReadBases"].apply(lambda x:len(x))
    return pileup_dataframe

data_set = sl2df
downsample(data_set, target = 55649)

rl2df = Somaticmitodf[Somaticmitodf.Source == "rl2"]
def downsample(pileup_dataframe = rl2df, target=55649):
    """
    input is a pileupdataframe with columns named NumOfReads and ReadBases
    """
    sum_of_lens = (sum(pileup_dataframe.NumOfReads.astype(int)))
    rand_indices = sorted(np.random.choice(range(0, sum_of_lens), size=(sum_of_lens - target), replace=False))[::-1]
    sum_of_read_lengths = [0]+summingPreviousLengths(list(pileup_dataframe.NumOfReads.astype(int)))
 
    new_column_values = [c for c in pileup_dataframe["ReadBases"]]
    for index in rand_indices:
        for i in range(len(sum_of_read_lengths)):
            if index < sum_of_read_lengths[i]:
                row_to_drop = i
                break
                
        term_to_drop = index - sum_of_read_lengths[i-1]
        
        rowstring = new_column_values[row_to_drop-1]
        new_rowstring_list = list(rowstring)
        new_rowstring_list[term_to_drop] = "_"
        new_rowstring = "".join(new_rowstring_list)
        new_column_values[row_to_drop-1] = new_rowstring
    pileup_dataframe["DownReadBases"] = new_column_values
    pileup_dataframe["DownReadBases"] = pileup_dataframe["DownReadBases"].apply(lambda x:x.replace("_",""))
    pileup_dataframe["DownNumOfReads"] = pileup_dataframe["DownReadBases"].apply(lambda x:len(x))
    return pileup_dataframe

data_set = rl2df
downsample(data_set, target = 55649)

rl1df = Somaticmitodf[Somaticmitodf.Source == "rl1"]
def downsample(pileup_dataframe = rl1df, target=55649):
    """
    input is a pileupdataframe with columns named NumOfReads and ReadBases
    """
    sum_of_lens = (sum(pileup_dataframe.NumOfReads.astype(int)))
    rand_indices = sorted(np.random.choice(range(0, sum_of_lens), size=(sum_of_lens - target), replace=False))[::-1]
    sum_of_read_lengths = [0]+summingPreviousLengths(list(pileup_dataframe.NumOfReads.astype(int)))
 
    new_column_values = [c for c in pileup_dataframe["ReadBases"]]
    for index in rand_indices:
        for i in range(len(sum_of_read_lengths)):
            if index < sum_of_read_lengths[i]:
                row_to_drop = i
                break
                
        term_to_drop = index - sum_of_read_lengths[i-1]
        
        rowstring = new_column_values[row_to_drop-1]
        new_rowstring_list = list(rowstring)
        new_rowstring_list[term_to_drop] = "_"
        new_rowstring = "".join(new_rowstring_list)
        new_column_values[row_to_drop-1] = new_rowstring
    pileup_dataframe["DownReadBases"] = new_column_values
    pileup_dataframe["DownReadBases"] = pileup_dataframe["DownReadBases"].apply(lambda x:x.replace("_",""))
    pileup_dataframe["DownNumOfReads"] = pileup_dataframe["DownReadBases"].apply(lambda x:len(x))
    return pileup_dataframe

data_set = rl1df
downsample(data_set, target = 55649)

sl1df = Somaticmitodf[Somaticmitodf.Source == "sl1"]
def downsample(pileup_dataframe = sl1df, target=55649):
    """
    input is a pileupdataframe with columns named NumOfReads and ReadBases
    """
    sum_of_lens = (sum(pileup_dataframe.NumOfReads.astype(int)))
    rand_indices = sorted(np.random.choice(range(0, sum_of_lens), size=(sum_of_lens - target), replace=False))[::-1]
    sum_of_read_lengths = [0]+summingPreviousLengths(list(pileup_dataframe.NumOfReads.astype(int)))
 
    new_column_values = [c for c in pileup_dataframe["ReadBases"]]
    for index in rand_indices:
        for i in range(len(sum_of_read_lengths)):
            if index < sum_of_read_lengths[i]:
                row_to_drop = i
                break
                
        term_to_drop = index - sum_of_read_lengths[i-1]
        
        rowstring = new_column_values[row_to_drop-1]
        new_rowstring_list = list(rowstring)
        new_rowstring_list[term_to_drop] = "_"
        new_rowstring = "".join(new_rowstring_list)
        new_column_values[row_to_drop-1] = new_rowstring
    pileup_dataframe["DownReadBases"] = new_column_values
    pileup_dataframe["DownReadBases"] = pileup_dataframe["DownReadBases"].apply(lambda x:x.replace("_",""))
    pileup_dataframe["DownNumOfReads"] = pileup_dataframe["DownReadBases"].apply(lambda x:len(x))
    return pileup_dataframe

data_set = sl1df
downsample(data_set, target = 55649)

sp1df = Somaticmitodf[Somaticmitodf.Source == "sp1"]
def downsample(pileup_dataframe = sp1df, target=55649):
    """
    input is a pileupdataframe with columns named NumOfReads and ReadBases
    """
    sum_of_lens = (sum(pileup_dataframe.NumOfReads.astype(int)))
    rand_indices = sorted(np.random.choice(range(0, sum_of_lens), size=(sum_of_lens - target), replace=False))[::-1]
    sum_of_read_lengths = [0]+summingPreviousLengths(list(pileup_dataframe.NumOfReads.astype(int)))
 
    new_column_values = [c for c in pileup_dataframe["ReadBases"]]
    for index in rand_indices:
        for i in range(len(sum_of_read_lengths)):
            if index < sum_of_read_lengths[i]:
                row_to_drop = i
                break
                
        term_to_drop = index - sum_of_read_lengths[i-1]
        
        rowstring = new_column_values[row_to_drop-1]
        new_rowstring_list = list(rowstring)
        new_rowstring_list[term_to_drop] = "_"
        new_rowstring = "".join(new_rowstring_list)
        new_column_values[row_to_drop-1] = new_rowstring
    pileup_dataframe["DownReadBases"] = new_column_values
    pileup_dataframe["DownReadBases"] = pileup_dataframe["DownReadBases"].apply(lambda x:x.replace("_",""))
    pileup_dataframe["DownNumOfReads"] = pileup_dataframe["DownReadBases"].apply(lambda x:len(x))
    return pileup_dataframe

data_set = sp1df
downsample(data_set, target = 55649)

alldfs = [rf1df, rf2df, rp1df, rp2df, sl2df, rl2df, rl1df, sl1df, sp1df]
DownsampledSomaticMitodf = pd.concat(alldfs)
DownsampledSomaticMitodf
