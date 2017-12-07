
# coding: utf-8

# In[ ]:


import Bio
from Bio import SeqIO
from Bio import GenBank
import sys
import os
import argparse
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pandas as pd
import re


# In[ ]:


user_input = input('Please enter your Genbank file path here and hit ENTER:')


# In[ ]:


#"for" loop so that all the genes are extracted. Their locations are taken and appended into a string
Coor = []
for rec in SeqIO.parse((user_input), "genbank"):
    if rec.features:
        for feature in rec.features:
            if feature.type == "gene":
                Coor.append(str(feature.location))


# In[ ]:


#The strings are joined together so that regex can be used. Here regex is used to slice out only the digits.
str1 = ''.join(Coor)
x = re.findall(r'\d+', str1)
x2 = list(map(int, x)) #The string is converted to a list with integers for later use in pandas DF.


# In[ ]:


#Parsing out the type"source" where the feature locations of the genome start and end can be found. 
CoorSE = []
for rec in SeqIO.parse((user_input), "genbank"):
    if rec.features:
        for feature in rec.features:
            if feature.type == "source":
                CoorSE.append(str(feature.location))


# In[ ]:


str2 = ''.join(CoorSE)#Turned into a string for regex.


# In[ ]:


#Regex is used to cut only the digits out from the start/end coordinates of the genome.
xx = re.findall(r'\d+', str2)#Coordinates are grabbed here
yy = re.findall(r'\d+', str2)[:1]#The start coordinate(0) is grabbed here
zz = re.findall(r'\d+', str2)[1:2]#The end coordiante is grabbed here
yyL = list(map(int, yy))#They are then converted into lists for the pandas data frame later.
zzL = list(map(int, zz))


# In[ ]:


starts = x2[::2] #All the start coordinates are grabbed from the original list of all coordinates and then seperated into two seperate list for "starts and stops".
stops = x2[1::2]


# In[ ]:



IGR_ends = [x-1 for x in starts]#First the IGR ends are calculated and -1 from the (starts) of genes
IGR_realends = IGR_ends + zzL#Then the last genome coordinate is added because if it was added before, the data frame would start at -1
IGR_starts = [x+1 for x in stops]#Starts of IGRs were calculated by adding +1 to the stops of genes.
IGR_realstarts = yyL + IGR_starts


# In[ ]:


#First pandas dataframe is made with two columns (start and stop) based on the IGR starts and stops variables created before.
df = pd.DataFrame({'start': IGR_realstarts})
df['stop'] = IGR_realends


# In[ ]:


#Boolean loop to intepret whether there is an intergenic region present.
Loopy = []
for start, stop in zip(list(map(int, df.start)), list(map(int, df.stop))):
    if start < stop:#If the start of an IGR is less than the end, this means there is an IGR present.
        Loopy.append(True)
    else:
        Loopy.append(False) #if not, it is not appended into the data frame.


# In[ ]:


IGRs = pd.Series(Loopy)#Series was created for the boolean loop so that data can rendered based on the previous df made.
df[Loopy]
ResultDF = df[Loopy]#Final result


# In[ ]:


IGRs = [] #Assigning one new column to the dataframe with a name being inserted incrementally.
num_IGR_len = len(ResultDF) +1
for x in range(1,num_IGR_len):
    IGRs.append('IGR_' + str(x)) 
idx = 0
ResultDF.insert(loc=idx, column='IGRs', value=IGRs)


# In[ ]:


CSVname = input('make a name for your file and add the extension .csv')
ResultDF.to.csv(CSVname) #Allows user to make a file name and export the df into a .csv file

