# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 21:30:29 2019

@author: tipputa
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 09:26:06 2016

@author: tipputa
"""

from Bio import SeqIO 
import pandas as pd

def getPosition(file, df): # without strand
    record = SeqIO.read(file, "genbank")
    target = df[str(record.id)]        
    
    for feature in record.features:
        if feature.type == "CDS":
            locus_tag = feature.qualifiers.get("locus_tag")
            pos_start = str(feature.location.start)
            target[target==locus_tag[0]] = pos_start
    df[str(record.id)] = target
    return(df)

def getPositionWithStrand(file, df): ## including strand information
    record = SeqIO.read(file, "genbank")
    print(str(record.id))
    target = list(df[str(record.id)])
    newList = []
    locusDict = {"-" : "-"}

    for feature in record.features:
        if feature.type == "CDS":
            locus_tag = feature.qualifiers.get("locus_tag")
            pos_start = feature.location.start
            strand = feature.location.strand     
            final = pos_start * strand
            locusDict.update({locus_tag[0]:final})

    while len(target) != 0:
        tmp = target.pop()
        newList.append(locusDict[tmp])

    newList.reverse()
    df[str(record.id)] = pd.Series(newList)
    return(df)    
    