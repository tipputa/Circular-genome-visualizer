# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30, 2019

@author: tipputa
"""

import numpy as np
from Bio import SeqIO #, SeqFeature
import os, math, json #, re
import pandas as pd
import visualizeUtil as visualizer
import timeRecord as recorder
import collections as cl
import createCircos

class Runs():
    def __init__(self, RootDir, gbDir):
        self.RootDir = RootDir
        self.gbDir = gbDir
        self.dataDir = RootDir + "data/"
        self.CircosDir = RootDir + "circos/" 
        self.blastDBDir = self.dataDir + "BlastDB/"
        self.summaryDir = self.dataDir + "summary/"
        self.blastResultDir = self.dataDir + "BlastResult/"    
        self.proteinFasDir = self.dataDir + "ProteinFasta/"    
        
    def writeJson(self, jsonRes):
        fname = self.RootDir + "text.json"
        with open(fname, "w") as f:
            json.dump(jsonRes, f, indent=4)
            
    def setJson(self, df, consensus_gsize, circosIN, rotated_circosIN,masterDict):
        res = cl.OrderedDict()
        res["consensus_genome_size"] = int(consensus_gsize)
        res["ortholog_info"] = []
        res["each_genome_info"] = []
        for i in range(0,len(df.index)):
            tmp = cl.OrderedDict()
            genes = []
            cIn = circosIN[df.index[i]]
            cIn_rotated = rotated_circosIN[i]
            tmp["genome_name"] = df.iloc[i,0]
            tmp["acc"] = df.iloc[i,1]
            tmp["genome_size"] = int(df.iloc[i,2])
            tmp["distance_from_consensus"] = int(df.iloc[i,3])
            tmp["order"] = i
            for j in range(0, len(cIn_rotated.index)):
                gene = cl.OrderedDict()
                gene["start"] = int(cIn.iloc[j, 1])
                gene["end"] = int(cIn.iloc[j, 2])
                gene["start_rotated"] = int(cIn_rotated.iloc[j, 1])
                gene["end_rotated"] = int(cIn_rotated.iloc[j, 2])
                gene["strand"] = int(cIn_rotated.iloc[j, 3])
                gene["locusTag"] = cIn_rotated.iloc[j, 4]
                gene["consensusId"] = int(Runs.convertTag2Color(cIn_rotated.iloc[j,4], masterDict))
                gene["angle"] = int(Runs.convertTag2Color(cIn_rotated.iloc[j,4], masterDict))
                
                genes.append(gene)
            tmp["genes"] = genes                
            res["each_genome_info"].append(tmp)                
        self.writeJson(res)
        
        
    def readGB(self, file, genome_size, CircosIN, df): ## gbファイルの読み込み、locusTagから角度への置換)
        record = SeqIO.read(file, "genbank")
        gsize = len(record.seq)
        step = gsize / 360
        genome_size.append(gsize)
        target = list(df[str(record.id)])
        cdsL = []
        newList = []
        locusDict = {"-" : "-"}
    
        for feature in record.features:
            if feature.type == "CDS":
                locus_tag = feature.qualifiers.get("locus_tag")
                pos_start = feature.location.start
                pos_end = feature.location.end
                strand = feature.strand
                tmp_out = pd.Series(["chr1", pos_start, pos_end, strand, locus_tag[0]])
                if tmp_out[2] - tmp_out[1] > 20000:
                    print("The position of first gene is modified")
                    tmp_out[2] = tmp_out[1] + 100
                    
                cdsL.append(tmp_out)
                cds_angle = round (feature.location.start / step)
                locusDict.update({locus_tag[0]:cds_angle})
    
        while len(target) != 0:
            tmp = target.pop()
            newList.append(locusDict[tmp])
    
        newList.reverse()
        df[str(record.id)] = pd.Series(newList)
        output = pd.DataFrame(cdsL)
        CircosIN.append(output)
        output.to_csv(self.RootDir + "/circos/data/" + str(record.id) + ".original.txt", sep = "\t",header=None, index=None)
        return(df)    
        
    def makeConsensusDict(df):
        rows, cols = df.shape
        colname = df.columns
        Dict = {}
        
        for col in colname[0:cols-3]:
            tmp = df.ix[df[col] != "-",[col,"Consensus"]]
            tmp.index = tmp[col]
            tmp=tmp[~tmp.duplicated([col])]  # Duplicated values are removed
            Dict.update(tmp["Consensus"])
        return(Dict)

    def convertTag2Color(tag, Dict):
            if tag in Dict:
                return int(Dict[tag])
            else:
                return -1

    def createCircosInput(RootDir, df_tmp, consensus_gsize, CircosIN, df_info, isRotate = True):
        temp_dict = Runs.makeConsensusDict(df_tmp)
        nrow, ncol = df_info.shape
        circosIn2 = []
        for Loop in range(0,nrow):
            tmp = CircosIN[df_info.index[Loop]]
            info = df_info.iloc[Loop]
            tmp = visualizer.rotatePosition(tmp, info, consensus_gsize, isRotate)
            circosIn2.append(tmp)
        return circosIn2, temp_dict

def runs(RootDir, gbDir):
    run = Runs(RootDir, gbDir)
    df = pd.read_csv(RootDir + '/data/all_blast_results.tsv', delimiter = "\t")
    df_num = df.copy()
    files = os.listdir(gbDir)
    
    genome_size = []
    CircosIN = []
    LoopCounter = 1
    print ("Start creating JSON file")
    NumGB = len([file for file in files if file.endswith(".gb")  or file.endswith('.gbk')])
    print ("Number of GenBank files:", NumGB)
    for file in files:
        if file.endswith(".gb") or file.endswith('.gbk'): 
            print("Reading: " + file + "  "+ str(LoopCounter) + "/" + str(NumGB)),
            LoopCounter += 1
            df_num = run.readGB(gbDir + file, genome_size, CircosIN, df_num)

    consensus_gsize = round(np.mean(genome_size))
    df_locusTag_aligned = pd.read_csv(RootDir + '\\data\\locusTag_aligned_angle.tsv', delimiter = "\t")
    df_info = pd.read_csv(RootDir + '\\data\\output_information.tsv', delimiter = "\t")
    df_info2 = pd.concat([pd.Series(files), df_info], axis=1)    
    sort_key = "Deviation (Aligned)"
    df_info_sorted = df_info2.sort_values(by=[sort_key], ascending=False)
    rotated_circosIn, master_dict = Runs.createCircosInput(RootDir, df_locusTag_aligned, consensus_gsize, CircosIN, df_info_sorted)
    run.setJson(df_info_sorted, consensus_gsize, CircosIN, rotated_circosIn, master_dict)
    
