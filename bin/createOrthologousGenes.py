#! /bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 11:20:44 2017

@author: tipputa
"""

from Bio import SeqIO
import os, time
import pandas as pd
import blastUtil, blast2tsv, createTable2
from tqdm import tqdm
import timeRecord as recorder
import gbParser as parser

class Run():
    def __init__(self, RootDir, gbDir):
        print("\nSetting:\n  Output directory: %s,\n  GenBank directory: %s\n\n" % (RootDir, gbDir))
        self.RootDir = RootDir
        self.gbDir = gbDir
        self.dataDir = RootDir + "data/"
        self.CircosDir = RootDir + "circos/" 
        self.blastDBDir = self.dataDir + "BlastDB/"
        self.summaryDir = self.dataDir + "summary/"
        self.blastResultDir = self.dataDir + "BlastResult/"    
        self.proteinFasDir = self.dataDir + "ProteinFasta/" 

        def mkdir(dirName):
            if not os.path.exists(dirName):
                os.makedirs(dirName)
            if not os.path.isdir(dirName):
                print("Dir name \"",dirName,"\" can't be used. Please check file names.\n")
                raise IOError

        mkdir(self.dataDir)
        mkdir(self.CircosDir)
        mkdir(self.CircosDir + "data")
        mkdir(self.blastDBDir)        
        mkdir(self.blastResultDir)        
        mkdir(self.proteinFasDir)
        mkdir(self.summaryDir)
         
    def makeBlastDB(self,gbDir): # for runs
        def makeblastdb(record, id="."):
            proteinFileName = self.proteinFasDir + id + ".fa"
            summaryName = self.summaryDir + id + ".tsv"
            parser.createProteinFasta(record, proteinFileName)
            parser.createTSV(record, summaryName)
        
            blastUtil.createDB(
                proteinFileName, self.blastDBDir + id + ".fa", self.blastDBDir + id + ".txt" ,dbType="prot")

        fileNames = []
        acc = []
        LoopCounter = 1
        files = os.listdir(gbDir)
 
        NumGB = len([file for file in files if file.endswith(".gb")  or file.endswith('.gbk')])
        print ("Number of GenBank files:", NumGB)
        for file in files:
            if file.endswith(".gb") or file.endswith('.gbk'): 
                record = SeqIO.read(gbDir + file, "genbank")
                print("Reading: " + file + " (" + record.id + ")  "+ str(LoopCounter) + "/" + str(NumGB)),
                LoopCounter += 1
                s = gbDir + file + "\t" + record.id + "\n"
                fileNames.append(s)            
                acc.append(record.id)
                makeblastdb(record, record.id)
        return(acc, fileNames)
 

    def makeAccList(self,gbDir): # for runsWOblast
        acc = []
        LoopCounter = 1    
        files = os.listdir(gbDir)

        NumGB = len([file for file in files if file.endswith(".gb")  or file.endswith('.gbk')])
        print ("Number of GenBank files:", NumGB)
        for file in files:
            if file.endswith(".gb") or file.endswith('.gbk'): 
                record = SeqIO.read(gbDir + file, "genbank")
                print("Reading: " + file + " (" + record.id + ")  "+ str(LoopCounter) + "/" + str(NumGB)),
                LoopCounter += 1
                acc.append(record.id)
        return(acc)
      
    
    def blastall(self, list):
        max = len(list) ** 2 
        print("Excute Protein Blast (" + str(max) + " patterns)")

        time.sleep(1)
        pbar = tqdm(total=max)
        for id in list:
            for id2 in list:
                queryFileName = self.proteinFasDir + id + ".fa"
                DBFileName = self.blastDBDir + id2 + ".fa"
                resultFileName = self.blastResultDir + id + "vs" + id2 + ".xml"
                stdout, stderr = blastUtil.blastp(queryFileName, DBFileName, resultFileName)
                pbar.update(1)
        pbar.close()
    
        for id in list:
            for id2 in list:
                forwardResult = self.blastResultDir + id + "vs" + id2 + ".xml"            
                reverseResult = self.blastResultDir + id2 + "vs" + id + ".xml"            
                targetSummary = self.summaryDir + id2 + ".tsv"
                outputFile = self.blastResultDir + id + "vs" + id2 + ".tsv"
                blast2tsv.main(forwardResult, reverseResult, targetSummary, outputFile)
    
    def afterBlast(self, list):
        for id in list:
            for id2 in list:
                forwardResult = self.blastResultDir + id + "vs" + id2 + ".xml"
                reverseResult = self.blastResultDir + id2 + "vs" + id + ".xml"
                targetSummary = self.summaryDir + id2 + ".tsv"
                outputFile = self.blastResultDir + id + "vs" + id2 + ".tsv"
                blast2tsv.main(forwardResult, reverseResult, targetSummary, outputFile)
    
    def cTable(self, pivot, targetList):
        D = createTable2.createTable(self.blastResultDir,pivot, targetList)
        DL = D.applymap(lambda x: x.locusTag)
        DL.to_csv(self.blastResultDir + pivot + "_locusTag.tsv", sep="\t", na_rep='-',index= None)
        return(DL)
        
    def makeTable(self, acc):
        a_str = map(str, acc)
        accs = ",".join(a_str)
        dfs = []
        for id in acc:
            df = self.cTable(id,accs)
            dfs.append(df)
        out = pd.DataFrame()        
        for id in dfs:
            out = pd.concat([out,id])        
        output = out.drop_duplicates()
        output2 = self.removeHyphen(output)
        return(output2)
    
    def removeHyphen(self, df):
        nrow, ncol = df.shape
        tag = []
        sum = 0
        for i in range(0, nrow):
            for j in range(0,ncol):
                if not df.ix[i,j] == "-":
                    break
                else:
                    sum = sum + 1
            if(sum == ncol):
                tag.append(df.ix[i].name)
            sum=0    
        df = df.drop(tag)
        return(df)



    def clustering(self, df):
        nrow, ncol = df.shape
        duprows = set()
        for col in range(0,ncol):
            df_tmp = df.iloc[:,col]
            df_tmp2 = df_tmp[~df_tmp.isin(["-"])]
            dup = set(df_tmp2[df_tmp2.duplicated(keep=False)].index)
            duprows = duprows.union(dup)
        allrows = set(range(0,nrow))
        targetRows = allrows.difference(duprows)
        df_res=df.iloc[list(targetRows),:]
        return df_res
    
    
    def convertLocusTag2StartPosition(self, gbDir, df):
        import getPositionFromGB as converter
        
        files = os.listdir(gbDir)
    
        for file in files:
            if file.endswith(".gb") or file.endswith('.gbk'): 
                converter.getPosition(gbDir + file, df)
        return(df)
    
    

 
def runs(timeRecorder, RootDir, gbDir):
    rec = recorder.timeRecord()        
    run = Run(RootDir, gbDir)
    print("1. Creating BLASTDB and proteinFasta\n")
    acc, fileNames = run.makeBlastDB(gbDir)
    output = ''.join(fileNames)
    with open(run.RootDir + "Input_GBs.txt", 'w') as f:
        f.write(output)
    rec.fin("Before BLAST", timeRecorder)
    print("\n2. Running BLAST\n")
    run.blastall(acc)
    rec.fin("BLASTexecute", timeRecorder)

    print("\n3. Parsing BLAST output\n")
    dfs=run.makeTable(acc)
    dfs.to_csv(run.dataDir + "all_blast_results.tsv", sep = "\t", index = None)
    print("Output blast result table (" + RootDir + "data/all_blast_results.tsv)")
    rec.fin("Making BLAST result tables", timeRecorder)


def runsWOblast(timeRecorder, RootDir, gbDir):
    rec = recorder.timeRecord()        
    run = Run(RootDir, gbDir)
    print("1. Skip: Creating BLASTDB and proteinFasta\n")
    print("2. Reading .gb files; Skip: Running BLAST\n")
    acc = run.makeAccList(gbDir)
    print("\n3. Parsing BLAST output\n")
    dfs=run.makeTable(acc)
    dfs.to_csv(run.dataDir + "all_blast_results.tsv", sep = "\t", index = None)
    print("Output blast result table (" + RootDir + "data/all_blast_results.tsv)")
    rec.fin("Making BLAST result tables", timeRecorder)


def clustering_woDuplicated(timeRecorder, RootDir, gbDir):
    rec = recorder.timeRecord()        
    run = Run(RootDir, gbDir)
    print("1. Creating BLASTDB and proteinFasta\n")
    acc, fileNames = run.makeBlastDB(gbDir)
    output = ''.join(fileNames)
    with open(run.RootDir + "Input_GBs.txt", 'w') as f:
        f.write(output)
    rec.fin("Before BLAST", timeRecorder)
    print("\n2. Running BLAST\n")
    run.blastall(acc)
    rec.fin("BLASTexecute", timeRecorder)

    print("\n3. Parsing BLAST output\n")
    dfs_tmp = run.makeTable(acc)
    dfs = run.clustering(dfs_tmp)
    dfs.to_csv(run.dataDir + "OrthologousGeneTag.tsv", sep = "\t")
    df_position = run.convertLocusTag2StartPosition(gbDir, dfs)
    df_position.to_csv(run.dataDir + "OrthologousGenePosition.tsv", sep = "\t")    
    print("Output blast result table (" + RootDir + "data/OrthologousGeneTag.tsv)")
    print("Output blast result table (" + RootDir + "data/OrthologousGenePosition.tsv)")
    rec.fin("Making BLAST result tables", timeRecorder)


def clustering_woDuplicated_afterBlast(timeRecorder, RootDir, gbDir):
    rec = recorder.timeRecord()
    run = Run(RootDir, gbDir)
    print("1. Skip: Creating BLASTDB and proteinFasta\n")
    print("2. Reading .gb files; Skip: Running BLAST\n")
    acc = run.makeAccList(gbDir)
    run.afterBlast(acc)
    print("\n3. Parsing BLAST output\n")
    dfs_tmp = run.makeTable(acc)
    dfs = run.clustering(dfs_tmp)
    dfs.to_csv(run.dataDir + "OrthologousGeneTag.tsv", sep = "\t")
    df_position = run.convertLocusTag2StartPosition(gbDir, dfs)
    df_position.to_csv(run.dataDir + "OrthologousGenePosition.tsv", sep = "\t")    
    print("Output blast result table (" + RootDir + "data/OrthologousGeneTag.tsv)")
    print("Output blast result table (" + RootDir + "data/OrthologousGenePosition.tsv)")
    rec.fin("Making BLAST result tables", timeRecorder)





