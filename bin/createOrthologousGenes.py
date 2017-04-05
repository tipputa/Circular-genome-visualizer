#! /bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 11:20:44 2017

@author: tipputa
"""

from Bio import SeqIO #, SeqFeature
import os, time, csv #, re
import pandas as pd
import blastUtil
import blast2tsv, createTable2
from tqdm import tqdm

def mkdir(dirName):
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    if not os.path.isdir(dirName):
        print("Dir name \"",dirName,"\" can't be used. Please check file names.\n")
        raise IOError

def locationToStr(locationObj):
    start = locationObj.start + 1
    end = locationObj.end
    before = "<" if "<" in str(locationObj.start) else ""
    after = ">" if ">" in str(locationObj.end) else ""

    if locationObj.strand == 1:
        return "%s%d..%s%d" % (before, start, after, end)
    elif locationObj.strand == -1:
        return "complement(%s%d..%s%d)" % (before, start, after, end)
    
def parseFeature(feature, num):
    locus_tag = feature.qualifiers.get("locus_tag", ["locus_" + str(num)])[0]
    location = locationToStr(feature.location)
    product = feature.qualifiers.get("product", ["product_" + str(num)])[0]
    function = feature.qualifiers.get("function", [""])[0]
    translation = feature.qualifiers.get("translation", [""])[0]
    return locus_tag, location, product, function, translation
                
def createProteinFasta(record, proteinFileName):
    num = 0
    with open(proteinFileName, "w") as f:
        for feature in record.features:
            num += 1
            if "translation" in feature.qualifiers.keys():
                locus_tag, location, product, function, translation = parseFeature(
                    feature, num)
                f.write(">%s %s %s\n" %
                    (locus_tag, product, record.name + ":" + location))
                f.write(feature.qualifiers["translation"][0] + "\n")
                
def createTSV(record, tsvFileName):
    csvWriter = csv.writer(
        open(tsvFileName, "w"), lineterminator="\n", delimiter='\t')
    Buffer = [("LocusTag", "sequence", "location",  "feature",
               "product", "function", "translation")]
    num = 0
    for feature in record.features:
        num += 1
        if feature.type in ["CDS", "rRNA", "tRNA"]:
            locus_tag, location, product, function, translation = parseFeature(feature, num)
            Buffer.append((locus_tag, record.name, location,
                               feature.type, product, function, translation))
    csvWriter.writerows(Buffer)

class timeRecord():
    start = time.time()
    def fin(self,str,timeRecorder):
        elapsed_time = round(time.time() - self.start,1)
        if elapsed_time > 60:
            elapsed_time_min = round(elapsed_time/60,1)
            timeRecordS = "Elapsed_time (" + str + "):{0}".format(elapsed_time_min) + "[min]"
        else:
            timeRecordS = "Elapsed_time (" + str + "):{0}".format(elapsed_time) + "[sec]"

        print("\n" + timeRecordS + "\n")
        timeRecorder.append(timeRecordS)
        self.start = time.time()
  
class Run():
    def __init__(self, RootDir, gbDir):
        print("\nSetting:\n  Root directory: %s,\n  GenBank directory: %s\n\n" % (RootDir, gbDir))
        self.RootDir = RootDir
        self.gbDir = gbDir
        self.dataDir = RootDir + "data/"
        self.CircosDir = RootDir + "circos/" 
        self.blastDBDir = self.dataDir + "BlastDB/"
        self.summaryDir = self.dataDir + "summary/"
        self.blastResultDir = self.dataDir + "BlastResult/"    
        self.proteinFasDir = self.dataDir + "ProteinFasta/"    
        mkdir(self.dataDir)
        mkdir(self.CircosDir)
        mkdir(self.CircosDir + "data")
        mkdir(self.blastDBDir)        
        mkdir(self.blastResultDir)        
        mkdir(self.proteinFasDir)
        mkdir(self.summaryDir)
      
    def makeblastdb(self, record, id="."):
        proteinFileName = self.proteinFasDir + id + ".fa"
        summaryName = self.summaryDir + id + ".tsv"
        createProteinFasta(record, proteinFileName)
        createTSV(record, summaryName)
    
        blastUtil.createDB(
            proteinFileName, self.blastDBDir + id + ".fa", dbType="prot")
    
    def blastall(self,list):
        max = len(list) ** 2 
        print("Excute Protein Blast (" + str(max) + " patterns)")
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
    
    def cTable(self,pivot, targetList):
        filterFunc = createTable2.createFilterFunction(score=0, identity=30, similarity=0, Qcoverage=0, Tcoverage=0, BBH=True)
        D = createTable2.createTable(self.blastResultDir,pivot, targetList)
        D.applymap(filterFunc)
        DL = D.apply(createTable2.LOCUSTAG, axis=1)
        DL.to_csv(self.blastResultDir + pivot + "_locusTag.tsv", sep="\t", na_rep='-',index= None)
        return(DL)
        
    def makeTable(self,acc):
        a_str = map(str,acc)
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
         
    def makeBlastDB(self,gbDir):
        fileNames = []
        acc = []
        LoopCounter = 1
    
        files = os.listdir(gbDir)
        print ("Number of GenBank files:", len([file for file in files if file.endswith(".gb")]))
        for file in files:
            if file.endswith(".gb"): 
                record = SeqIO.read(gbDir + file, "genbank")
                print(str(LoopCounter),":", record.id)            
                LoopCounter += 1
                s = gbDir + file + "\t" + record.id + "\n"
                fileNames.append(s)            
                acc.append(record.id)
                self.makeblastdb(record, record.id)
        return(acc, fileNames)

    def makeAccList(self,gbDir):
        acc = []
        LoopCounter = 1
    
        files = os.listdir(gbDir)
        print ("Number of GenBank files:", len([file for file in files if file.endswith(".gb")]))
        for file in files:
            if file.endswith(".gb"): 
                record = SeqIO.read(gbDir + file, "genbank")
                print(str(LoopCounter),":", record.id)            
                LoopCounter += 1
                acc.append(record.id)
        return(acc)
 
 
def runs(timeRecorder, RootDir, gbDir):
    rec = timeRecord()        
    run = Run(RootDir, gbDir)
    print(" Creating BLASTDB and proteinFasta\n\n")
    acc, fileNames = run.makeBlastDB(gbDir)
    output = ''.join(fileNames)
    with open(run.RootDir + "Input_GBs.txt", 'w') as f:
        f.write(output)
    rec.fin("Before BLAST", timeRecorder)
    print("\n Running BLAST\n")
    run.blastall(acc)
    rec.fin("BLASTexecute", timeRecorder)

    print("\n  Parsing BLAST output\n")
    dfs=run.makeTable(acc)
    dfs.to_csv(run.dataDir + "all_blast_results.tsv", sep = "\t", index = None)
    print("Output blast result table (./data/all_blast_results.tsv)\n")
    rec.fin("Making BLAST result tables", timeRecorder)

    return(dfs)

def runsWOblast(timeRecorder, RootDir, gbDir):
    rec = timeRecord()        
    run = Run(RootDir, gbDir)
    acc = run.makeAccList(gbDir)
    print("\n  Parsing BLAST output\n")
    dfs=run.makeTable(acc)
    dfs.to_csv(run.dataDir + "all_blast_results.tsv", sep = "\t", index = None)
    print("Output blast result table (./data/all_blast_results.tsv)\n")
    rec.fin("Making BLAST result tables", timeRecorder)

    return(dfs)

    