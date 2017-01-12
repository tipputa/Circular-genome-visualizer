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
        elapsed_time = time.time() - self.start
        timeRecordS = str + "elapsed_time :{0}".format(elapsed_time) + "[sec]"
        print(timeRecordS)
        timeRecorder.append(timeRecordS)
        
        self.start = time.time()
  
class Run():
    def __init__(self, binDir, RootDir, gbDir):
        print("Setting: Bin %s,\n Root %s,\n gb %s" % (binDir, RootDir, gbDir))
        self.binDir = binDir
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
        print("Excute Protein Blast (",max,"patterns)")
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
        return(output)
        
    def check(self,a,b):
        for i in range(0,len(a)):
            if not b[i] == "-":
                if a[i] == "-":
                    a[i] = b[i]
                if isinstance(a[i],list):
                    if not b[i] in a[i]:
                        a[i].append(b[i])
                elif isinstance(b[i], list):
                    if not a[i] in b[i]:
                        b[i].append(a[i])
                else:
                    if(a[i]!=b[i]):
                        a[i] = [a[i],b[i]]
        return(a)
    
    #       "{0:03d}".format(1)
    
    def checkDup(self,df,i):
        a = df.ix[:,i]
        nrow,  = a.shape 
        pool = []
        dup = []
        for i in range(0,nrow):
            if not a[i] == "-":
                if isinstance(a[i], list):
                    b = a[i].copy()
                    while b:
                        tmp = b.pop()
                        if tmp in pool:
                            if not tmp in dup:
                                dup.append(tmp)
                        else:
                            pool.append(tmp)
                else:
                    if a[i] in pool:
                        if not a[i] in dup:
                        #    print("duplicated",a[i])
                            dup.append(a[i])
                    else:
                        pool.append(a[i])
        return(dup)
    
    def rmDuplication(self,d, dup, col):
        while dup:
            tmp = dup.pop()
            tag = d.ix[d.ix[:,col]==tmp,].index
            target = d.ix[tag].copy()
            new = self.rmDup(target)
            d = d.drop(tag)
            d = d.append(new)        
        return(d)
        
    def rmDup(self,a):
        nrow, col = a.shape
        out = a.ix[0]
        for i in range(1,nrow):
            self.check(out,a.ix[i])
        return(out)
    
    def removeDuplications(self,df):    
        nrow, ncol = df.shape
        for i in range(0,ncol):
            dup = self.checkDup(df,i)
            df = self.rmDuplication(df, dup, i)
            
        def removeHyphen(df):
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
        df = removeHyphen(df)
        return(df)
        
        
    def checkList(self,df):
        nrow, ncol = df.shape
        rem = []
        for i in range(0,nrow):
            for j in range(0,ncol):
                if isinstance(df.ix[i,j],list):
                    if not df.ix[i].name in rem:
                        rem.append(df.ix[i].name)
                    break
        return(rem)
    
            
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

        
def runs(timeRecorder, binDir, RootDir, gbDir):
    rec = timeRecord()        
    run = Run(binDir, RootDir, gbDir)
    acc, fileNames = run.makeBlastDB(gbDir)
    output = ''.join(fileNames)
    with open(run.RootDir + "Input_GBs.txt", 'w') as f:
        f.write(output)
    rec.fin("Before BLAST, ", timeRecorder)
    run.blastall(acc)
    rec.fin("BLAST, ", timeRecorder)

    dfs=run.makeTable(acc)

    rec.fin("making BLAST result tables, ", timeRecorder)

    dfs_rmDup = run.removeDuplications(dfs)
    Checks = run.checkList(dfs_rmDup)
    dfs_Duplications = dfs_rmDup.ix[Checks]
    dfs_rmDupAll = dfs_rmDup.drop(Checks)
    dfs_rmDup.to_csv(run.dataDir + "input_including_duplicates.tsv", sep = "\t", index = None)
    dfs_Duplications.to_csv(run.dataDir + "duplicates.tsv", sep = "\t", index = None)
    dfs_rmDupAll.to_csv(run.dataDir + "input_rmDup.tsv", sep = "\t", index = None)

    print("Number of removed gene clusters : ",len(Checks))    

    rec.fin("After BLAST, ", timeRecorder)    
    return(dfs_rmDupAll)

def runsWOblast(timeRecorder, binDir, RootDir, gbDir):
    rec = timeRecord()        
    run = Run(binDir, RootDir, gbDir)
    acc, fileNames = run.makeBlastDB(gbDir)
    output = ''.join(fileNames)
    with open(run.RootDir + "Input_GBs.txt", 'w') as f:
        f.write(output)
    rec.fin("Before BLAST, ", timeRecorder)
  #  run.blastall(acc)
   # rec.fin("BLAST, ", timeRecorder)

    dfs=run.makeTable(acc)

    rec.fin("making BLAST result tables, ", timeRecorder)

    dfs_rmDup = run.removeDuplications(dfs)
    Checks = run.checkList(dfs_rmDup)
    dfs_Duplications = dfs_rmDup.ix[Checks]
    dfs_rmDupAll = dfs_rmDup.drop(Checks)
    dfs_rmDup.to_csv(run.dataDir + "input_including_duplicates.tsv", sep = "\t", index = None)
    dfs_Duplications.to_csv(run.dataDir + "duplicates.tsv", sep = "\t", index = None)
    dfs_rmDupAll.to_csv(run.dataDir + "input_rmDup.tsv", sep = "\t", index = None)
    
    print("Number of removed gene clusters : ",len(Checks))    

    rec.fin("After BLAST, ", timeRecorder)    
    return(dfs_rmDupAll)

if __name__ == '__main__':
    recordAll = timeRecord()
    rec = timeRecord()
    
    timeRecorder = []
    binDir = "/Users/tipputa/Google/1_study/spyder/bin/"
    RootDir = "/Users/tipputa/Google/1_study/1_circos_for_paper/circos_testcase/170105_Hpylori_allPython/"
  #  gbDir = "/Users/tipputa/Google/1_study/1_circos_for_paper/circos_testcase/160422_Hpylori_all/gb/"
    gbDir = "/Users/tipputa/Google/1_study/1_circos_for_paper/circos_testcase/170105_Hpylori_allPython/gb/"
    
 #   runs(timeRecorder, binDir, RootDir, gbDir)
    runsWOblast(timeRecorder, binDir, RootDir, gbDir)

    pd.Series(timeRecord).to_csv(RootDir +   "data/Calculation_times.txt", header = None, index = None)
    