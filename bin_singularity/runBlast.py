#! /bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 11:20:44 2017

@author: tipputa
"""

from Bio import SeqIO
import os, time, sys
import pandas as pd
import blastUtil, blast2tsv, createTable2
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
      
    
    def blastall(self, list, binDir):
        command=" ".join(["qsub ", binDir + "/qblast_all.sh ",self.blastDBDir,self.proteinFasDir, self.blastResultDir, ",".join(list)])
        for id in list:
            command=" ".join(["qsub ", binDir + "/qblast_all.sh ",binDir,id,self.proteinFasDir,self.blastDBDir, self.blastResultDir, ",".join(list)])
            os.system(command)
 
def runs(timeRecorder, RootDir, gbDir, binDir):
    rec = recorder.timeRecord()        
    run = Run(RootDir, gbDir)
    print("1. Creating BLASTDB and proteinFasta\n")
    acc, fileNames = run.makeBlastDB(gbDir)
    output = ''.join(fileNames)
    with open(run.RootDir + "Input_GBs.txt", 'w') as f:
        f.write(output)
    rec.fin("Before BLAST", timeRecorder)
    print("\n2. Running BLAST\n")
    run.blastall(acc, binDir)
    rec.fin("BLASTexecute", timeRecorder)


