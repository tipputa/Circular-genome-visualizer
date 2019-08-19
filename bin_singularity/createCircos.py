# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 11:20:44 2016

@author: tipputa
"""

import numpy as np
from Bio import SeqIO #, SeqFeature
import os, math #, re
import pandas as pd
import visualizeUtil as visualizer
import timeRecord as recorder

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
                tmp_out = pd.Series(["chr1", pos_start, pos_end, locus_tag[0]])
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
    
    def calcEachConsensus(self, tmp): # それぞれのclusterごとに、consensusを計算、0-360と-180~180の両方。全て"-"の場合error.
        def calcSummary(tmp_num, std_threathold = 5, median_range = 10): # 各clusterのaverage, std, median, medianMean, Finalを出力
            count = len(tmp_num)
            tmp_mean = round(tmp_num.mean())
            tmp_std = round(float(tmp_num.std()),2)
            med = int(count / 2 + 0.5) - 1
            sorted = tmp_num[tmp_num.argsort(order=False)]
            tmp_median = sorted[med]
            tmp_medianMean = round(tmp_num[(tmp_num < tmp_median + median_range) & (tmp_num > tmp_median - median_range)].mean())
            if tmp_std > std_threathold:
                Final = tmp_medianMean
            else:
                Final = tmp_mean
            L = pd.Series([tmp_mean, tmp_std, tmp_median, tmp_medianMean, Final])
            return(L)
        
        def oppositeDirection(tmp_num):
            tmp_num[tmp_num > 180] -= 360
            return(tmp_num)
            
        def correctOppositeDirection(tmp_num):
            tmp_num[tmp_num < 0] += 360
            return(tmp_num)
        
        tmp_num = tmp[tmp!="-"]
        tmp_len = len(tmp_num)
        L = pd.Series()
        if tmp_len > 1:
            L1 = calcSummary(tmp_num)
            L2 = calcSummary(oppositeDirection(tmp_num))
            if L1[1] > L2[1]:
                L2 = correctOppositeDirection(L2)
                L = L2
            else:
                L = L1
            L[5] = tmp_len        
        elif tmp_len == 1:
            Val = int(tmp_num)
            L = pd.Series([Val,0,Val,Val,Val,1])
        else:
            print("Error : there are not values, please check input file.")
            L = pd.Series([0,0,0,0,0,0])        
        return(L)
    

    def calcConsensus(self, df_num):
        tmp_out = []
        df_tmp = df_num.ix[:,self.cols[0:self.num_col]]
        for row in range(0,self.num_row):
            tmp = df_tmp.ix[row]
            tmp_out.append(self.calcEachConsensus(tmp))
    
        out = pd.DataFrame(tmp_out)
        out.columns = ["Mean","Std","Median","MedianMean","Consensus","Count"]
        df_cons = pd.concat([df_tmp,out],axis=1)
        return(df_cons)

    def calcDeviation(self, keys, consensus):
        result = math.sqrt(np.sum((keys - consensus)**2)/len(keys)-1)
        return(round(result,3))            
    
    def calcDeviations(self, df_num):
        Deviations = []
        consensus = df_num["Consensus"]
        for col in self.cols[0:self.num_col]:
            tmp = df_num[col]
            if df_num.dtypes[col]=="int64":
                df_num[col]=df_num[col].astype(object)
            res = self.calcDeviation(tmp[tmp!="-"], consensus[tmp!="-"])
            Deviations.append(res)
        return(Deviations)
    
    def AlignMinDeviation(self, df_cons):
        df_cons = df_cons.copy()
        def calcRotate(df_cons, mins = 0, maxs = 360, step = 10, strandON = 1):
            def prepareAlign(df_keys, strand = 0, angle = 0):
                df_key = df_keys.copy()
                if strand == 0:
                    df_key = df_key + angle
                else:
                    df_key = 360 - df_key + angle
                df_key[df_key >= 360] -= 360
                df_key[df_key < 0] += 360
                return(df_key)

            Deviations = []
            Strand = []
            Angle = []
            for col in self.cols[0:self.num_col]:
                df_comp = df_cons.ix[df_cons[col] != "-", [col, "Consensus"]]
                df_comp.columns = ["key","Consensus"]
                consensus = df_comp["Consensus"]
                key = df_comp["key"]
                out_strand = 0
                out_angle = 0
                devi = 0
                MinDevi = 10000
                out_df = pd.Series()
                for strand in (0, strandON):
                    for angle in range(mins,maxs,step):
                        df_prepare = prepareAlign(key, strand, angle)
                        devi = self.calcDeviation(df_prepare, consensus)
                        if MinDevi > devi:
                            out_strand = strand
                            out_angle = angle
                            out_df = df_prepare
                            MinDevi = devi
                Deviations.append(MinDevi)
                Strand.append(out_strand)
                Angle.append(out_angle)
                df_cons.ix[df_cons[col] !="-",col] = out_df
                print(col + " = strand: " + str(out_strand) + ", angle: " + str(out_angle))
            df_output = self.calcConsensus(df_cons)
            return(df_output, Deviations, Strand, Angle)
            
        print("6.1 First Alignment (rotaging by 10 degrees)")
        print("# strand: original(0) or reverse complement(1), angle: rotated degree")
        df1, tmpDev, tmpStr, tmpAng = calcRotate(df_cons)
        print("\n6.2 Second Alignment (rotaging by 1 degree)")
        print("# strand: original(0) or reverse complement(1), angle: rotated degree")
        df2, tmpDev, tmpStr2, tmpAng2 = calcRotate(df1, -20, 20, 1, 0)
        outAng = pd.Series(tmpAng) + pd.Series(tmpAng2)
        outAng[outAng > 360] -= 360
        outAng[outAng < 0] += 360
        return(df2, tmpDev, tmpStr, outAng)
        
def runs(df, timeRecorder, RootDir, gbDir):
    rec = recorder.timeRecord()        
    run = Runs(RootDir, gbDir)
    df_num = df.copy()
    files = os.listdir(gbDir)
    
    genome_size = []
    CircosIN = []
    LoopCounter = 1
    print ("4. Reading GenBank files again")
    NumGB = len([file for file in files if file.endswith(".gb")  or file.endswith('.gbk')])
    print ("Number of GenBank files:", NumGB)
    for file in files:
        if file.endswith(".gb") or file.endswith('.gbk'): 
            print("Reading: " + file + "  "+ str(LoopCounter) + "/" + str(NumGB)),
            LoopCounter += 1
            df_num = run.readGB(gbDir + file, genome_size, CircosIN, df_num)

    run.num_row, run.num_col = df_num.shape
    num_row, num_col = df_num.shape
    cols = df_num.columns
    run.cols = df_num.columns
    AccNo = list(cols[0:num_col])
    print ("\n5. Calculating consensus and deviation\n")
    df_cons = run.calcConsensus(df_num)
    Dev_original = run.calcDeviations(df_cons)
    rec.fin("Before alignment", timeRecorder)
    print ("6. Alignment by minimizing the deviation\n")
    df_cons2, Dev, Strand, Angle = run.AlignMinDeviation(df_cons)
    rec.fin("Alignment", timeRecorder)

    print ("\n7. file output\n")        
    df_cons.to_csv(RootDir + "data/output_unaligned_consensus.tsv", sep = "\t")    
    df_cons2.to_csv(RootDir + "data/output_aligned_consensus.tsv", sep = "\t")
    consensus_gsize = round(np.mean(genome_size))
    df_info = pd.concat([pd.Series(AccNo),pd.Series(genome_size),pd.Series(Strand),Angle,pd.Series(Dev), pd.Series(Dev_original)], axis = 1)
    df_info.columns = ["AccNo", "Genome_size", "Strand", "Angle", "Deviation (Aligned)", "Deviation (Original)"]
    df_info.to_csv(RootDir + "data/output_information.tsv", sep = "\t", index = None)
    df_locusTag = pd.concat([df, df_cons[["Consensus", "Std", "Count"]]], axis = 1)
    df_locusTag_aligned = pd.concat([df, df_cons2[["Consensus", "Std", "Count"]]], axis = 1)
    df_locusTag.to_csv(RootDir + "data/locusTag_unaligned_angle.tsv", sep = "\t", index = None)
    df_locusTag_aligned.to_csv(RootDir + "data/locusTag_aligned_angle.tsv", sep = "\t", index = None)
    df_locusTag_aligned_core = df_locusTag_aligned[df_locusTag_aligned["Count"]==num_col]
    df_locusTag_aligned_core.to_csv(RootDir + "data/locusTag_aligned_angle_coreGenes.tsv", sep = "\t", index = None)
    
    print ("8. Running circos\n")
    os.chdir(RootDir)
    
    def visualize(dfs, sort_key, tag, df_info, isRotate = True):        
        print("Creating circos input (" + RootDir + "./circos/data/circos_" + tag + "_*.txt)")
        visualizer.createCircosInput(RootDir, dfs, consensus_gsize, CircosIN, df_info, tag, isRotate)                      
        visualizer.createCircosConf(RootDir, df_info, sort_key, tag)
        print("Run circos; suffix: "  + tag )
        os.system('singularity exec /usr/local/biotools/c/circos:0.69.5--pl5.22.0_0 circos -conf '+ RootDir +'circos/circos.conf -silent')
        os.system("mv circos.png " + RootDir + "circos_" + tag + ".png")
        os.system("rm circos.svg")
        print("Output circos_" + tag + ".png, circos_" + tag + "_info.txt, and RingOrder_" + tag + "_df.tsv\n")

    visualize(df_locusTag, "Deviation (Original)", "unaligned", df_info, False)
    visualize(df_locusTag_aligned, "Deviation (Aligned)", "aligned", df_info)
    visualize(df_locusTag_aligned, "Deviation (Original)", "aligned_UnalignedRingOrder", df_info)
    visualize(df_locusTag_aligned_core, "Deviation (Aligned)", "aligned_core", df_info)

    rec.fin("Running Circos", timeRecorder)
