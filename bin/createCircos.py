# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 11:20:44 2016

@author: tipputa
"""

import numpy as np
from Bio import SeqIO #, SeqFeature
import os, math #, re
import pandas as pd
import createConfigs as Conf

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
        tmp_num = tmp[tmp!="-"]
        tmp_len = len(tmp_num)
        L = pd.Series()
        if tmp_len > 1:
            L1 = self.calcSummary(tmp_num)
            L2 = self.calcSummary(self.oppositeDirection(tmp_num))
            if L1[1] > L2[1]:
                L2 = self.correctOppositeDirection(L2)
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
    
    def calcSummary(self, tmp_num, std_threathold = 5, median_range = 10):
        count = len(tmp_num)
        tmp_mean = round(tmp_num.mean())
        tmp_std = round(float(tmp_num.std()),2)
        med = int(count / 2) - 1
        sorted = tmp_num[tmp_num.argsort(order=False)]
        tmp_median = sorted[med]
        tmp_medianMean = round(tmp_num[(tmp_num < tmp_median + median_range) & (tmp_num > tmp_median - median_range)].mean())
        if tmp_std > std_threathold:
            Final = tmp_medianMean
        else:
            Final = tmp_mean
        L = pd.Series([tmp_mean, tmp_std, tmp_median, tmp_medianMean,Final])
        return(L)
    
    def oppositeDirection(self, tmp_num):
        tmp_num[tmp_num > 180] -= 360
        return(tmp_num)
        
    def correctOppositeDirection(self, tmp_num):
        tmp_num[tmp_num < 0] += 360
        return(tmp_num)
    
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
        
    
    def calcDeviation(self,keys, consensus):
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
        def calcRotate(df_cons,mins=0,maxs=360,step=10,strandON=1):
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
                for strand in (0,strandON):
                    for angle in range(mins,maxs,step):
                        df_prepare = self.prepareAlign(key, strand, angle)
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
                print("col: " + col + ", strand: " + str(out_strand) + ", angle: " + str(out_angle))
            df_output = self.calcConsensus(df_cons)
            return(df_output, Deviations, Strand, Angle)
        print("\n First Alignment (rotaging by 10 degrees)\n")
        df1, tmpDev, tmpStr, tmpAng = calcRotate(df_cons)
        print("\n Second Alignment (rotaging by 1 degree)\n")
        df2, tmpDev, tmpStr2, tmpAng2 = calcRotate(df1,-20,20,1,0)
        outAng = pd.Series(tmpAng) + pd.Series(tmpAng2)
        outAng[outAng > 360] -= 360
        outAng[outAng < 0] += 360
        return(df2, tmpDev, tmpStr, outAng)
        
    def prepareAlign(self, df_keys, strand = 0, angle = 0):
        df_key = df_keys.copy()
        if strand == 0:
            df_key = df_key + angle
        else:
            df_key = 360 - df_key + angle
        df_key[df_key >= 360] -= 360
        df_key[df_key < 0] += 360
        return(df_key)

    def rotatePosition(self, circosIN, info, consensus_gsize):
        out = circosIN.copy()
        ratio = consensus_gsize / info["Genome_size"]
        step_size = round(consensus_gsize / 360)
        out[1] = round(out[1] * ratio)
        out[2] = round(out[2] * ratio)
        tmp = out[1].copy()
        if info["Strand"] == 0:
            if info["Angle"] != 0:
                out[1] = out[1] + step_size * info["Angle"]
                out[2] = out[2] + step_size * info["Angle"]
                out.ix[out[1] > consensus_gsize, 1] -= consensus_gsize 
                out.ix[out[2] > consensus_gsize, 2] -= consensus_gsize
        else:
            if info["Angle"] == 0:
                out[1] = consensus_gsize - out[2]
                out[2] = consensus_gsize - tmp
            else:            
                out[1] = consensus_gsize - out[2] + step_size * info["Angle"]
                out[2] = consensus_gsize - tmp + step_size * info["Angle"]
                out.ix[out[1] > consensus_gsize, 1] -= consensus_gsize 
                out.ix[out[2] > consensus_gsize, 2] -= consensus_gsize
        out.ix[out[1] > out[2],2] = consensus_gsize
        out[1] = out[1].astype(int)
        out[2] = out[2].astype(int)
        return(out)
     
    def AlignGsize(self, circosIN, info, consensus_gsize):
        out = circosIN.copy()
        ratio = consensus_gsize / info["Genome_size"]
        out[1] = round(out[1] * ratio)
        out[2] = round(out[2] * ratio)
        out.ix[out[1] > out[2],2] = consensus_gsize
        out[1] = out[1].astype(int)
        out[2] = out[2].astype(int)
        return(out)
                
            
    def makeConsensusDict1(self, df):
        rows, cols = df.shape
        colname = df.columns
        Dict = {}
        
        for col in colname[0:cols-3]:
            tmp = df.ix[df[col] != "-",[col,"Consensus"]]
            tmp.index = tmp[col]
            tmp=tmp.drop_duplicates([col]) # First duplicated value is remaining
            Dict.update(tmp["Consensus"]) 
        return(Dict)
        
    def makeConsensusDict2(self, df):
        rows, cols = df.shape
        colname = df.columns
        Dict = {}
        
        for col in colname[0:cols-3]:
            tmp = df.ix[df[col] != "-",[col,"Consensus"]]
            tmp.index = tmp[col]
            tmp=tmp[~tmp.duplicated([col])]  # Duplicated values are removed
            Dict.update(tmp["Consensus"])
        return(Dict)
            
    def convertTag2Color(self, circosIN, Dict):
        out = circosIN.copy()
        L = []
        count1 = 0
        count2 = 0
        for tag in out[3]:
            if tag in Dict:
                L.append("fill_color=hue{0:03d}".format(int(Dict[tag])))
                count1 += 1
            else:
                L.append("fill_color=white")
                count2 += 1
        out[3] = pd.Series(L)    
        return(out, count1, count2)
                
    def createCircosInput(self, df_tmp, name_key,consensus_gsize, CircosIN, df_info):
        temp_dict = self.makeConsensusDict2(df_tmp)
        each_file_info = []
        for Loop in range(0,self.num_col):
            tmp = CircosIN[Loop]
            info = df_info.ix[Loop]
            tmp = self.rotatePosition(tmp, info, consensus_gsize)
            tmp, cRemaining, cOmitted = self.convertTag2Color(tmp,temp_dict)
            s = info[0] + ", start bp: " + str(tmp[1].min()) + ", End bp: " + str(tmp[2].max()) + ", Number of visualized genes: " + str(cRemaining) + ", Number of removal genes: " + str(cOmitted)
            print(s)
            each_file_info.append(s + "\n")
            file_name = self.CircosDir + "data/circos_" + name_key + "_" + info[0] +".txt"
            tmp.to_csv(file_name, sep = "\t",header=None, index=None)
            Loop += 1
        
        s = ''.join(each_file_info)
        with open(self.RootDir + "circos_" + name_key + "_info.txt", 'w') as f:
            f.write(s)
    def createCircosInput2(self, df_tmp, name_key,consensus_gsize, CircosIN, df_info):
        temp_dict = self.makeConsensusDict2(df_tmp)
        each_file_info = []
        for Loop in range(0,self.num_col):
            tmp = CircosIN[Loop]
            info = df_info.ix[Loop]
            tmp = self.AlignGsize(tmp, info, consensus_gsize)
            tmp, cRemaining, cOmitted = self.convertTag2Color(tmp,temp_dict)
            s = info[0] + ", start bp: " + str(tmp[1].min()) + ", End bp: " + str(tmp[2].max()) + ", Number of visualized genes: " + str(cRemaining) + ", Number of removal genes: " + str(cOmitted)
            print(s)
            each_file_info.append(s + "\n")
            file_name = self.CircosDir + "data/circos_" + name_key + "_" + info[0] +".txt"
            tmp.to_csv(file_name, sep = "\t",header=None, index=None)
            Loop += 1
        
        s = ''.join(each_file_info)
        with open(self.RootDir + "circos_" + name_key + "_info.txt", 'w') as f:
            f.write(s)

             
    def createCircosConf(self, sort_key, name_key, df_info, out = "circos.conf"):
        df_info_sorted = df_info.sort_values(by=[sort_key], ascending=False)
        Acc_sorted = df_info_sorted.ix[:,0]
        Count, = Acc_sorted.shape
        step = round(0.9/Count,2)
        r1 = 0.99
        r0 = 0.99 - step
        head = Conf.header
        output = []
        output.append(head)
        for acc in Acc_sorted:
            s = "<highlight>\nfile = data/circos_" + name_key + "_" + acc + ".txt\nr1 = " + str(r1) + "r\nr0 = " + str(r0) + "r\n</highlight>\n"         
            output.append(s)
            r1 = r0
            r0 = round(r1 - step,2)
        tail = Conf.tail
        output.append(tail)
        output_write = ''.join(output)
                    
        f = open(self.CircosDir  + out, 'w')
        f.write(output_write)
        f.close
        df_info_sorted.to_csv(self.RootDir + "RingOrder_" + name_key + "_df.tsv", sep = "\t", index = None)      
        
    def createKaryotype(self, max):
        st = []
        st.append("chr\t-\tchr1\tchr1\t0\t"+ str(int(max)) + "\thue000\n")
        step = int(round(max/360))
        for i in range(0,359):
            stmp = "band\tchr1\t" + str(i+1) + "\t" + str(i+1) + "\t" + str(i*step) + "\t" + str((i+1)*step) + "\thue{0:03d}\n".format(i) 
            st.append(stmp)
        i=359
        stmp = "band\tchr1\t" + str(i+1) + "\t" + str(i+1) + "\t" + str(i*step) + "\t" + str(int(max)) + "\thue{0:03d}\n".format(i) 
        st.append(stmp)
        s = ''.join(st)
        with open(self.CircosDir + "data/karyotype.txt", 'w') as f:
            f.write(s) 


def runs(df, timeRecorder, RootDir, gbDir):
    run = Runs(RootDir, gbDir)
    df_num = df.copy()
    files = os.listdir(gbDir)

    genome_size = []
    CircosIN = []
    LoopCounter = 1
    NumGB = len([file for file in files if file.endswith(".gb")])
    print ("Number of GenBank files:", NumGB)
    for file in files:
        if file.endswith(".gb"): 
            print("Reading: " + file + "  "+ str(LoopCounter) + "/" + str(NumGB)),
            LoopCounter += 1
            df_num = run.readGB(gbDir + file, genome_size, CircosIN, df_num)

    run.num_row, run.num_col = df_num.shape
    num_row, num_col = df_num.shape
    cols = df_num.columns
    run.cols = df_num.columns
    AccNo = list(cols[0:num_col])
    df_cons = run.calcConsensus(df_num)
    Dev_original = run.calcDeviations(df_cons)
    df_cons2, Dev, Strand, Angle = run.AlignMinDeviation(df_cons)
    df_cons.to_csv(RootDir + "data/output_unaligned_consensus.tsv", sep = "\t")    
    df_cons2.to_csv(RootDir + "data/output_aligned_consensus.tsv", sep = "\t")
        
    consensus_gsize = round(np.mean(genome_size))
    run.createKaryotype(consensus_gsize)
    df_info = pd.concat([pd.Series(AccNo),pd.Series(genome_size),pd.Series(Strand),Angle,pd.Series(Dev), pd.Series(Dev_original)], axis = 1)
    df_info.columns = ["AccNo", "Genome_size", "Strand", "Angle", "Deviation (Aligned)", "Deviation (Original)"]
    df_info.to_csv(RootDir + "data/output_information.tsv", sep = "\t", index = None)
    df_locusTag = pd.concat([df, df_cons[["Consensus", "Std", "Count"]]], axis = 1)
    df_locusTag_aligned = pd.concat([df, df_cons2[["Consensus", "Std", "Count"]]], axis = 1)
    df_locusTag.to_csv(RootDir + "data/locusTag_unaligned_angle.tsv", sep = "\t", index = None)
    df_locusTag_aligned.to_csv(RootDir + "data/locusTag_aligned_angle.tsv", sep = "\t", index = None)
    df_locusTag_aligned_core = df_locusTag_aligned[df_locusTag_aligned["Count"]==num_col]
    df_locusTag_aligned_core.to_csv(RootDir + "data/locusTag_aligned_angle_coreGenes.tsv", sep = "\t", index = None)
    
    os.chdir(RootDir + "circos/")
    print("\n\nCreate Circos input (./circos/data/).")
    run.createCircosInput2(df_locusTag, "unaligned", consensus_gsize, CircosIN, df_info)
    run.createCircosConf("Deviation (Original)", "unaligned", df_info)
    print("Run circos for unaligned dataset.")
    os.system('circos')     
    os.system('mv circos.png ../circos_unaligned.png')
    print("Output circos_unaligned.png, circos_unaligned_info.txt, and RingOrder_unaligned_df.tsv")
    
    print("\nCreate Circos input (./circos/data/).")
    run.createCircosInput(df_locusTag_aligned, "aligned", consensus_gsize, CircosIN, df_info)
    run.createCircosConf("Deviation (Aligned)", "aligned", df_info)
    print("Run circos for aligned dataset")
    os.system('circos')     
    os.system('mv circos.png ../circos_aligned.png')
    print("Output circos_aligned.png, circos_aligned_info.txt, and RingOrder_aligned_df.tsv")

    print("\nCreate Circos input (./circos/data/).")
    run.createCircosInput(df_locusTag_aligned, "aligned_UnalignedRingOrder", consensus_gsize, CircosIN, df_info)
    run.createCircosConf("Deviation (Original)", "aligned_UnalignedRingOrder", df_info)
    print("Run circos for aligned dataset (Unaligned ring order)")
    os.system('circos')     
    os.system('mv circos.png ../circos_aligned_UnalignedRingOrder.png')
    print("Output circos_aligned_UnalignedRingOrder.png, circos_aligned_UnalignedRingOrder_info.txt, and RingOrder_aligned_UnalignedOrder_df.tsv")

    print("\nCreate Circos input (./circos/data/).")
    run.createCircosInput(df_locusTag_aligned_core, "aligned_core", consensus_gsize, CircosIN, df_info)
    run.createCircosConf("Deviation (Aligned)", "aligned_core", df_info)
    print("Run circos for aligned core gene dataset")
    os.system('circos')     
    os.system('mv circos.png ../circos_aligned_core.png')
    print("Output circos_aligned_core.png, circos_aligned_core_info.txt, and RingOrder_aligned_core_df.tsv")

    