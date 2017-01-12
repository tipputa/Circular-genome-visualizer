# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 11:20:44 2016

@author: tipputa
"""

import numpy as np
from Bio import SeqIO #, SeqFeature
import os, math, time #, re
import pandas as pd

class Runs():
    def __init__(self, binDir, RootDir, gbDir):
        self.binDir = binDir
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
                    print("Error: the position of first gene is modified")
                    tmp_out[2] = tmp_out[1] + 100
                    
                cdsL.append(tmp_out)
                cds_angle = round (feature.location.start / step)
                locusDict.update({locus_tag[0]:cds_angle})
    
        while len(target) != 0:
            tmp = target.pop()
            newList.append(locusDict[tmp])
    
        newList.reverse()
        df[str(record.id)] = pd.Series(newList)
        CircosIN.append(pd.DataFrame(cdsL))
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
        Loop = 1
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
        print("Start: " + str(Loop) + " Loop")
        Loop += 1
        df1, tmpDev, tmpStr, tmpAng = calcRotate(df_cons)
        print("Start: " + str(Loop) + " Loop")
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
        
        
            
    def makeConsensusDict(self, df):
        rows, cols = df.shape
        colname = df.columns
        Dict = {}
        
        for col in colname[0:cols-1]:
            tmp = df.ix[df[col] != "-",[col,"Consensus"]]
            tmp.index = tmp[col]
            Dict.update(tmp["Consensus"])
        return(Dict)
            
    def convertTag2Color(self, circosIN, Dict):
        out = circosIN.copy()
        L = []
        count = 0
        for tag in out[3]:
            if tag in Dict:
                L.append("fill_color=hue{0:03d}".format(int(Dict[tag])))
            else:
                L.append("fill_color=white")
                count += 1
        out[3] = pd.Series(L)
    
    
        return(out, count)
                
    def createCircosInput(self, df_tmp, name_key,consensus_gsize, CircosIN, df_info):
        CircosIN_rotated = []
        temp_dict = self.makeConsensusDict(df_tmp)
        for Loop in range(0,self.num_col):
            tmp = CircosIN[Loop]
            info = df_info.ix[Loop]
            tmp = self.rotatePosition(tmp, info, consensus_gsize)
            tmp, counter = self.convertTag2Color(tmp,temp_dict)
            print (info[0] + " = start posi.: " + str(tmp[1].min()) + ",End posi.: " + str(tmp[2].max()) + ", Number of removal genes:" + str(counter))
            file_name = self.CircosDir + "data/circos_" + name_key + "_" + info[0] +".txt"
            CircosIN_rotated.append(tmp)
            tmp.to_csv(file_name, sep = "\t",header=None, index=None)
            Loop += 1
    
    def createCircosConf(self, sort_key, name_key, df_info, out = "circos.conf"):
        Acc_sorted = df_info.sort_values(by=[sort_key], ascending=False).ix[:,0]
        Count, = Acc_sorted.shape
        step = 0.9/Count
        r1 = 0.99
        r0 = 0.99 - step
        f = open(self.binDir+"circosConf_head.txt")
        head = f.read()
        f.close
        output = []
        output.append(head)
        for acc in Acc_sorted:
            s = "<highlight>\nfile = data/circos_" + name_key + "_" + acc + ".txt\nr1 = " + str(r1) + "r\nr0 = " + str(r0) + "r\n</highlight>\n"         
            output.append(s)
            r1 = r0
            r0 = r1 - step
        f2 = open(self.binDir+"circosConf_tail.txt")
        tail = f2.read()
        f2.close()
        output.append(tail)
        output_write = ''.join(output)
            
        
        f3 = open(self.CircosDir  + out, 'w')
        f3.write(output_write)
        f3.close
    
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
        
def runs(df, timeRecorder, binDir, RootDir, gbDir):
    run = Runs(binDir, RootDir, gbDir)
    df_num = df.copy()
    files = os.listdir(gbDir)

    genome_size = []
    CircosIN = []
    LoopCounter = 1
    print (len([file for file in files if file.endswith(".gb")]))
    for file in files:
        if file.endswith(".gb"): 
            print(str(LoopCounter) + "..."),
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
    df_cons2.to_csv(RootDir + "data/output_consensus.tsv", sep = "\t")
        
    consensus_gsize = round(np.mean(genome_size))
    run.createKaryotype(consensus_gsize)
    df_info = pd.concat([pd.Series(AccNo),pd.Series(genome_size),pd.Series(Strand),Angle,pd.Series(Dev), pd.Series(Dev_original)], axis = 1)
    df_info.columns = ["AccNo", "Genome_size", "Strand", "Angle", "Deviation (Aligned)", "Deviation (Original)"]
    df_info.to_csv(RootDir + "data/output_information.tsv", sep = "\t", index=None)
    df_locusTag_aligned = pd.concat([df, df_cons2["Consensus"]], axis = 1)
    df_locusTag = pd.concat([df, df_cons["Consensus"]], axis = 1)

    os.chdir(RootDir + "circos/")
    run.createCircosInput(df_locusTag, "unaligned", consensus_gsize, CircosIN, df_info)
    run.createCircosConf("Deviation (Original)", "unaligned", df_info)
    os.system('circos')     
    os.system('mv circos.png ../circos_unaligned.png')
    run.createCircosInput(df_locusTag_aligned, "aligned", consensus_gsize, CircosIN, df_info)
    run.createCircosConf("Deviation (Aligned)", "aligned", df_info)
    os.system('circos')     
    os.system('mv circos.png ../circos_aligned.png')

    os.chdir(RootDir + "circos/")
    os.system('circos') 
if __name__ == '__main__':

    timeRecorder = []
    binDir = "/Users/tipputa/Google/1_study/spyder/bin/"
    RootDir = "/Users/tipputa/Google/1_study/1_circos_for_paper/circos_testcase/170105_Hpylori_allPython/"
  #  gbDir = "/Users/tipputa/Google/1_study/1_circos_for_paper/circos_testcase/160422_Hpylori_all/gb/"
    gbDir = "/Users/tipputa/Google/1_study/1_circos_for_paper/circos_testcase/170105_Hpylori_allPython/gb/"
   # df = pd.read_csv(RootDir + '/data/input_rmDup.tsv', delimiter = "\t")
    df = pd.read_csv("/Users/tipputa/Google/1_study/spyder/input.tsv", delimiter = "\t")
    runs(df, timeRecorder, binDir, RootDir, gbDir)