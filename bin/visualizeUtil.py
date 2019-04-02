# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 19:00:27 2017

@author: ac144809
"""

import pandas as pd
import createConfigs as Conf


def rotatePosition(circosIN, info, consensus_gsize, isRotate = True):
    out = circosIN.copy()
    ratio = consensus_gsize / info["Genome_size"]
    step_size = round(consensus_gsize / 360)
    out[1] = round(out[1] * ratio)
    out[2] = round(out[2] * ratio)
    if not isRotate:
       out[1] = out[1].astype(int)
       out[2] = out[2].astype(int)
       return(out)
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
        
def convertTag2Color(circosIN, Dict):
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
            
def createCircosInput(RootDir, df_tmp, consensus_gsize, CircosIN, df_info, name_key = "test", isRotate = True):
    createKaryotype(consensus_gsize, RootDir)
    temp_dict = makeConsensusDict(df_tmp)
    nrow, ncol = df_info.shape
    each_file_info = []
    for Loop in range(0,nrow):
        tmp = CircosIN[Loop]
        info = df_info.ix[Loop]
        tmp = rotatePosition(tmp, info, consensus_gsize, isRotate)
        tmp, counter = convertTag2Color(tmp,temp_dict)
        s = info[0] + " = Start position: " + str(tmp[1].min()) + ", End position: " + str(tmp[2].max()) + ", Number of removed genes:" + str(counter)
        print(s)
        each_file_info.append(s+"\n")
        file_name = RootDir + "circos/data/circos_" + name_key + "_" + info[0] +".txt"
        tmp.to_csv(file_name, sep = "\t",header=None, index=None)
        Loop += 1
 
        s = ''.join(each_file_info)
        with open(RootDir + "circos_" + name_key + "_info.txt", 'w') as f:
            f.write(s)

def createCircosConf(RootDir, df_info, sort_key = None, name_key = "test", out = "circos.conf"):
    if sort_key is None:
        df_info_sorted = df_info
    else:
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
    with open(RootDir + "/circos/" + out, 'w') as f:
        f.write(output_write)                
    df_info_sorted.to_csv(RootDir + "RingOrder_" + name_key + "_df.tsv", sep = "\t", index = None)      
    
def createKaryotype(Gsize, RootDir):
    st = []
    st.append("chr\t-\tchr1\tchr1\t0\t"+ str(int(Gsize)) + "\thue000\n")
    step = int(round(Gsize/360))
    for i in range(0,359):
        stmp = "band\tchr1\t" + str(i+1) + "\t" + str(i+1) + "\t" + str(i*step) + "\t" + str((i+1)*step) + "\thue{0:03d}\n".format(i) 
        st.append(stmp)
    i=359
    stmp = "band\tchr1\t" + str(i+1) + "\t" + str(i+1) + "\t" + str(i*step) + "\t" + str(int(Gsize)) + "\thue{0:03d}\n".format(i) 
    st.append(stmp)
    s = ''.join(st)
    with open(RootDir + "circos/data/karyotype.txt", 'w') as f:
        f.write(s) 
        
def readInfo(RootDir, df_name, CircosIN, minCount = None):
    df_info = pd.read_csv(df_name, sep = "\t")  
    df_locusTag = pd.read_csv(RootDir + '/data/locusTag_aligned_angle.tsv', delimiter = "\t")
    if minCount is not None:
        df_locusTag = df_locusTag[df_locusTag["Count"]>=minCount]
    consensus_genomesize = df_info.ix[:,1].mean()
    accList = df_info.ix[:,0]
    for acc in accList:
        tmp = pd.read_csv(RootDir + "/circos/data/" + acc + ".original.txt", header = None,delimiter = "\t")
        CircosIN.append(tmp)
    return(df_info, df_locusTag, consensus_genomesize)
