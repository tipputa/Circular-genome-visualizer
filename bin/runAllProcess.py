#! /bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 11:20:44 2017

@author: tipputa
"""

import pandas as pd
import sys, time
import createOrthologousGenes as createOrtho
import createCircos as createCircos

class timeRecord():
    start = time.time()

    def fin(self,str,timeRecorder):
        elapsed_time = time.time() - self.start
        timeRecordS = str + "elapsed_time :{0}".format(elapsed_time) + "[sec]"
        print(timeRecordS)
        timeRecorder.append(timeRecordS)
        
        self.start = time.time()

if __name__ == '__main__':
    recordAll = timeRecord()
    rec = timeRecord()
    
    timeRecorder = []
    if len(sys.argv)==4:
        binDir = sys.argv[1]
        RootDir = sys.argv[2]
        gbDir = sys.argv[3]
    else:
        binDir = "/Users/tipputa/Google/1_study/spyder/bin/"
        RootDir = "/Users/tipputa/Google/1_study/1_circos_for_paper/circos_testcase/170105_Hpylori_allPython/"
  #      gbDir = "/Users/tipputa/Google/1_study/1_circos_for_paper/circos_testcase/160422_Hpylori_all/gb/"
        gbDir = "/Users/tipputa/Google/1_study/1_circos_for_paper/circos_testcase/170105_Hpylori_allPython/gb/"

    createOrtho.runs(timeRecorder, binDir, RootDir, gbDir)
#    df = createOrtho.runsWOblast(timeRecorder, binDir, RootDir, gbDir)
    rec.fin("After All, ", timeRecorder)
    df = pd.read_csv(RootDir + '/data/input_rmDup.tsv', delimiter = "\t")
    createCircos.runs(df, timeRecorder, binDir, RootDir, gbDir)
    recordAll.fin("All process:, ", timeRecorder)

    pd.Series(timeRecord).to_csv(RootDir + "Calculation_times.txt", header = None, index = None)
    