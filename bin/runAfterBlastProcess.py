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
        elapsed_time = round(time.time() - self.start,1)
        if elapsed_time > 60:
            elapsed_time_min = round(elapsed_time/60,1)
            timeRecordS = "Elapsed_time (" + str + "):{0}".format(elapsed_time_min) + "[min]"
        else:
            timeRecordS = "Elapsed_time (" + str + "):{0}".format(elapsed_time) + "[sec]"

        print("\n" + timeRecordS + "\n")
        timeRecorder.append(timeRecordS)
        self.start = time.time()
  
useage = """

Error: please input results and genbank directories.

Usage: 
   python runAfterBlastProcess.py <output directory> <input directory (GenBank files)>
   e.g: python runAllProcess.py ~/test/ ~/test/gb/
"""

if __name__ == '__main__':
    timeRecorder = []
    recordAll = timeRecord()

    if len(sys.argv)==3:
        RootDir = sys.argv[1] + "/"
        gbDir = sys.argv[2] + "/"
        
    else:
       print(useage)
       quit()
        
    createOrtho.runsWOblast(timeRecorder, RootDir, gbDir)
    df = pd.read_csv(RootDir + '/data/all_blast_results.tsv', delimiter = "\t")
    createCircos.runs(df, timeRecorder, RootDir, gbDir)
    recordAll.fin("All process", timeRecorder)
    print("\n\nFin.")
    pd.Series(timeRecorder).to_csv(RootDir + "Calculation_times.txt", header = None, index = None)
