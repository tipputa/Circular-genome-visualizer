#! /bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 11:20:44 2017

@author: tipputa
"""

import pandas as pd
import sys, time
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
   python runOriginalClusterFile.py <output directory> <input directory (GenBank files)> <path to the cluster file>
   e.g: runOriginalClusterFile.py ~/test/ ~/test/gb/ ~/test/data/all_blast_results.txt
"""

if __name__ == '__main__':
    timeRecorder = []
    recordAll = timeRecord()

    if len(sys.argv)==4:
        RootDir = sys.argv[1] + "/"
        gbDir = sys.argv[2] + "/"
        Cluster = sys.argv[3] 
    else:
       print(useage)
       quit()
        
    df = pd.read_csv(Cluster, delimiter = "\t")
    createCircos.runs(df, timeRecorder, RootDir, gbDir)
    recordAll.fin("All process", timeRecorder)
    print("\n\nFin.")
    pd.Series(timeRecorder).to_csv(RootDir + "Calculation_times.txt", header = None, index = None)
