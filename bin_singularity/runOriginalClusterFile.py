#! /bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 11:20:44 2017

@author: tipputa
"""

import pandas as pd
import sys, os
import createCircos as createCircos
import timeRecord as recorder

useage = """

Error: please input the output directory, the GenBank directory, and the original cluster file.

Usage: 
   python runOriginalClusterFile.py <output directory> <input directory (GenBank files)> <path to the cluster file>
   e.g: runOriginalClusterFile.py ~/test/ ~/test/gb/ ~/test/data/all_blast_results.txt
"""

if __name__ == '__main__':
    timeRecorder = []
    recordAll = recorder.timeRecord()

    if len(sys.argv)==4:
        RootDir = sys.argv[1] + "/"
        gbDir = sys.argv[2] + "/"
        Cluster = sys.argv[3] 
    else:
       print(useage)
       quit()
       
    os.chdir(RootDir)    
    df = pd.read_csv(Cluster, delimiter = "\t")
    createCircos.runs(df, timeRecorder, RootDir, gbDir)
    recordAll.fin("All process", timeRecorder)
    print("\n\nFin.")
    pd.Series(timeRecorder).to_csv(RootDir + "Calculation_times.txt", header = None, index = None)
