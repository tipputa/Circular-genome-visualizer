#! /bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 11:20:44 2017

@author: tipputa
"""

import pandas as pd
import sys, os
import createOrthologousGenes as createOrtho
import createCircos as createCircos
import timeRecord as recorder

useage = """

Error: please input the absolute PATH of output and genbank directories .

Usage: 
   python runAllProcess.py <output directory> <input directory (GenBank files)>
   e.g: python runAllProcess.py ~/study/ ~/study/gb/
"""

if __name__ == '__main__':
    timeRecorder = []
    recordAll = recorder.timeRecord()

    if len(sys.argv)==3:
        RootDir = sys.argv[1] + "/"
        gbDir = sys.argv[2] + "/"
        
    else:        
        print(useage)
        quit()
    
    os.chdir(RootDir)
    createOrtho.runs(timeRecorder, RootDir, gbDir)
    df = pd.read_csv(RootDir + '/data/all_blast_results.tsv', delimiter = "\t")
    createCircos.runs(df, timeRecorder, RootDir, gbDir)
    recordAll.fin("All process", timeRecorder)
    print("\n\nFin.")
    pd.Series(timeRecorder).to_csv(RootDir + "Calculation_times.txt", header = None, index = None)
