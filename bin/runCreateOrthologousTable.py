# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 21:21:31 2019

@author: tipputa
"""

import pandas as pd
import sys, os
import createOrthologousGenes as createOrtho
import timeRecord as recorder

useage = """

Error: please input the absolute PATH of output and genbank directories .

Usage: 
   python runCreateOrthologousTable.py <output directory> <input directory (GenBank files)> <optional; True = after blast>
   e.g: python runCreateOrthologousTable.py ~/study/ ~/study/gb/
"""


if __name__ == '__main__':
    timeRecorder = []
    recordAll = recorder.timeRecord()
    afterBlast = False
    
    if len(sys.argv)==4:
        RootDir = sys.argv[1] + "/"
        gbDir = sys.argv[2] + "/"
        afterBlast = True
  
        
    elif len(sys.argv)==3:
        RootDir = sys.argv[1] + "/"
        gbDir = sys.argv[2] + "/"
        
    else:        
        print(useage)
        quit()
    
    os.chdir(RootDir)
    if (afterBlast):
        createOrtho.clustering_woDuplicated(timeRecorder, RootDir, gbDir)
    else:
        createOrtho.clustering_woDuplicated_afterBlast(timeRecorder, RootDir, gbDir)
