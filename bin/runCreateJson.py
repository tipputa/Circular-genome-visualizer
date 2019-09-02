# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 19:00:27 2017

@author: ac144809
"""

import os, sys
import createJson

useage="""
python runCreateJson.py  <output directory> <gbDir>
"""

if __name__ == '__main__':
    
    args = len(sys.argv)
    RootDir = sys.argv[1]
    gbDir = sys.argv[2]

    os.system("cd " + RootDir)
    createJson.runs(RootDir, gbDir)
