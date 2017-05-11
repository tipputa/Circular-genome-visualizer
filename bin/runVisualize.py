# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 19:00:27 2017

@author: ac144809
"""

import os, sys
import visualizeUtil as visualizer
  
useage="""
python runVisualize.py  <output directory> <data info file> <option; key word for output; default:"test"> <option; minimum number of genes in each cluster> <option; sorting column name>
"""

if __name__ == '__main__':
    
    args = len(sys.argv)
    if args == 6:
        RootDir = sys.argv[1] + "/"
        df_name = sys.argv[2]
        tag = sys.argv[3]
        min_genes = sys.argv[4]
        sort_key = sys.argv[5]
    elif args == 5:
        RootDir = sys.argv[1] + "/"
        df_name = sys.argv[2]
        tag = sys.argv[3]
        min_genes = [4]
        sort_key = None
    elif args == 4:
        RootDir = sys.argv[1] + "/"
        df_name = sys.argv[2]
        tag = sys.argv[3]
        sort_key = None
        min_genes = None
    elif args == 3:
        RootDir = sys.argv[1] + "/"
        df_name = sys.argv[2]
        tag = "test"
        sort_key = None
        min_genes = None
    else:
        print(useage)
        quit()
    
    CircosIN = []
    df_info, df_locusTag, Gsize = visualizer.readInfo(RootDir, df_name, CircosIN, min_genes)
    os.chdir(RootDir + "circos/")
    visualizer.createCircosInput(RootDir, df_locusTag, Gsize, CircosIN, df_info, tag)
    visualizer.createCircosConf(RootDir, df_info,sort_key, tag)
    print("Run circos :"  + tag )
    os.system('circos')
    cmd = "mv circos.png ../circos_" + tag + ".png"     
    os.system(cmd)


