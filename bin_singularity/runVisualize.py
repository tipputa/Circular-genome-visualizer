# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 19:00:27 2017

@author: ac144809
"""

import os, sys
import visualizeUtil as visualizer

useage="""
python runVisualize.py  <output directory> <configuration file> <option; key word for output; default:"test"> <option; the minimum number of genes in each cluster; default: 1> <option; sorting column name; default: None>
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
        min_genes = sys.argv[4]
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

    os.system("cd " + RootDir)
    CircosIN = []
    df_info, df_locusTag, consensus_gsize = visualizer.readInfo(RootDir, df_name, CircosIN, min_genes)

    def visualize(dfs, sort_key, tag):        
        print("\nCreate Circos input (" + RootDir + "./circos/data/circos_" + tag + "_*.txt)")
        visualizer.createCircosInput(RootDir, dfs, consensus_gsize, CircosIN, df_info, tag)
        visualizer.createCircosConf(RootDir, df_info, sort_key, tag)
        print("circos suffix: "  + tag )
        os.system('singularity exec /usr/local/biotools/c/circos\:0.69.2--0 circos -conf '+ RootDir +'circos/circos.conf')
        os.system("mv circos.png " + RootDir + "circos_" + tag + ".png")
        os.system("rm circos.svg")
        print("Output circos_" + tag + ".png, circos_" + tag + "_info.txt, and RingOrder_" + tag + "_df.tsv\n")

    visualize(df_locusTag, sort_key, tag)
