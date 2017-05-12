#! /usr/bin/env python
# coding:utf8

# from getEntry import GetEntry
# import searchAccesion

from Bio.Blast.Applications import NcbiblastxCommandline
import os



def createDB(fastaFileName, dbFileName, logName, dbType="prot"):
    toolPath = "makeblastdb"
    command = " ".join(
        [toolPath, "-in", fastaFileName, "-dbtype", dbType, "-out", dbFileName, "-logfile", logName])
    os.system(command)


def blastp(queryFileName, dbFileName, outFileName):
    comline = NcbiblastxCommandline(query=queryFileName, db=dbFileName, out=outFileName,
                                    cmd="blastp", evalue=0.001, outfmt=5, num_alignments=1)
    stdout, stderr = comline()
    return stdout, stderr


