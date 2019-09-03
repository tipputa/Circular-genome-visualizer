#! /usr/bin/env python
# coding:utf8

# from getEntry import GetEntry
# import searchAccesion

from Bio.Blast.Applications import NcbiblastxCommandline
import os



def createDB(fastaFileName, dbFileName, logName, dbType="prot"):
    toolPath = "singularity exec /usr/local/biotools/b/blast\:2.5.0--boost1.64_2 makeblastdb"
    command = " ".join(
        [toolPath, "-in", fastaFileName, "-dbtype", dbType, "-out", dbFileName, "-logfile", logName])
    os.system(command)


def blastp(queryFileName, dbFileName, outFileName):
    toolPath = 'singularity exec /usr/local/biotools/b/blast\:2.5.0--boost1.64_2 blastp'
    command = " ".join([toolPath, "-query", queryFileName, "-db", dbFileName, "-out", outFileName,
        "-evalue 0.00001", "-outfmt 5", "-num_alignments 1"])
    os.system(command)
