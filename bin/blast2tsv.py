
# from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import sys
import pandas as pd
import csv

def getIdentity(hsp):
    return hsp.identities * 100.0 / hsp.align_length


def getSimilarity(hsp):
    return hsp.positives * 100.0 / hsp.align_length


def getQueryCoverage(hsp, queryLength):
    return 100.0 * (hsp.query_end - hsp.query_start + 1) / queryLength


def getTargetCoverage(hsp, targetLength):
    return 100.0 * (hsp.sbjct_end - hsp.sbjct_start + 1) / targetLength

def formatFloat(f):
    return "%.1f" % (f)

def formatEval(Eval):
    return "0" if Eval==0 else "%.1e" % (Eval)

def showBestHit(record, revResult, tsDataFrame):
    queryName = getLocusTag(record.query)
    queryLength = record.query_length
    if len(record.alignments) > 0:
        bestAlignment = record.alignments[0]
        bestHSP = bestAlignment.hsps[0]
        targetLength = bestAlignment.length

        locusTag = getLocusTag(bestAlignment.hit_def)

        product, sequence, location = tuple(tsDataFrame.loc[locusTag,["product","sequence","location"]])

        if revResult[locusTag] == queryName:
            BBH = 1
        else:
            BBH = 0
        return [queryName, locusTag, product, sequence, location] + \
            [formatEval(bestHSP.expect), str(bestHSP.bits)] + \
            list(map(formatFloat,[getIdentity(bestHSP), getSimilarity(bestHSP), 
                getQueryCoverage(bestHSP, queryLength), getTargetCoverage(bestHSP, targetLength), BBH]))
    else:
        return [queryName,"-","-","-","-"] + ["-", "0", "0", "0", "0", "0", "0"]


def getLocusTag(definition):
    return definition.split()[0]


def getReverseHit(fileName):
    def getBestHit(record):
        bestAlignment = record.alignments[0]
        return getLocusTag(bestAlignment.hit_def)

    hitDict = {}
    result_handle = open(fileName)
    records = NCBIXML.parse(result_handle)
    for record in records:
        queryLT = getLocusTag(record.query)
        if len(record.alignments) > 0:
            targetLT = getBestHit(record)
            hitDict[queryLT] = targetLT
        else:
            hitDict[queryLT] = ""
    return hitDict



def main(forwardResult, reverseResult, targetSummary, outputFile):
    result_handle = open(forwardResult)
    records = NCBIXML.parse(result_handle)
    revResult = getReverseHit(reverseResult)
    tsDataFrame=pd.read_csv(targetSummary, sep="\t", index_col="LocusTag")

    Buffer=[]
    Buffer.append(["query", "target","product","sequence", "location","Eval", "score", "identity", "similarity", "query_coverage","target_coverage", "BBH"])

    for record in records:

        Buffer.append(showBestHit(record, revResult, tsDataFrame))

    with open(outputFile,"w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(Buffer)

if __name__ == '__main__':
    import os.path
    ROOT = sys.argv[1]
    PROJECT = sys.argv[2]
    SELF = sys.argv[3]
    OTHER =sys.argv[4]

    # run   python blastParser.py /home/ytanizaw/mica LH LOOC260 vini

    forwardResult = os.path.join(ROOT,"projects", PROJECT, "data", SELF, "BlastResult", OTHER + ".xml")
    reverseResult = os.path.join(ROOT,"projects", PROJECT, "data", OTHER, "BlastResult", SELF + ".xml")
    targetSummary = os.path.join(ROOT,"projects", PROJECT, "data", OTHER, "summary.tsv")
    outputFile = os.path.join(ROOT,"projects", PROJECT, "data", SELF, "BlastResult", OTHER + ".tsv")

    print(forwardResult, reverseResult, targetSummary, outputFile)
    main(forwardResult, reverseResult, targetSummary, outputFile)


