#from Bio.Blast import NCBIWWW

import sys, os, os.path
import pandas as pd
#import numpy as np
#import pathUtil




class BlastResult:

    def __init__(self, query, locusTag, product, sequence, location, Eval, score, identity, similarity, query_coverage, target_coverage, BBH):
        self.query = query
        self.locusTag = locusTag
        self.product = product
        self.sequence = sequence
        self.location = location
        self.Eval = float(Eval) if locusTag != "-" else "-"
        self.score = score
        self.identity = identity
        self.similarity = similarity
        self.query_coverage = query_coverage
        self.target_coverage = target_coverage
        self.BBH = BBH
        self.visible = True if locusTag != "-" else False

    def summary(self):
        if self.visible:
            return "%s\n%s\nScore : %.1f\nE-val : %.1e\nIdentity : %.1f\nQ-Coverage : %.1f\nT-Coverage : %.1f" % \
                (self.locusTag, self.product, self.score, self.Eval, self.identity, self.query_coverage, self.target_coverage)
        else:
            return "-"

def getResultDF(tsvFileName, targetName):

    D=pd.read_csv(tsvFileName, sep="\t")
    D["result"]=[BlastResult(*tuple(D.loc[Index])) for Index in D.index]
    D["strain"]=targetName
    return D[["query","strain","result"]]





def getSummaryDataFrame(tsvFileName):

    df = pd.read_csv(tsvFileName, index_col=0,sep="\t")
    return df[["sequence","location","feature","product"]]

def createTable(dir, pivot, targetList):
    targets = targetList.split(",")
    D = pd.DataFrame(None, columns=["query","strain","result"])
    for target in targets:
        resultTsv = dir + pivot + "vs" + target + ".tsv"
        d = getResultDF(resultTsv, target)
        D = D.append(d)
    DP = D.pivot(index="query",columns="strain")
    DP.columns=[x[1] for x in DP.columns]
    DP = DP[targets]
    return DP



# below are metrics functions
def Zvalue(dfRow):
    '''For standardization by row'''

    S=pd.Series([x.score for x in dfRow if x.visible])
    STDEV=S.std(ddof=0)
    MEAN=S.mean()

    if len(S) == 1 or STDEV==0:
#         dfRow[dfRow > 0] = 0
        return [0 if x.visible else "-" for x in dfRow ]
    else:
        return [(x.score -MEAN)/STDEV if x.visible else "-" for x in dfRow ]

def ratio(dfRow):
    '''score ratio by row : score / Max(score)'''

    S=pd.Series([x.score for x in dfRow if x.visible])
    if len(S)>0:
        MAX=S.max()
    else:
        MAX=1
    return [x.score/MAX if x.visible else "-" for x in dfRow ]

def LOCUSTAG(dfRow):
    '''show locus tag'''

    return [x.locusTag if x.visible else "-" for x in dfRow ]

def getScore(blastResult):
    return blastResult.score if blastResult.visible else "-"

def getLT(blastResult):
    return blastResult.locusTag if blastResult.visible else "-"

def getSummary(blastResult):
    return blastResult.summary()

#below are filter function
# def isNotBBH(blastResult):
#     return blastResult.BBH == 0



def createFilterFunction(score=0, identity=0, similarity=0, Qcoverage=0, Tcoverage=0, BBH=False):
    def _filterFunc(blastResult):
        if blastResult.visible:
            if BBH is True and blastResult.BBH == 0:
                blastResult.visible = False
                return
            if blastResult.query_coverage < Qcoverage:
                blastResult.visible = False
                return
            if blastResult.target_coverage < Tcoverage:
                blastResult.visible = False
                return
            if blastResult.score < score:
                blastResult.visible = False
                return
            if blastResult.identity < identity:
                blastResult.visible = False
                return
            if blastResult.similarity < similarity:
                blastResult.visible = False
                return
    return _filterFunc

if __name__ == '__main__':

     #if len(sys.argv)>=2:
     project = sys.argv[1]
     pivot = sys.argv[2]
     targetList = sys.argv[3]
     name = pivot
    #else:
	#project="TO"
	#pivot="1002"
	#targetList="1000,1001,1002,1003,16,2165,4_3,AG30,ATCC14917,AY01,DmCS_001,EGD-AQ4,FMNP01,JDM1,Lp90,NC8,P-8,ST-III,UCMA3037,WCFS1,WJL,ZJ316"
	#name="1002"
    # Change the line below to change filter threshold
     filterFunc = createFilterFunction(score=0, identity=30, similarity=0, Qcoverage=0, Tcoverage=0, BBH=True)
    # filterFunc = createFilterFunction(score=0, identity=30, similarity=0, Qcoverage=60, Tcoverage=60, BBH=False)
    # filterFunc = createFilterFunction(score=0, identity=30, similarity=0, Qcoverage=0, Tcoverage=0, BBH=False)
     D = createTable(project, pivot, targetList)

     summaryTSVFile = os.path.join(pathUtil.getGenomeDataDir(project, pivot), "summary.tsv")
     DS = getSummaryDataFrame(summaryTSVFile)

     analysisPath=pathUtil.getProjectRoot(project) + "/analysis/"
     if not os.path.exists(analysisPath):
         os.makedirs(analysisPath)


     D.applymap(filterFunc)

 #    DS.join(D.applymap(getSummary), how='inner').to_csv(analysisPath + name + "_table_data.tsv", sep="\t", na_rep='-',index_label="locusTag")
 #    DS.join(D.apply(Zvalue, axis=1), how='inner').to_csv(analysisPath + name + "_table_Zvalue.tsv", sep="\t", na_rep='-',index_label="locusTag")
 #    DS.join(D.apply(ratio, axis=1), how='inner').to_csv(analysisPath + name + "_table_ratio.tsv", sep="\t", na_rep='-',index_label="locusTag")
     DS.join(D.apply(LOCUSTAG, axis=1), how='inner').to_csv(analysisPath + name + "_table_locusTag.tsv", sep="\t", na_rep='-',index_label="locusTag")


