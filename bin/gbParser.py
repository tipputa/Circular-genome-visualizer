#! /bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 11:20:44 2017

@author: tipputa
"""
import csv

def locationToStr(locationObj):
    start = locationObj.start + 1
    end = locationObj.end
    before = "<" if "<" in str(locationObj.start) else ""
    after = ">" if ">" in str(locationObj.end) else ""

    if locationObj.strand == 1:
        return "%s%d..%s%d" % (before, start, after, end)
    elif locationObj.strand == -1:
        return "complement(%s%d..%s%d)" % (before, start, after, end)
    
def parseFeature(feature, num):
    locus_tag = feature.qualifiers.get("locus_tag", ["locus_" + str(num)])[0]
    location = locationToStr(feature.location)
    product = feature.qualifiers.get("product", ["product_" + str(num)])[0]
    function = feature.qualifiers.get("function", [""])[0]
    translation = feature.qualifiers.get("translation", [""])[0]
    return locus_tag, location, product, function, translation
                
def createProteinFasta(record, proteinFileName):
    num = 0
    with open(proteinFileName, "w") as f:
        for feature in record.features:
            num += 1
            if "translation" in feature.qualifiers.keys():
                locus_tag, location, product, function, translation = parseFeature(
                    feature, num)
                f.write(">%s %s %s\n" %
                    (locus_tag, product, record.name + ":" + location))
                f.write(feature.qualifiers["translation"][0] + "\n")
                
def createTSV(record, tsvFileName):
    csvWriter = csv.writer(
        open(tsvFileName, "w"), lineterminator="\n", delimiter='\t')
    Buffer = [("LocusTag", "sequence", "location",  "feature",
               "product", "function", "translation")]
    num = 0
    for feature in record.features:
        num += 1
        if feature.type in ["CDS", "rRNA", "tRNA"]:
            locus_tag, location, product, function, translation = parseFeature(feature, num)
            Buffer.append((locus_tag, record.name, location,
                               feature.type, product, function, translation))
    csvWriter.writerows(Buffer)
