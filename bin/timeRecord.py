#! /bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 11:20:44 2017

@author: tipputa
"""
import time

class timeRecord():
    start = time.time()
    def fin(self,str,timeRecorder):
        elapsed_time = round(time.time() - self.start,1)
        if elapsed_time > 60:
            elapsed_time_min = round(elapsed_time/60,1)
            timeRecordS = "Elapsed time (" + str + "): {0}".format(elapsed_time_min) + " [min]"
        else:
            timeRecordS = "Elapsed time (" + str + "): {0}".format(elapsed_time) + " [sec]"

        print("\n" + timeRecordS + "\n")
        timeRecorder.append(timeRecordS)
        self.start = time.time()
