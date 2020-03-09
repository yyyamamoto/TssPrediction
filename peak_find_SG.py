# -*- coding: utf-8 -*-
#
"""
Author: Tosei Hiratsuka (2020.3.9)

READ manual.txt BEFORE USE.

Input: The scoring data of the base sequence
Output: The position and value of the detected peak

At the command prompt, enter the following.
    python3 peak_find_SG.py [Scored by IGI1] [Scored by IGI2] [Scored by PRI]...

For example,
    python3 peak_find_SG.py Chr1_scan_IGI200_60_forward.txt Chr1_scan_IGI750_450_forward.txt Chr1_scan_PRI200_60-750_450_forward.txt


The peak is filtered by
・ Derivation of SG filter
・ PRI threshold
・ IGI threshold
・ Integration of peaks at close distances

The scored file is assumed to be generated by chrom_scan.py.


If a file with the same name as the output file already exists, it will be overwritten without warning.
If you want to change the output directory, enter the variable ("save_dir") below.
"""

save_dir = r""
#The directory for saving output files
#If blank, the output file will be in the directory where this script is located.



import os
import sys
import re
import numpy as np 
from scipy import signal
import time

t1=time.time()
#Get start time.

if save_dir == "":
    save_dir = os.path.dirname(os.path.abspath(__file__)) + os.sep
    #Get the absolute path of the directory where this script is located.

N = 151
#Smoothing width (odd number)
D = 1
#Order used for smoothing by SG filter.

th_PRI = 0
th_IGI = 0
#The threshold for peak detection.

tougou = 100
#The width for peak integration (bp).

def infoExt(scanFile):
    """
    To extract information from the file scored by chrom_scan.py.
    :param scanFile: The file scored by chrom_scan.py.
    :return: Chromosome number, Direction, Score name.
    """

    with open(scanFile, "rt") as sf:
        for line in sf:
            if re.match("Chromosome", line):
                chrom = re.search("([cC]hr\d+)\D?", line).group(1)
            elif re.match("Score", line):
                scName = re.search("Score:\s(.+)\n", line).group(1)
            elif re.match("Strand", line):
                strand = re.search("Strand:\s(.+)\n", line).group(1)
            elif re.match("\d+", line):
                break
    return chrom, scName, strand

def scoreExt(scanFile):
    """
    To generate an array of scores from the file scored by chrom_scan.py.
    Also returns the top location.
    :param scanFile: The file scored by chrom_scan.py.
    :return: numpy.ndarray of locations and values.
    """

    start = ""
    score = []
    with open(scanFile, "rt") as sf:
        for line in sf:
            if re.match("\d", line):
                line = line.replace("\n", "").split("\t")
                score.append(float(line[1]))
                if start == "":
                    start = int(line[0])
    score = np.array(score)
    return start, score

def sgfilter(ndarray, window, degree):
    """
    Apply an SG-filter to the input data and return the smoothed value.
    :param ndarray: numpy.ndarray of raw data.
    :param window: Smoothing width (odd number)
    :param degree: Order used for smoothing
    :return: numpy.ndarray of smoothed data.
    """

    y = ndarray
    #raw data
    ys = signal.savgol_filter(y, window, degree)
    #smoothed data
    return ys

def peakFind(rawPRI, IGI1, IGI2, window, degree, threshold_PRI, threshold_IGI, start):
    """
    Apply an SG-filter to the input data and return the position of the peak
    :param rawPRI: Raw data of PRI value (numpy.ndarray)
    :param IGI1: Smoothed IGI (numpy.ndarray) Part 1
    :param IGI2: Smoothed IGI (numpy.ndarray) Part 2
    :param window: Smoothing width (odd number)
    :param degree: Order used for smoothing
    :param threshold_PRI: PRI threshold at peak detection
    :param threshold_IGI: IGI threshold at peak detection
    :param start: First PRI location (return value of scoreExt)
    :return: numpy.ndarray of peak positions and values
    """

    ps = signal.savgol_filter(rawPRI, window, degree)
    #smoothed data
    pd = signal.savgol_filter(rawPRI, window, degree, deriv=1)
    #First derivative of smoothed data
    del rawPRI

    q = int((window - 1) / 2)
    #Peak search range (one side)

    peakI = []
    #To store indexes of peaks.
    for i in range(q*2, ps.size - 1 - q*2 + 1):
        #Search for peaks only in indexes with valid smoothing values.
        if pd[i] >= 0 and pd[i+1] < 0:
            #"i" is the provisional position of peak.
            peak_index = ps[i-q : i+q+1].argmax() + i - q
            #Put the maximum value in the smoothed range as the peak.
            #The index of peak in original array.
            if ps[peak_index] > threshold_PRI and IGI1[peak_index] > threshold_IGI and IGI2[peak_index] > threshold_IGI:
                if len(peakI) == 0 or peak_index != peakI[-1]:
                    peakI.append(peak_index)
                    #If the index of peak does not overlap with the previous one, add it to the list.
    del IGI1, IGI2
    peak = ps[peakI]
    del ps
    peakPos = np.array(peakI) + start
    return peakPos, peak

def peakUnite(peak_x, peak_y, window):
    """
    Input: peak position
    Output: merges and returns peaks at close distances
    :param peak_x: Peak position (generated by peakFind)
    :param peak_y: Peak value (generated by peakFind)
    :param window: Range to combine peaks
    :return: Integrated peak position (numpy.ndarray)
    """

    resultX = []
    resultY = []
    #To store integrated peak positions.

    for pos in peak_x:
        posUp = pos - window
        posDown = pos + window
        partX = peak_x[(posUp <= peak_x) & (peak_x <= posDown)]
        partY = peak_y[np.where((posUp <= peak_x) & (peak_x <= posDown))]
        #Extract peaks within a distance of "window" bases from peak.
        
        if partX.size > 0:
            peak = partY.max()
            peakX = partX[partY.argmax()]
            #Choose only the maximum value of peaks within the range.
            if not peakX in resultX:
                resultX.append(peakX)
                resultY.append(peak)
                #When there is no overlap, add the peak to the list.
    return np.array(resultX), np.array(resultY)

def save_generate(save_directory, chrom, orientation):
    """
    To generate output files (forward and reverse).
    :param save_directory: The directory for saving output files.
    :param chrom: The chromosome number.
    :param orientation: forward or reverse.
    :return: Names of output files.
    """
    
    saveFile = save_directory + "{}_peak_{}.txt".format(chrom, orientation)
    return saveFile

def saveWrite(save_file, chrom, strand, window, degree, threshold_PRI, threshold_IGI, integration, *arg):
    """
    To write basic information to the output file.
    :param save_file: The output file.
    :param chrom: Chromosome number.
    :param strand: The chain.（forward or reverse）
    :param window: The width for smoothing.
    :param degree: Order used for smoothing
    :param threshold_PRI: PRI threshold at peak detection
    :param threshold_IGI: IGI threshold at peak detection
    :param *arg: Others that you want to write to the output file. These are written at the beginning of the file.
    """
    
    with open(save_file, "wt") as save:
        for item in arg:
            print(item, file=save)
        print("\nChromosome:", chrom, file=save)
        print("Score:", score, file=save)
        print("Strand:", strand, file=save)
        print("Window Size:", window, file=save)
        print("Degree:", degree, file=save)
        print("PRI threshold: {}\nIGI threshold: {}\nPeak Integration: {}".format(threshold_PRI, threshold_IGI, integration), file=save)
        print("\nPosition", "Score", sep="\t", file=save)

inputs = sys.argv
if len(inputs) < 4 or (len(inputs) - 1) % 3 != 0:
    print('Input scored chromosome file, (one PRI file, two IGI files) * n as command line arguments.')
    quit()

PRI_F = []
PRI_R = []
IGI_F = []
IGI_R = []
for f in inputs[1:]:
    fn = os.path.basename(f)
    if "PRI" in fn:
        if "forward" in fn:
            PRI_F.append(f)
        else:
            PRI_R.append(f)
    elif "IGI" in fn:
        if "forward" in fn:
            IGI_F.append(f)
        else:
            IGI_R.append(f)

flst = []
for i in range(len(PRI_F)):
    flst.append(PRI_F[i])
    flst.append(IGI_F[i*2])
    flst.append(IGI_F[i*2+1])
for i in range(len(PRI_R)):
    flst.append(PRI_R[i])
    flst.append(IGI_R[i*2])
    flst.append(IGI_R[i*2+1])

del inputs, PRI_F, PRI_R, IGI_F, IGI_R

for i in range(0, len(flst), 3):
    pri = flst[i]
    igi1 = flst[i+1]
    igi2 = flst[i+2]
    print("-"*50)
    print(os.path.basename(pri))
    print(os.path.basename(igi1))
    print(os.path.basename(igi2))
    
    chrom, score, ori = infoExt(pri)
    save = save_generate(save_dir, chrom, ori)
    saveWrite(save, chrom, ori, N, D, th_PRI, th_IGI, tougou, __file__, pri, igi1, igi2)

    start, PRI = scoreExt(pri)
    start, IGI1 = scoreExt(igi1)
    start, IGI2 = scoreExt(igi2)
    print('"scoreExt complete"')

    IGI1 = sgfilter(IGI1, N, D)
    IGI2 = sgfilter(IGI2, N, D)
    print('"S-G filter complete"')

    peakPos, peak = peakFind(PRI, IGI1, IGI2, N, D, th_PRI, th_IGI, start)
    del PRI, IGI1, IGI2
    print('"peakFind complete"')

    peakPos, peak = peakUnite(peakPos, peak, N)
    print('"peakUnite complete"')

    for pos, p in zip(peakPos, peak):
        with open(save, "at") as s:
            print(pos, p, sep="\t", file=s)

t2=time.time()
t=t2-t1
if t > 60:
    print('time:'+str(t/60)+'(min)')
else:
    print('time:'+str(t)+'(s)')
#To display the time required for processing.