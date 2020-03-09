# -*- coding: utf-8 -*-

#
"""
Author: Tosei Hiratsuka (2020.3.9)

READ manual.txt BEFORE USE.

Input: FASTA, IGI, PRI
Output: The numeric sequence converted by IGI or PRI

At the command prompt, enter the following.
    python3 chrom_scan.py [FASTA file] [Score file 1] [Score file 2]...

For example,
    python3 chrom_scan.py Chr1_for_test.con IGI200_60_forChr1-2.txt IGI750_450_forChr1-2.txt PRI200_60-750_450_forChr1-2.txt 5UTR-IGI750_450_forChr1-2.txt


If a file with the same name as the output file already exists, it will be overwritten without warning.
If you want to change the output directory, enter the variable ("save_dir") below.
"""

save_dir = r""
#The directory for saving output files
#If blank, the output file will be in the directory where this script is located.



import os
import time
import sys
import re

t1=time.time()
#Get start time.

if save_dir == "":
    save_dir = os.path.dirname(os.path.abspath(__file__)) + os.sep
    #Get the absolute path of the directory where this script is located.


def score_name_extract(score_table_file):
    """
    To extract the name of score from the IGI or PRI file.
    :param score_table_file: The IGI or PRI file.
    :return: The name of score.
    """

    fn=os.path.basename(score_table_file)
    #score_name = re.search("([IGI|PRI].+)\.txt", fn).group(1)
    score_name = re.search("([IGI|PRI|5UTR].+)\_for", fn).group(1)
    return score_name

def chrom_name_extract(chrom_file):
    """
    To extract the chromosome number from the FASTA file.
    :param chrom_file: FASTA file.
    :return: The chromosome number.
    """
    
    fn = os.path.basename(chrom_file)
    chrom = re.search("([cC]hr\d+)\D.*", fn).group(1)
    return chrom

def save_generate(save_directory, chrom, score_name):
    """
    To generate output files (forward and reverse).
    :param save_directory: The directory for saving output files.
    :param chrom: The chromosome number.
    :param score_name: The name of score (IGI or PRI).
    :return: Names of output files.
    """
    
    save_f = save_directory + "{}_scan_{}_forward.txt".format(chrom, score_name)
    save_r = save_f.replace("forward", "reverse")
    return save_f, save_r

def Reverse(sequence):
    """
    To generate the complementary sequence.
    :param sequence: The base sequence.
    :return: The complementary sequence to the input sequence.
    """
    
    basedict={'A':'T','T':'A','G':'C','C':'G'}
    #塩基対のディクショナリ
    rsequence=''
    for base in sequence[::-1]:
        rsequence+=basedict.get(base,'X')
    #The tail of "sequence" is equal to the head of "rsequence".
    #This function outputs "X" for mixed base.
    return rsequence

def index_dict_generate(score_table_file):
    """
    To extract the dictionary and the number of bases from PRI or IGI file.
    :param score_table_file: IGI or PRI file.
    :return: The dictionary object (key: base sequence, value: IGI or PRI score).
    :return: The number of bases.
    """
    
    with open(score_table_file,'rt') as fi:
        indexlist=[line.replace('\n','').split('\t') for line in fi
                   if re.match("[ATGC][ATGC]+", line)]
        #If the base sequence is in the line, generate the list object separating by a tab character.

    bp=len(indexlist[0][0])
    #The number of bases.

    indexdict=dict(indexlist)
    #{"base sequence": "IGI or PRI score"}
    del indexlist

    return indexdict, bp

def seq_ext(sequence_file):
    """
    To extract the whole of sequence from the FASTA file.
    :param sequence_file: FASTA file.
    :return: The whole of sequence.
    """
    
    seq = ""
    with open(sequence_file, "rt") as sf:
        for line in sf:
            if not re.match(">", line):
                seq += line.replace("\n", "")
    return seq

def saveWrite(save_file, chrom, score, strand, *arg):
    """
    To write basic information to the output file.
    :param save_file: The output file.
    :param chrom: Chromosome number.
    :param score: The name of IGI or PRI score.
    :param strand: The chain.（forward or reverse）
    :param *arg: Others that you want to write to the output file. These are written at the beginning of the file.
    """
    
    with open(save_file, "wt") as save:
        for item in arg:
            print(item, file=save)
        print("\nChromosome:", chrom, file=save)
        print("Score:", score, file=save)
        print("Strand:", strand, file=save)
        print("\nPosition", "Score", sep="\t", file=save)


argvs = sys.argv
if len(argvs) < 3:
    print('Input one chromosome base sequence file and one or more PRI or IGI files with command line arguments.')
    quit()
    #If no more than two files are entered, then exit.

scoreFiles = []

for f in argvs[1:]:
    if re.search("PRI", os.path.basename(f)) or re.search("IGI", os.path.basename(f)):
        scoreFiles.append(f)
    else:
        seqFile = f

for scfile in scoreFiles:
    scname = score_name_extract(scfile)
    scDict, splen = index_dict_generate(scfile)
    #"splen" is the width to separate array when scoring. (split length)
    if splen % 2 == 0:
        offset = int(splen / 2)
    else:
        offset = int((splen +1 ) / 2)

    print("--------------------------------")
    print(os.path.basename(scfile))
    print(os.path.basename(seqFile))

    chrom = chrom_name_extract(seqFile)
    saveF, saveR = save_generate(save_dir, chrom, scname)
    saveWrite(saveF, chrom, scname, "forward", __file__, seqFile, scfile)
    saveWrite(saveR, chrom, scname, "reverse", __file__, seqFile, scfile)

    seq = seq_ext(seqFile)

    F_lst = []
    R_lst = []
    #To store scores of sequence.
    for i in range(len(seq) - splen + 1):
        part = seq[i : i + splen]
        F_lst.append(scDict.get(part,"0"))

        rpart = Reverse(part)
        R_lst.append(scDict.get(rpart, "0"))
    del seq, scDict

    with open(saveF, "at") as fsave:
        for pos, fresult in enumerate(F_lst, offset):
            print(pos, fresult, sep="\t", file=fsave)
    del F_lst

    with open(saveR, "at") as rsave:
        for pos, rresult in enumerate(R_lst, offset):
            print(pos, rresult, sep="\t", file=rsave)
    del R_lst


t2=time.time()
t=t2-t1
if t > 60:
    print('time:'+str(t/60)+'(min)')
else:
    print('time:'+str(t)+'(s)')
#To display the time required for processing.