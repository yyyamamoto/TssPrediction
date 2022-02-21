# -*- coding: utf-8 -*-

#
"""
塩基の出現数を記録したファイルを２つ読み込み、IGIとPRI（またはFUI）を算出する。
PRIの算出時には、先に読み込まれたファイルのIGIから後に読み込まれたファイルのIGIが引かれる。
出現数はTSS_count.pyでの作成を想定している。
"""

from time import time
import os
from sys import argv
import re
import math

t1=time()
#開始時刻の取得

save_dir = os.path.dirname(__file__) + os.sep
#出力したファイルを保存するディレクトリ

def read_count_file(file_path):
    """
    塩基の出現数を数えたファイルを読み込み、塩基の組み合わせのリストと出現数のリストを返す
    """

    base_list, count_list = [], []

    with open(file_path, "r") as cf:
        for line in cf:
            if re.match("[ATGC]+\s\d+\n?", line):
                temp = line.replace('\n','').split('\t')
                base_list.append(temp[0])
                count_list.append(int(temp[1]))

    return base_list, count_list

def calc_igi(count_list):
    """
    塩基の出現数のリストを受け取り、算出したIGIのリストを返す
    """

    count_sum = sum(count_list)

    igi_list = []
    for c in count_list:
        if c == 0:
            igi_list.append(0)
        else:
            igi_list.append(math.log10(c / count_sum) - math.log10(1 / (4 ** 8)))
    return igi_list

def calc_pri(igi_list1, igi_list2):
    """
    ２つのIGIのリストを受け取り、算出したPRIのリストを返す
    """

    pri_list = []
    for igi1, igi2 in zip(igi_list1, igi_list2):
        if igi1 == 0 or igi2 == 0:
            pri_list.append(0)
        else:
            pri_list.append(igi1 - igi2)
    return pri_list


def saveWrite(save_file, *arg):
    """
    使用したファイル名を出力ファイルに書き込む
    """

    with open(save_file, "wt") as save:
        for item in arg:
            print(item, file=save)
        print("", file=save)


argvs=argv
#コマンドライン引数から処理するファイルを取得

if len(argvs) != 3:
    print('２つのファイルをコマンドライン引数で指定してください')
    quit()

igi_list_list = []
for count_file in argvs[1:]:
    base_list, count_list = read_count_file(count_file)

    igi_list = calc_igi(count_list)
    del count_list
    igi_list_list.append(igi_list)

    save_file = save_dir + "IGI_" + os.path.basename(count_file)
    #saveWrite(save_file, __file__, count_file)
    #with open(save_file, "a") as save:
    with open(save_file, "w") as save:
        for base, igi in zip(base_list, igi_list):
            print(base, igi, sep="\t", file=save)

pri_list = calc_pri(igi_list_list[0], igi_list_list[1])
del igi_list_list

save_file = save_dir + "PRI.txt"
#saveWrite(save_file, __file__, argvs[1], argvs[2])
#with open(save_file, "a") as save:
with open(save_file, "w") as save:
    for base, pri in zip(base_list, pri_list):
        print(base, pri, sep="\t", file=save)


t2=time()
t=t2-t1
if t > 60:
    print('time:'+str(t/60)+'(min)')
else:
    print('time:'+str(t)+'(s)')
#処理に要した時間を出力
