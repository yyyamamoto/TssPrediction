# -*- coding: utf-8 -*-

#
"""
塩基配列を端から８塩基ごとに区切り、各塩基の組み合わせの出現回数を数える
入力する塩基配列は-1000から始まるものでなければならない
塩基配列はファイルの各行の末尾になければならない
"""

#----------------------------設定ここから----------------------------
region_list=[ [-200, -60], [-750, -450] ]
#探索領域のリスト（TSSより上流に限る）
#----------------------------設定ここまで----------------------------


import os
save_dir = os.path.dirname(__file__) + os.sep
#出力したファイルを保存するディレクトリ

combination_file = os.path.dirname(__file__) + os.sep + r"8bp_all.txt"
#区切る塩基の全ての組み合わせを書いたファイル

bp=8
#区切る塩基数の指定


from os.path import basename
from sys import argv
import re

def chrom_ext(promoter_file):
    """
    ファイル名から染色体番号を返す
    :param file: ファイル名(染色体番号がfor_chr1のような形式で記されているものを想定）
    :return: 染色体番号(例:chr1)
    """
    fname=basename(promoter_file)

    temp = re.search("[cC]hr(\d+)\D", fname)
    if temp != None:
        chrom = temp.group(1)
        return "chr" + chrom
    else:
        return ""

def tss_count(promoter_file,combination_file,bp,start,end, save_directory):
    """
    塩基の組み合わせがプロモーターの特定の領域に何回ずつ出現するかをファイルに出力する
    :param promoter_file: プロモーター配列のファイル
    :param combination_file: bp個の塩基の組み合わせを記したファイル
    :param bp: 塩基数
    :param start: 探索領域の始点（上流側）
    :param end: 探索領域の終点（下流側）
    :param save_directory: 出力ファイルの保存場所
    """

    fn_promoter=promoter_file
    fn_8bp=combination_file

    chrom=chrom_ext(fn_promoter)
    fn_save = save_directory + r"for_" + chrom+"_"+str(start)+'_'+str(end)+r'.txt'
    #raw文字列でファイルを指定

    save=open(fn_save,'w')

    print(__file__, file=save)
    print(fn_promoter, file=save)
    print(fn_8bp, file=save)
    save.write(str(start)+':'+str(end)+'\n')
    #使用したファイル名を書き込む


    f8=open(fn_8bp,'r')

    baseslist = []
    for line in f8:
        bases = line.replace("\n", "")
        baseslist.append(bases)


    promlist=[]
    #空のリストを作成

    fp=open(fn_promoter,'r')

    for line in fp:
        line = line.replace('\n','')
        line = filter(lambda w: len(w) > 0, re.split(r'\s', line))
        prom = list(line)[-1]
        promlist.append(prom)
    #プロモーターファイルから1行ずつ読み込み
    #改行文字の削除
    #行の末尾の配列部分を抽出
    #promolistの末尾に追加

    fp.close()


    seqlist=[prom[1000+start:(1000+start)+(-1*start+end)+1] for prom in promlist]
    #プロモーターの配列から探索領域を抽出し、seqlistに追加

    partlist=[seq[i:i+bp] for seq in seqlist for i in range(len(seqlist[0])-bp+1)]
    #seqlistの要素について、先頭から８塩基ずつ抜き出してリスト化

    del promlist
    del seqlist


    num={}
    #空のディクショナリを作成

    for part in partlist:
        num[part]=num.get(part,0)+1
    #キーに配列、要素に配列の数を追加
    #キーが存在しなければ新たに追加され、要素は0+１となる
    #既にキーが存在すれば、既存の要素に１足した値が要素に代入される

    for bases in baseslist:
        save.write(bases + '\t' + str(num.get(bases,0)) +'\n')
    #塩基の組み合わせリストの順で数を出力
    #ディクショナリは並び順が一定でないため

    save.close()

argvs=argv
#コマンドライン引数から処理するファイルを取得

if len(argvs)<2:
    print(' -1000から始まるプロモーター配列ファイルをコマンドライン引数で指定してください')
    quit()

print("="*50)
print("split:", bp)
for pfile in argvs[1:]:
    print("-"*50)
    print(basename(pfile))
    for item in region_list:
        start=item[0]
        end=item[1]
        tss_count(pfile,combination_file,bp,start,end, save_dir)
        print(start,":",end,"completed")
print("="*50)
