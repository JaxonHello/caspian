import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import heapq
import argparse
import math
import matplotlib as mpl
import hdbscan
import pandas as pd
from sklearn import metrics
import datetime
###########################################
#modify the parametr:metricList,dir,res,hicname
metricList=["euclidean","manhattan","chebyshev"]
dir="TAD_results/SimulationData"
res = 40000
hicname="4noise"
hicfile=dir+"/"+hicname+".hic"

mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42
def boundaryPlot(labels):
    n = len(labels)
    boundary = np.zeros(n)
    i = 0
    label = -1
    start = 0
    while i < n:
        if labels[i] == label:
            boundary[i] = start
        else:
            start = i
            label = labels[i]
            boundary[i] = i
        i = i + 1
    return boundary

def readTAD(tadfile):
    #tads = "/home/ghaiyan/project/CASPIAN/evaluate_TADS/GM12878/chr19_5kb/TAD/{}.txt".format(tadsname)
    f = open(tadfile)
    line=f.readline()
    start=[]
    end=[]
    while line:
        line = line.split()
        start1 = int(line[0])
        end1 = int(line[2])
        start.append(start1)
        end.append(end1)
        line=f.readline()
    f.close()
    return start, end

def tadQuality(tadFile,hic):
    """TAD quality"""
    n = len(hic)
    tad = np.loadtxt(tadFile)
    intra = 0
    intra_num = 0
    for n in range(len(tad)):
        for i in range(int(tad[n,0]),int(tad[n,2]+1)):
            for j in range(int(tad[n,0]),int(tad[n,2]+1)):
                intra = intra + hic[i,j]
                intra_num = intra_num + 1

    if intra_num!=0:
        intra = intra / intra_num
        print("intra TAD: %0.3f" % intra)
    else:
        intra = 0
    
    inter = 0
    inter_num = 0
    for n in range(len(tad) - 1):
        for i in range(int(tad[n,0]),int(tad[n,2]+1)):
            for j in range(int(tad[n+1,0]),int(tad[n+1,2]+1)):
                inter = inter + hic[i,j]
                inter_num = inter_num + 1
    if inter_num!=0:
        inter = inter / inter_num
        print("inter TAD: %0.3f" % inter)
    else:
        inter = 0
    print("quality: %0.3f" % (intra - inter))
    quality=(intra - inter)/len(tad)
    return quality

def getLabel(hic,start, end):
    n = len(hic)
    labels = np.zeros(n)
    for i in range(n):
        labels[i] = 0
    for j in range(len(start)+1):
        s=start[j-1]
        m=end[j-1]
        labels[s]=2
        labels[m]=2
        for k in range(s+1,m):
            labels[k]=1
    return labels

def evalTAD(hic,TAD1,TAD2):
    start1, end1 = readTAD(TAD1)
    start2, end2 = readTAD(TAD2)
    print(len(start1))
    print(len(start2))
    label1=getLabel(hic,start1,end1)
    label2=getLabel(hic,start2,end2)
    AMI=metrics.adjusted_mutual_info_score(label1, label2)
    RI=metrics.rand_score(label1, label2)
    AR=metrics.adjusted_rand_score(label1, label2)
    HS=metrics.homogeneity_score(label1, label2)
    VMS=metrics.v_measure_score(label1, label2)
    FMS=metrics.fowlkes_mallows_score(label1, label2)
    return AMI, RI, AR, HS, VMS, FMS

def calcuDist(arr, e):
    size = len(arr)
    idx = 0
    val = abs(e - arr[idx][0])
    
    for i in range(1, size):
        val1 = abs(e - arr[i][0])
        if val1 < val:
            idx = i
            val = val1
    if arr[idx][0] < e and arr[idx][1] < e:
        return e - arr[idx][1]
    elif arr[idx][0] < e and arr[idx][1] > e:
        return 0
    else:
        return e - arr[idx][0]

def getlist(tadfile,ctcf):
    #tad = "/home/ghaiyan/project/CASPIAN/evaluate_TADS/GM12878/chr19_5kb/TAD/{}.txt".format(name)  
    #tad = "/home/ghaiyan/project/CASPIAN/evaluate_TADS/GM12878/chr19_5kb/TAD/compare/{}.txt".format(name)
    distances = []
    with open(tadfile) as tad:
        for num, line in enumerate(tad):
            line = line.split()
            start = int(line[1])
            end = int(line[2])
            dist_start = calcuDist(ctcf, start)
            dist_end = calcuDist(ctcf, end)
            if abs(dist_start) <= abs(dist_end):
                distances.append(dist_start)
            else:
                distances.append(dist_end)
        tad.close()
    return list(set(distances))

def getctcf(factorname,chr):
    filename = "/home/ghaiyan/project/CASPIAN/chip-seq/GM12878/{}.bed".format(factorname)
    ctcf=[]
    with open(filename, 'r') as file_to_read:
        for i, line in enumerate(file_to_read):
            line = line.strip().split()
            chrname="chr"+str(chr)
            if line[0] == chrname:
                ctcf.append([int(line[1]), int(line[2])])
        file_to_read.close()
    return ctcf

def getCount(tadlist):
    count=0
    i=0
    for i in range(len(tadlist)):

        if(abs(tadlist[i])<20000):
            count=count+1
            i=i+1
        else:
            i=i+1
    print(count)
    countratio=count/len(tadlist)
    return count,countratio



print("start to read Hi-C matrix")
hic=np.loadtxt(hicfile)
for metricName in metricList:
    starttime = datetime.datetime.now()
    clusterer = hdbscan.HDBSCAN(metric=metricName)
    clusterer.fit(hic)
    label=clusterer.labels_
    print("labels:",clusterer.labels_)
    print("maxlable:",clusterer.labels_.max())
    
    np.savetxt(dir+"/{}_{}.label".format(hicname,metricName),label)
    boundaries=boundaryPlot(label)
    print("boundaries:",boundaries)
    
    np.savetxt(dir+"/{}_{}.label".format(hicname,metricName),boundaries)
    tadfile=dir+"/{}_{}.tad".format(hicname,metricName)
    i = 0
    with open(tadfile, "w") as out:
        while i < len(boundaries):
            if boundaries[i] < i:
                start = i - 1
                while i<len(boundaries) and boundaries[i] == start:
                    end = i
                    i = i + 1
                if end-start>=(2):
                    startbin = start * res
                    endbin = end * res
                    out.write("\t".join((str(start), str(startbin), str(end), str(endbin))) + "\n")
                else:
                    start=start-1
            i = i + 1
        out.close()
    endtime = datetime.datetime.now()
    print (endtime - starttime)
    tadsname=str(hicname)+str(metricName)
    n = len(hic)
    print(n)
    start, end = readTAD(tadfile)
    print("length of TADs",len(start))
    lentad=len(start)
    tad_label=start+end
    tad_label.sort()
    #print(tad_label)

    palette=sns.color_palette("bright",10)
    #print(labels)
    plt.figure(figsize=(10.5,10))
    start1=0
    end1=300
    sns.heatmap(data=hic[start1:end1, start1:end1], robust=True,cmap="OrRd")

    for i in range(0,lentad):
        if start1<start[i]<end1 and start1<end[i]<end1:
            #print(start[i])
            plt.hlines(y=start[i]-start1,xmin=start[i]-start1,xmax=end[i]-start1)
            plt.vlines(x=end[i]-start1,ymin=start[i]-start1,ymax=end[i]-start1)
    plt.title('TAD boundary')
    plt.savefig(dir+'/{}_{}_{}_{}.pdf'.format(hicname,metricName,start1,end1), format='pdf', bbox_inches='tight')
    plt.show()
    
    

    qualityFile=dir+"/quality-metrics.txt"
    with open(qualityFile,'a+') as f:
            f.write("\t".join((str(hicname),str(metricName),str(lentad),str((endtime - starttime))))+'\n')
            f.close()





