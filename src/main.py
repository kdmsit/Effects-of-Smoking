import numpy as np
import pandas as pd
from numpy.linalg import matrix_rank
from scipy.stats import f
import matplotlib.pyplot as plt
from tabulate import tabulate

def readGeneData(inputfilePath):
    data = pd.read_csv(inputfilePath, sep=" ").values
    genedata = []
    GeneSymbol = []
    for probe in data:
        g = []
        gene = probe[0].split('\t')
        GeneSymbol.append(gene[49])
        gene.remove(gene[0])
        for i in range(len(gene)):
            if (i < 48):
                g.append(float(gene[i]))
        genedata.append(g)
    return genedata,GeneSymbol

def generateModelMatrix():
    A = []
    A_hat = []
    for i in range(4):
        for j in range(12):
            if i == 0:
                A.append([0, 0, 1, 0])
                A_hat.append([0, 1, 1, 0])
            elif i == 1:
                A.append([1, 0, 0, 0])
                A_hat.append([1, 0, 1, 0])
            elif i == 2:
                A.append([0, 0, 0, 1])
                A_hat.append([0, 1, 0, 1])
            elif i == 3:
                A.append([0, 1, 0, 0])
                A_hat.append([1, 0, 0, 1])
    return A,A_hat

def calculatepvalue(genedata):
    P_VALUE = []
    for probe in genedata:
        a = np.matmul(np.matmul(A,np.linalg.pinv(np.matmul(np.transpose(A),A))),np.transpose(A))
        b = np.matmul(np.matmul(A_hat, np.linalg.pinv(np.matmul(np.transpose(A_hat), A_hat))), np.transpose(A_hat))
        c=(np.matmul(np.matmul(np.transpose(probe),np.subtract(a,b)),probe))/(np.matmul(np.matmul(np.transpose(probe),np.subtract(np.identity(48),a)),probe))
        dfn=48-matrix_rank(A)
        dfd=matrix_rank(A)-matrix_rank(A_hat)
        fstat=c*(dfn)/(dfd)
        pvalue=f.cdf(fstat,dfd,dfn)
        P_VALUE.append(1-pvalue)
    return P_VALUE

def plothistogram(P_VALUE):
    plt.hist(P_VALUE)
    print("histogram.png is saved in home directory.")
    plt.savefig('histogram.png')
    plt.show()

def shortlistGeneSymbols(P_VALUE,GeneSymbol):
    gene_shortlist = []
    for i in range(len(P_VALUE)):
        g = []
        if (P_VALUE[i] <= 0.05):
            g.append(GeneSymbol[i])
            g.append(P_VALUE[i])
            gene_shortlist.append(g)
    return gene_shortlist

def readCancerGeneList():
    xenobioticMetabolism=[]
    freeRadicalResponse=[]
    dnaRepair=[]
    nkCellCytotoxicity=[]
    xenobioticMetabolismfile = "../data/XenobioticMetabolism1.txt"
    nkCellCytotoxicityfile = "../data/NKCellCytotoxicity.txt"
    freeRadicalResponsefile = "../data/FreeRadicalResponse.txt"
    dnaRepairfile = "../data/DNARepair1.txt"
    f = open(xenobioticMetabolismfile, "r")
    k=0
    for line in f:
        if(k>=2):
            xenobioticMetabolism.append(line.strip())
        k=k+1
    f.close()
    f = open(nkCellCytotoxicityfile, "r")
    k = 0
    for line in f:
        if (k >= 2):
            nkCellCytotoxicity.append(line.strip())
        k = k + 1
    f.close()
    f = open(freeRadicalResponsefile, "r")
    k = 0
    for line in f:
        if (k >= 2):
            freeRadicalResponse.append(line.strip())
        k = k + 1
    f.close()
    f = open(dnaRepairfile, "r")
    k = 0
    for line in f:
        if (k >= 2):
            dnaRepair.append(line.strip())
        k = k + 1
    f.close()
    return xenobioticMetabolism,freeRadicalResponse,dnaRepair,nkCellCytotoxicity

def findintersectionlist(xenobioticMetabolism, freeRadicalResponse, dnaRepair, nkCellCytotoxicity,gene_shortlist):
    xenobioticMetabolism_intersectlist = []
    freeRadicalResponse_intersectlist = []
    dnaRepair_intersectlist = []
    nkCellCytotoxicity_intersectlist = []
    for gene in gene_shortlist:
        if gene[0] in xenobioticMetabolism:
            xenobioticMetabolism_intersectlist.append(gene)
        if gene[0] in freeRadicalResponse:
            freeRadicalResponse_intersectlist.append(gene)
        if gene[0] in dnaRepair:
            dnaRepair_intersectlist.append(gene)
        if gene[0] in nkCellCytotoxicity:
            nkCellCytotoxicity_intersectlist.append(gene)
    return xenobioticMetabolism_intersectlist,freeRadicalResponse_intersectlist,dnaRepair_intersectlist,nkCellCytotoxicity_intersectlist

def findGeneValueAvg(geneIntersectionList):
    geneEffect = []
    i = 0
    for gene in geneIntersectionList:
        geneid = gene[0]
        for j in range(len(GeneSymbol)):
            if(GeneSymbol[j]==geneid):
                geneAvg = [geneid]
                gdata = genedata[j]
                geneAvg.append(np.mean(gdata[0:11]))
                geneAvg.append(np.mean(gdata[12:23]))
                geneAvg.append(np.mean(gdata[24:35]))
                geneAvg.append(np.mean(gdata[36:47]))
                # geneAvg.append(np.median(gdata[0:11]))
                # geneAvg.append(np.median(gdata[12:23]))
                # geneAvg.append(np.median(gdata[24:35]))
                # geneAvg.append(np.median(gdata[36:47]))
                geneEffect.append(geneAvg)
        i = i + 1
    return geneEffect

def calGeneEffect(geneEffectlist):
    womenSMvsNSM_uplist = []
    womenSMvsNSM_downlist = []
    menSMvsNSM_uplist = []
    menSMvsNSM_downlist = []
    for i in range(len(geneEffectlist)):
        gene=geneEffectlist[i]
        geneid=gene[0]
        if(gene[1]>gene[2]):
            menSMvsNSM_downlist.append(geneid)
        else:
            menSMvsNSM_uplist.append(geneid)
        if (gene[3] > gene[4]):
            womenSMvsNSM_downlist.append(geneid)
        else:
            womenSMvsNSM_uplist.append(geneid)
    return womenSMvsNSM_uplist,womenSMvsNSM_downlist,menSMvsNSM_uplist,menSMvsNSM_downlist


if __name__ == "__main__":
    inputfilePath = "../data/Raw Data_GeneSpring.txt"
    genedata, GeneSymbol=readGeneData(inputfilePath)
    xenobioticMetabolism, freeRadicalResponse, dnaRepair, nkCellCytotoxicity=readCancerGeneList()
    A, A_hat=generateModelMatrix()
    P_VALUE = calculatepvalue(genedata)
    print("Part-1:P values Generated using 2-way ANOVA framework:")
    print(P_VALUE)
    print("\n")
    print("Part-2: Draw the Histogram of p-values:")
    plothistogram(P_VALUE)
    print("\n")
    print("Part-4: shortlist rows:")
    gene_shortlist =shortlistGeneSymbols(P_VALUE,GeneSymbol)
    print(gene_shortlist)
    print("\n")
    xenobioticMetabolism_intersectlist, \
    freeRadicalResponse_intersectlist, \
    dnaRepair_intersectlist, \
    nkCellCytotoxicity_intersectlist=\
        findintersectionlist(xenobioticMetabolism, freeRadicalResponse, dnaRepair, nkCellCytotoxicity,gene_shortlist)
    print("Part-5: Intersect with the following gene lists: Xenobiotic metabolism, Free Radical Response, DNA Repair, Natural Killer Cell Cytotoxicity:")
    print("Number of Xenobiotic metabolism Genes    "+str(len(xenobioticMetabolism_intersectlist)))
    print("Xenobiotic metabolism Genes are:")
    print(xenobioticMetabolism_intersectlist)
    print("\n")
    print("Number of Free Radical Response Genes    " + str(len(freeRadicalResponse_intersectlist)))
    print("Free Radical Response Genes are:")
    print(freeRadicalResponse_intersectlist)
    print("\n")
    print("Number of DNA Repair Genes   " + str(len(dnaRepair_intersectlist)))
    print("DNA Repair Genes are:")
    print(dnaRepair_intersectlist)
    print("\n")
    print("Number of Natural Killer Cell Cytotoxicity Genes  " + str(len(nkCellCytotoxicity_intersectlist)))
    print("Natural Killer Cell Cytotoxicity Genes are:")
    print(nkCellCytotoxicity_intersectlist)
    print("\n")
    print("\n")
    geneEffect_xenobioticMetabolism=findGeneValueAvg(xenobioticMetabolism_intersectlist)
    geneEffect_dnaRepair = findGeneValueAvg(dnaRepair_intersectlist)
    geneEffect_nkCellCytotoxicity = findGeneValueAvg(nkCellCytotoxicity_intersectlist)




    print("For Xenobiotic metabolism Genes::")
    print(tabulate(geneEffect_xenobioticMetabolism,
                   headers=['Geneid', 'Men_Non_Smoker', 'Men_Smoker', 'Women_Non_Smoker', 'Women_Smoker']))
    womenSMvsNSM_uplist, womenSMvsNSM_downlist, menSMvsNSM_uplist, menSMvsNSM_downlist = calGeneEffect(
        geneEffect_xenobioticMetabolism)
    print("Women Smokers vs non-Smokers Up Genes::", list(dict.fromkeys(womenSMvsNSM_uplist)))
    print("Women Smokers vs non-Smokers Down Genes::", list(dict.fromkeys(womenSMvsNSM_downlist)))
    print("Men Smokers vs non-Smokers Up Genes::", list(dict.fromkeys(menSMvsNSM_uplist)))
    print("Men Smokers vs non-Smokers Down Genes::", list(dict.fromkeys(menSMvsNSM_downlist)))
    print("\n")
    woman_Uplist=womenSMvsNSM_uplist
    woman_Downlist=womenSMvsNSM_downlist
    man_Uplist =menSMvsNSM_uplist
    man_Downlist =menSMvsNSM_downlist

    print("For DNA Repair Genes::")
    print(tabulate(geneEffect_dnaRepair,
                   headers=['Geneid', 'Men_Non_Smoker', 'Men_Smoker', 'Women_Non_Smoker', 'Women_Smoker']))
    womenSMvsNSM_uplist, womenSMvsNSM_downlist, menSMvsNSM_uplist, menSMvsNSM_downlist = calGeneEffect(
        geneEffect_dnaRepair)
    print("Women Smokers vs non-Smokers Up Genes::",list(dict.fromkeys(womenSMvsNSM_uplist)))
    print("Women Smokers vs non-Smokers Down Genes::",list(dict.fromkeys(womenSMvsNSM_downlist)))
    print("Men Smokers vs non-Smokers Up Genes::", list(dict.fromkeys(menSMvsNSM_uplist)))
    print("Men Smokers vs non-Smokers Down Genes::", list(dict.fromkeys(menSMvsNSM_downlist)))
    print("\n")
    woman_Uplist.extend(womenSMvsNSM_uplist)
    woman_Downlist.extend(womenSMvsNSM_downlist)
    man_Uplist.extend(menSMvsNSM_uplist)
    man_Downlist.extend(menSMvsNSM_downlist)

    print("For Natural Killer Cell Cytotoxicity Genes::")
    print(tabulate(geneEffect_nkCellCytotoxicity,
                   headers=['Geneid', 'Men_Non_Smoker', 'Men_Smoker', 'Women_Non_Smoker', 'Women_Smoker']))
    womenSMvsNSM_uplist, womenSMvsNSM_downlist, menSMvsNSM_uplist, menSMvsNSM_downlist = calGeneEffect(
        geneEffect_nkCellCytotoxicity)
    print("Women Smokers vs non-Smokers Up Genes::", list(dict.fromkeys(womenSMvsNSM_uplist)))
    print("Women Smokers vs non-Smokers Down Genes::", list(dict.fromkeys(womenSMvsNSM_downlist)))
    print("Men Smokers vs non-Smokers Up Genes::", list(dict.fromkeys(menSMvsNSM_uplist)))
    print("Men Smokers vs non-Smokers Down Genes::", list(dict.fromkeys(menSMvsNSM_downlist)))
    woman_Uplist.extend(womenSMvsNSM_uplist)
    woman_Downlist.extend(womenSMvsNSM_downlist)
    man_Uplist.extend(menSMvsNSM_uplist)
    man_Downlist.extend(menSMvsNSM_downlist)
    print("\n")
    print("\n")
    print("For OverAll All types of Genes::")
    print("Women Smokers vs non-Smokers Up Genes::", list(dict.fromkeys(woman_Uplist)))
    print("Women Smokers vs non-Smokers Down Genes::", list(dict.fromkeys(woman_Downlist)))
    print("Men Smokers vs non-Smokers Up Genes::", list(dict.fromkeys(man_Uplist)))
    print("Men Smokers vs non-Smokers Down Genes::", list(dict.fromkeys(man_Downlist)))


