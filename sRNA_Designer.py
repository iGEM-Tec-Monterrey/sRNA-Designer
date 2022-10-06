#La presente función calcula el Average base energy score de una secuencia con base en el número de estructuras locales mínimas"
#PAQUETES PRINCIPALES: barriers y RNAsubopt"
#NO TOCAR POR FAVOR"

"""
Para referencias futuras visite:"
https://www.tbi.univie.ac.at/RNA/Barriers/barriers.1.html
https://anaconda.org/bioconda/barriers
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8244778/
https://www.sciencedirect.com/science/article/pii/S1476927115302371#bib0045
https://www.researchgate.net/publication/260950912_Basin_Hopping_Graph_A_computational_framework_to_characterize_RNA_folding_landscapes
https://academic.oup.com/bioinformatics/article/30/14/2009/2390590
https://sci-hub.yncjkj.com/10.1016/j.sbi.2006.05.010  << PARA CONSULTA GENERAL
https://pubmed.ncbi.nlm.nih.gov/10070264/
https://pubmed.ncbi.nlm.nih.gov/20192764/
https://www.tbi.univie.ac.at/RNA/RNAsubopt.1.html
https://academic.oup.com/nar/article/36/suppl_2/W70/2505776?login=true
"""

import RNA
import subprocess
from subprocess import Popen, PIPE, run
import sys
import contextlib
from contextlib import redirect_stdout
import os
import pandas as pd
import time
import matplotlib.pyplot as plt
import re
from nupack import *

"IMPORTADOR DE FUNCIONES"
from ale_score import *
from complement import *
from local_min import *
from loop_counter import *
from seq_score import *
from u_score import *
from abe_score import *
from ale_score import *



"Secuencia del HFQ Binding Domain"
sHFQ='UUUCUGUUGGGCCAUUGCAUUGCCACUGAUUUUCCAACAUAUAAAAAGACAAGCCCGAACAGUCGUCCGGGCUUUUUUUCUCGAG'
sHFQ_DNA='tttctgttgggccattgcattgccactgattttccaacatataaaaagacaagcccgaacagtcgtccgggctttttttctcgag'

"Rango energético de RNAsubopt"
dRangoE=5


def FinalScoreM1(ArchivoGen, iVentana, sHFQ, dRangoE):
    Gen = open(ArchivoGen)
    Gen = Gen.read()
    Gen=  Gen.upper()
    Gen = Gen.replace(" ","")
    Gen = Gen.replace("\n","")
    Gen = Gen.replace("T","U")

    "Necesitamos crear la ventana deslizante e ir depositando en un Pandas la info de los sRNAs"
    dfsRNA=pd.DataFrame(columns=['mRNA', 'sRNA', 'SEQscore', 'ABEscore', 'ALEscore','Uscore', 'FinalScore'])
    for i in range(len(Gen)-iVentana+1):
        mRNA=Gen[i:i+iVentana]
        sRNA=complement(mRNA)
        sRNA_HFQ=str(sRNA)+str(sHFQ)
        dUscore=u_score(sRNA_HFQ)
        dSEQscore=seq_score(sRNA_HFQ)
        dfMinLocal=local_min(sRNA_HFQ, dRangoE)
        dABEscore=abe_score(dfMinLocal)
        dALEscore=ALEscore(dfMinLocal)
        dFinalScore=dSEQscore+dUscore-dABEscore-dALEscore
        dfsRNA=dfsRNA.append({'mRNA':mRNA, 'sRNA':sRNA_HFQ, 'SEQscore':dSEQscore, 'ABEscore':dABEscore, 'ALEscore':dALEscore, 'Uscore':dUscore, 'FinalScore':dFinalScore}, ignore_index=True)

    return(dfsRNA)

"dfsRNA=FinalScoreM1('ampR.fa', 24, sHFQ,5)"
"print(dfsRNA)"
"sRNA_Final=dfsRNA.to_csv('sRNA-final5.csv')"


"dfSRNA=(FinalScore('Gen2.txt',30))"
"dfSRNA.to_csv('Final.csv')"

"MODELO II"
"El modelo propuesto por (Vazquez-Anderson et al., 2017) no solo comprende a la secuencia del sRNA, sino que también considera"
"la accesibilidad del mRNA y las condiciones intracelulares. Es importante señalar que el modelo presentado por los colaboradores"
"fué sometido a un proceso de mejoramiento basado en evidencia experimental"

"TÉRMINO 1"
"Cálculo de la energía libre de hibridación mRNA-sRNA -dGAsT"
def sRNA_mRNA_Energía(sRNA, mRNA):
    fSeq=open('SeqFxHibridacion.fa', 'w+')
    fSeq.write(sRNA)
    fSeq.write('&')
    fSeq.write(mRNA)
    fSeq.close()
    os.system('RNAcofold -a -d2 --noLP --output-format=D < SeqFxHibridacion.fa > cofold-hibridacion.csv')
    dfEnergía=pd.read_csv('cofold-hibridacion.csv')
    os.remove('cofold-hibridacion.csv')
    os.remove('SeqFxHibridacion.fa')
    dDeltaGibbs=dfEnergía.iat[0, 6]
    return float(dDeltaGibbs)

"TÉRMINOS 2 Y 3"
"Cálculo de la energía de la región target y del sRNA -dGTF, -dGAsF"
def deltaGIndividual(seq):
    fSeqDG=open('sDeltaGIND.fa', 'w')
    fSeqDG.write('>Secuencia-Temporal')
    fSeqDG.write('\n')
    fSeqDG.write(seq)
    fSeqDG.close()
    os.system("RNAfold -p -d2 --noLP < sDeltaGIND.fa > sDeltaGIND-fold.out")
    fDeltaGIND=open('sDeltaGIND-fold.out', 'r')
    sRenglonDG=fDeltaGIND.readlines()[3]
    aListDG=re.findall(r"[-+]?(?:\d*\.\d+|\d+)", sRenglonDG)
    dDeltaG=float(aListDG[0])
    os.remove('sDeltaGIND.fa')
    os.remove('sDeltaGIND-fold.out')
    return dDeltaG

"TÉRMINO 4"
"Cálculo del factor de accesibilidad"
"TÉRMINO 4"
"Cálculo del factor de accesibilidad"
def FactorAccesibilidad(mRNA):
    a = Strand(mRNA, name='a')
    my_model = Model(material='RNA', celsius=37)
    pairing = Complex([a])
    set1 = ComplexSet(strands=[a], complexes=SetSpec(max_size=1, include=[pairing]))
    complex_results1 = complex_analysis(set1, model=my_model, compute=['pairs'])
    pairing_result = complex_results1[pairing]
    "plt.imshow(pairing_result.pairs.to_array())"
    "plt.xlabel('Índice')"
    "plt.ylabel('Índice')"
    "plt.title('Gráfica de probabilidad de apareamiento individual')"
    "plt.colorbar()"
    "plt.clim(0, 1)"

    aMatrizProbabilidad=pairing_result.pairs.to_array()
    aList=[]
    for i in range(len(mRNA)):
        for j in range(i+1,len(mRNA)):
            aList.append(aMatrizProbabilidad[i][j])

    aFactorDisp = [1-x for x in aList]

    dSumFactDisp=sum(aFactorDisp)
    return dSumFactDisp


def FinalScoreMod1(ArchivoGen, iVentana, sHFQ):
    Gen = open(ArchivoGen)
    Gen = Gen.read()
    Gen=  Gen.upper()
    Gen = Gen.replace(" ","")
    Gen = Gen.replace("\n","")
    Gen = Gen.replace("T","U")

    "Necesitamos crear la ventana deslizante e ir depositando en un Pandas la info de los sRNAs"
    dfsRNA=pd.DataFrame(columns=['mRNA', 'sRNA', 'SEQscore', 'ABEscore', 'ALEscore','Uscore', 'Score Modelo 1'])
    for i in range(len(Gen)-iVentana+1):
        mRNA=Gen[i:i+iVentana]
        sRNA=complement(mRNA)
        dUscore=u_score(sRNA+sHFQ)
        dSEQscore=seq_score(sRNA+sHFQ)
        dfMinLocal=local_min(sRNA+sHFQ)
        dABEscore=abe_score(dfMinLocal)
        dALEscore=ALEscore(dfMinLocal)
        dFinalScore=dSEQscore+dUscore-dABEscore-dALEscore
        dfsRNA=dfsRNA.append({'mRNA':mRNA, 'sRNA':sRNA, 'SEQscore':dSEQscore, 'ABEscore':dABEscore, 'ALEscore':dALEscore, 'Uscore':dUscore, 'FinalScore':dFinalScore}, ignore_index=True)
    
    output=dfsRNA.to_csv('Score6.csv', index=True)
    return(dfsRNA)

"FinalScoreMod1('cloR.fa',24,sHFQ)"

def FinalScoreMod2(ArchivoGen, iVentana):
    Gen = open(ArchivoGen)
    Gen = Gen.read()
    Gen=  Gen.upper()
    Gen = Gen.replace(" ","")
    Gen = Gen.replace("\n","")
    Gen = Gen.replace("T","U")

    "Necesitamos crear la ventana deslizante e ir depositando en un Pandas la info de los sRNAs"
    dfsRNA=pd.DataFrame(columns=['ELibre Hibrid', 'ELibre RegTarget', 'ELibre Fold sRNA', 'Factor de accesibilidad', 'Delta G overall', 'Score modelo 2'])
    for i in range(len(Gen)-iVentana+1):
        mRNA=Gen[i:i+iVentana]
        sRNA=complement(mRNA)
        dELibreHibrid=sRNA_mRNA_Energía(sRNA, mRNA)
        dELibreTarget=deltaGIndividual(mRNA)
        dELibresRNA=deltaGIndividual(sRNA)
        dFactAccesib=FactorAccesibilidad(mRNA)
        dDGOverall=dELibreHibrid-dELibreTarget-dELibresRNA
        dScoreMod2=dFactAccesib*dELibreTarget+dELibreTarget*dELibreHibrid+dELibreHibrid+dELibreTarget+dFactAccesib
        dfsRNA=dfsRNA.append({'ELibre Hibrid':dELibreHibrid, 'ELibre RegTarget':dELibreTarget, 'ELibre Fold sRNA':dELibresRNA, 'Factor de accesibilidad':dFactAccesib, 'Delta G overall':dDGOverall, 'Score modelo 2':dScoreMod2}, ignore_index=True)
    output=dfsRNA.to_csv('Score5.csv', index=True)
    return output

"FinalScoreMod2('cloR.fa',22)"

"Función normalizadora"
def Normalizar(dEntrada, dMin, dMax):
    dNorm=(dEntrada-dMin)/(dMax-dMin)
    return dNorm

"FUNCIÓN INTEGRADORA FINAL MODELOS 1 Y 2"
def ScoreIntegrador(ArchivoGen, iVentana, sHFQ, dRangoE):
    Gen = open(ArchivoGen)
    Gen = Gen.read()
    Gen=  Gen.upper()
    Gen = Gen.replace(" ","")
    Gen = Gen.replace("\n","")
    Gen = Gen.replace("T","U")

    "Necesitamos crear la ventana deslizante e ir depositando en un Pandas la info de los sRNAs"
    dfsRNA=pd.DataFrame(columns=['mRNA', 'sRNA', 'SEQscore', 'ABEscore', 'ALEscore','Uscore', 'FinalScoreM1','ELibre Hibrid', 'ELibre RegTarget', 'ELibre Fold sRNA', 'Factor de accesibilidad', 'Delta G overall', 'FinalScoreM2'])
    for i in range(len(Gen)-iVentana+1):
        mRNA=Gen[i:i+iVentana]
        sRNA=complement(mRNA)
        if i==0:
            mRNAtF=Gen[i:(i+iVentana+1)]
        
        if i != 0 and i != (len(Gen)-iVentana):
            mRNAtF=Gen[(i-1):(i+iVentana+1)]
        
        if i==(len(Gen)-iVentana):
            mRNAtF=Gen[(i-1):(i+iVentana)]

        sRNA_HFQ=str(sRNA)+str(sHFQ)
        dUscore=u_score(sRNA_HFQ)
        dSEQscore=seq_score(sRNA_HFQ)
        dfMinLocal=local_min(sRNA_HFQ, dRangoE)
        dABEscore=abe_score(dfMinLocal)
        dALEscore=ALEscore(dfMinLocal)
        dFinalScoreM1=dSEQscore+dUscore-dABEscore-dALEscore
        dELibreHibrid=sRNA_mRNA_Energía(sRNA, mRNA)
        dELibreTarget=deltaGIndividual(mRNAtF)
        dELibresRNA=deltaGIndividual(sRNA_HFQ)
        dFactAccesib=FactorAccesibilidad(mRNA)
        dDGOverall=dELibreHibrid-dELibreTarget-dELibresRNA
        dScoreMod2=dFactAccesib*dELibreTarget+dELibreTarget*dELibreHibrid+dELibreHibrid+dELibreTarget+dFactAccesib
        
        dfsRNA=dfsRNA.append({'mRNA':mRNA, 'sRNA':sRNA_HFQ, 'SEQscore':dSEQscore, 'ABEscore':dABEscore, 'ALEscore':dALEscore, 'Uscore':dUscore, 'FinalScoreM1':dFinalScoreM1, 'ELibre Hibrid':dELibreHibrid, 'ELibre RegTarget':dELibreTarget, 'ELibre Fold sRNA':dELibresRNA, 'Factor de accesibilidad':dFactAccesib, 'Delta G overall':dDGOverall, 'FinalScoreM2':dScoreMod2}, ignore_index=True)
    dMaxM1=dfsRNA['FinalScoreM1'].max()
    dMinM1=dfsRNA['FinalScoreM1'].min()
    dMaxM2=dfsRNA['FinalScoreM2'].max()
    dMinM2=dfsRNA['FinalScoreM2'].min()
    dfsRNA['Score1Norm']=dfsRNA['FinalScoreM1'].apply(Normalizar, args=(dMinM1, dMaxM1))
    dfsRNA['Score2Norm']=dfsRNA['FinalScoreM2'].apply(Normalizar, args=(dMinM2, dMaxM2))
    dfsRNA['FINAL']=dfsRNA['Score1Norm']+dfsRNA['Score2Norm']

    print(dfsRNA)
    output=dfsRNA.to_csv('ScoreFN.csv', index=True)
    return output

print("sRNA Designer v1.0")
print("This python script analyzes a defined mRNA sequence, generates a dataset of candidates sRNA sequences and then applies a score series")
print("Make sure you have all files in the same directory, or indicate in which folder is your mRNA sequence: ")
fGen=str(input("Type the file in which your gene is located in FASTA format: "))
iOligoLen=int(input("Type the length of the desired sRNA sequences: "))
iAboveMFE=int(input("Type the range of free enrgy (in kcal/mol) above the MFE for local minima calculation: "))

ScoreIntegrador(fGen, iOligoLen, sHFQ, iAboveMFE)
