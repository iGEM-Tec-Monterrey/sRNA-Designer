import os
import pandas as pd
import subprocess

"Generador de estructuras locales mínimas"
"La siguiente fx genera un archivo de datos con los mínimos locales"
"Estos datos serán analizados por ALEscore y ABEscore"
def local_min(sRNA,dRangoE):
    sRangoE="-e "+str(dRangoE)
    if os.path.isfile('RNAsubopt-test.out')== True:
        os.remove('RNAsubopt-test.out')
    if os.path.isfile('barriers-test.out')== True:
        os.remove('barriers-test.out')
    txtsRNA=open('temp.txt', 'w')
    txtsRNA.write(sRNA)
    txtsRNA.close()
    "Creamos lista de esructuras subóptimas"
    pSubopt=subprocess.Popen(['RNAsubopt', sRangoE, '--infile=temp.txt', '--outfile=RNAsubopt-test.out', '-d2', '--noLP', '-s'])
    (output, err)=pSubopt.communicate()
    pSubopt_status=pSubopt.wait()
    print("RNAsubopt finalizado")

    "Llamamos a Barriers para encontrar mínimos locales"
    os.system("barriers --rates -G RNA-noLP --max 50 --minh 0.1 --quiet < RNAsubopt-test.out > barriers-test.out")
    
    "Neceistamos extraer los valores de la energía de los mínmos locales"
    fLocalMinimos=open("barriers-test.out", "r")
    aMatrizBarriers=[]
    for line in fLocalMinimos:
        stripped_line=line.strip()
        line_list=stripped_line.split()
        aMatrizBarriers.append(line_list)
    fLocalMinimos.close()

    "La matriz aMatrizBarriers tiene dentro nuestros datos a partir de la segunda fila"
    "Creamos una nueva con solo los valores de la energía en mínimos locales"
    dfMinLocal=pd.DataFrame(columns=['Estructura P/P', 'deltaG'])

    for i in range(1, len(aMatrizBarriers)):
        dfMinLocal.loc[len(dfMinLocal.index)]=[aMatrizBarriers[i][1], aMatrizBarriers[i][2]]
    
    return dfMinLocal