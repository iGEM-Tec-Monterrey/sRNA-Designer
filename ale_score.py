"Luego creamos la función para iterar entre el archivo de pandas"
from loop_counter import *
def ALEscore(dfMinLocal):
    dSumatoria=0
    for i in range(len(dfMinLocal)):
        if float(loop_counter(dfMinLocal.iat[i,0]))>0 and float(dfMinLocal.iat[i,1]) != 0:
            dSumatoria=dSumatoria+loop_counter(dfMinLocal.iat[i, 0])/float(dfMinLocal.iat[i, 1])
        else:
            dSumatoria=dSumatoria+0

    dALEscore=dSumatoria/len(dfMinLocal)

    return dALEscore