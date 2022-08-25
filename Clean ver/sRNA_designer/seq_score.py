"La siguiente función calcula la Sequence Score de un string dado"
def seq_score(sRNA):
    "Primero necesitamos crear un diccionario con la información acerca de la matriz de puntos"
    aDiccionarioPuntosNonGenic={
        "AA":-0.170,
        "UA":-0.302,
        "GA":-0.146,
        "CA":-0.205,
        "AU":-0.191,
        "UU":-0.075,
        "GU":-0.119,
        "CU":-0.106,
        "AG":0.293,
        "UG":0.214,
        "GG":0.193,
        "CG":-0.091,
        "AC":0.187,
        "UC":0.208,
        "GC":0.051,
        "CC":0.211
    }

    aDiccionarioPuntos={
        "AA":-0.021,
        "UA":-0.071,
        "GA":0.109,
        "CA":-0.002,
        "AU":-0.003,
        "UU":0.2,
        "GU":-0.250,
        "CU":0.033,
        "AG":-0.129,
        "UG":0.081,
        "GG":0.143,
        "CG":-0.084,
        "AC":-0.149,
        "UC":0.067,
        "GC":-0.154,
        "CC":0.247
    }

    "'Tonces 'ora contabilizamos la sumatoria de cada uno de los dinucleótidos ubicados en el sRNA"
    dSumDinucleótidos=0
    aMatrizTest=[]
    for i in range(len(sRNA)-1):
        dSumDinucleótidos=dSumDinucleótidos+aDiccionarioPuntos[sRNA[i:i+2]]
        aMatrizTest.append(aDiccionarioPuntos[sRNA[i:i+2]])
        
    dSeqScore=dSumDinucleótidos/len(sRNA)
    return dSeqScore