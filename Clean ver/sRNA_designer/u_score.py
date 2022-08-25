def u_score(sRNA):
    Ucounter = 0
    Utcounter = 0
    sRNAlenght = len(sRNA)
    sRNAlenght = sRNAlenght - 10
    sRNA10 = sRNA[sRNAlenght:len(sRNA)]
    
    for i in sRNA10:
        if i == "U":
            Ucounter = Ucounter + 1
    
    for i in sRNA:
        if i == "U":
            Utcounter = Utcounter + 1

    Uscore = ((Ucounter/10)/(Utcounter/len(sRNA)))       
    return(Uscore)