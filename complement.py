def complement(n):
    complemento = []
    for i in n:
        if i == "A":
            complemento.append("U")
        elif i == "U":
            complemento.append("A")
        elif i == "G":
            complemento.append("C")
        elif i == "C":
            complemento.append("G")
            
    complemento = "".join(complemento)
    inverso = []
    i = len(complemento)
    while i > 0: 
        inverso += complemento[i-1]
        i = i - 1 
    complemento = inverso
    complemento = "".join(complemento)
    return(complemento)