def abe_score(local_min):
    return local_min['deltaG'].astype('float').mean()/len(local_min.iat[0,0])