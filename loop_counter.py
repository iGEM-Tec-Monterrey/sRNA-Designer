"Analizador de loops- ALEscore"
"Esta fórmula analiza el número de loops dentro de las estructuras con minimos locales"
"Primero necesitamos generar una funcion que cuente loops, donde sEPP es la est punto paréntesis"
def loop_counter(sEPP):
    import re
    string = sEPP
    count = 0
    while len(re.findall('\(\.+\)', string))>0:
        inner_loops = re.findall('\(\.+\)', string)
        count = count + len(inner_loops)
        string = re.sub('\(\.+\)','', string)
        while len(re.findall('\(\)', string))>0:
            string = re.sub('\(\)','', string)
    
    return count