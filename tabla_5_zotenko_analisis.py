from __future__ import division
import networkx as nx
import numpy as np

def ldata(archive):
    f=open(archive)
    data=[]
    for line in f:
        line=line.strip()
        col=line.split()
        data.append(col)
    return data

# archivos de redes
file_y2h = "./dataset/yeast_Y2H_curado.txt"
file_apms = "./dataset/yeast_AP-MS_curado.txt"
file_lit = "./dataset/yeast_LIT_curado.txt"
file_reg = "./dataset/yeast_LIT_Reguly_curado.txt"

# archivos de pares totales
pares_y2h = "./dataset/pares_totales_y2h.txt"
pares_apms = "./dataset/pares_totales_apms.txt"
pares_lit = "./dataset/pares_totales_lit.txt"
pares_reg = "./dataset/pares_totales_reg.txt"

# cargo los alpha y beta de cada red. La 1ra tupla tiene los valores simulados, la 2da los valores de las regresiones
# siempre el primero es alpha y el 2do beta
alpha_beta_y2h = [(0.009, 0.11), (0.0206, 0.1606)]
alpha_beta_apms = [(0.035, 0.19), (0.0406, 0.2221)]
alpha_beta_lit = [(0.04, 0.2), (0.0724, 0.2803)]
alpha_beta_reg = [(0.04, 0.14), (0.0485, 0.0575)]

# acomodo las cosas para iterar
files_redes = [file_y2h, file_apms, file_lit, file_reg]
files_pares = [pares_y2h, pares_apms, pares_lit, pares_reg]
nombres = ["red Y2H", "red AP-MS", "red LIT", "red REG"]
alphas_y_betas = [alpha_beta_y2h, alpha_beta_apms, alpha_beta_lit, alpha_beta_reg]

# lista de esenciales:
lista_esenciales_raw = ldata("./dataset/Essential_ORFs_paperHe_curado.txt")
lista_esenciales = []

for l in lista_esenciales_raw:
    lista_esenciales.append(l[0])

set_esenciales = set(lista_esenciales)

# defino una funcion que calcula la P_E de un nodo perteneciente a un grafo, con valores de alpha y beta establecidos por un indice "red", donde ademas hay que decir si son simulados o experimentales:

# POR DEFECTO usa los valores simulados
def P_E(nodo, grafo, red, simulado=True):

    k = grafo.degree(nodo)
    if simulado == True:
        alpha = alphas_y_betas[red][0][0]
        beta = alphas_y_betas[red][0][1]
    if simulado == False:
        alpha = alphas_y_betas[red][1][0]
        beta = alphas_y_betas[red][1][1]
    p = 1 - ( (1 - alpha)**k )*(1 - beta)
    
    return p

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

## en el caso experimental me estan quedando alpha negativos???

# aca defino si trabajo usando los simulados o los experimentales
es_simulado = True

if es_simulado == True:
    print("Usando parametros SIMULADOS")
    
    results = open("./numero_esperado_simulados.txt", "w")
    results.write("Usando parametros simulados:\n\n")

if es_simulado == False:
    print("Usando parametros EXPERIMENTALES")
    
    results = open("./numero_esperado_experimentales.txt", "w")
    results.write("Usando parametros experimentales:\n\n")

results.write("\t\t\tReal\tModelado\n")
    
for i in range(len(nombres)):
    # defino el grafo para pedir el grado de los nodos
    lista = ldata(files_redes[i])
    G = nx.Graph()
    G.add_edges_from(lista)
    
    # selecciono la lista de pares TOTALES correspondiente a la red:
    lista_pares = ldata(files_pares[i])
    
    mismo_tipo = []
        
    # separo CUALES pares son del mismo tipo:
    for l in lista_pares:
        set_l = set(l)
        interseccion = set_l.intersection(set_esenciales)
        
        if (len(interseccion) == 0) or (len(interseccion) == 2):
            mismo_tipo.append(l)

    # los cuento:
    cant_mismo_tipo = len(mismo_tipo)

    # teniendo los pares del mismo tipo, para cada nodo calculo P_E
    proba_mismo_tipo = []

    for m in mismo_tipo:
        pmt = P_E(m[0], G, i, es_simulado)*P_E(m[1], G, i, es_simulado) + (1 - P_E(m[0], G, i, es_simulado))*(1 - P_E(m[1], G, i, es_simulado))
        proba_mismo_tipo.append(pmt)
        
    # la suma de todas las probabilidades de que cada par sea del mismo tipo:
    numero_esperado = sum(proba_mismo_tipo)
    
    print("=======")
    print(nombres[i])
    print("=======")
    print("Numero esperado = %d" % (numero_esperado))
    print("Numero real = %d" % cant_mismo_tipo)
    
    results.write("%s\t\t%d\t\t%d\n" % (nombres[i], cant_mismo_tipo, numero_esperado))
    
results.close()
    
