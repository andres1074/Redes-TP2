from __future__ import division
import networkx as nx
import matplotlib.pylab as plt
import numpy as np
from collections import Counter
from random import choice


# Funcion para cargar datos
def ldata(archive):
    f=open(archive)
    data=[]
    for line in f:
        line=line.strip()
        col=line.split()
        data.append(col)
    return data

# Files and locations #################
file_bin = "./dataset/yeast_Y2H_curado.txt"
file_mul = "./dataset/yeast_AP-MS_curado.txt"
file_lit = "./dataset/yeast_LIT_copy.txt"
file_reg = "./dataset/yeast_LIT_Reguly_curado.txt"

# Lista proteinas esenciales
lista_esenciales_raw = ldata("./dataset/Essential_ORFs_paperHe_curado.txt")
lista_esenciales = []

for l in lista_esenciales_raw:
    lista_esenciales.append(l[0])
########################################

# Funciones ##############################################
# 1.Funcion definir subconjunto de esenciales
def subconjunto_ess(lista):
    sublista = []
    for x in lista:
        for y in lista_esenciales:
            if (x[0] == y):
                sublista.append(x)
    return sublista

# OBS: hay algo muy poco eficiente, y es que la funcion mira la interseccion entre la red y los esenciales cada vez que se la llama. Para todas las iteraciones sobre una misma red esto es redundante

# 2.Funcion remover no esenciales 
def remover_no_esenciales(file):
    lista = ldata(file)
    G = nx.Graph()
    G.add_edges_from(lista)
    
    dict_grados = dict(G.degree())
    lista_grados = sorted(dict_grados.items(), key=lambda kv: kv[1])
    
    esenciales_red = subconjunto_ess(lista_grados)

    grado_esenciales = []
    for e in esenciales_red:
        grado_esenciales.append(e[1])    

    # counter es un diccionario que tiene como keys los grados para los que hay nodos esenciales, y como values la cantidad de nodos esenciales para ese grado
    counter = Counter(grado_esenciales)
    #print(counter)
    
    nodos_bin = G.nodes()

    grados_que_tienen_esenciales = list(counter.keys())
    cuantos_hay_que_remover = list(counter.values())
    
    for indice in range(len(counter)):
        grado_check = grados_que_tienen_esenciales[indice]
        #print("grado check = %d" % (grado_check))
        #print("================")
            
        # agrupo los nodos que tienen ese grado y que no son esenciales
        nodos_de_interes = []
        for i in dict_grados:
            if G.degree(i) == grado_check:
                nodos_de_interes.append(i)
                
        # me falta descartar los esenciales
        set_nodos_de_interes = set(nodos_de_interes)
        set_esenciales = set(lista_esenciales)
        nodos_de_interes = list(set_nodos_de_interes - set_esenciales)
        
        removidos = 0

        while removidos < cuantos_hay_que_remover[indice]:
            #print("cuantos hay que remover = %d" % (cuantos_hay_que_remover[indice]))
            #print("removidos = %d" % (removidos))
            try:
                nodo_azar = choice(nodos_de_interes)
                G.remove_node(nodo_azar)
                removidos = removidos + 1
                
                # luego de remover un nodo, actualizo mi conjunto de nodos de interes
                nodos_de_interes = []
                for i in dict_grados:
                    if G.degree(i) == grado_check:
                        nodos_de_interes.append(i)
                        
                # me falta descartar los esenciales
                set_nodos_de_interes = set(nodos_de_interes)
                set_esenciales = set(lista_esenciales)
                nodos_de_interes = list(set_nodos_de_interes - set_esenciales)
            except:
                grado_check = grado_check - 1
                
                # agrupo los nodos que tienen ese grado y que no son esenciales
                nodos_de_interes = []
                for i in dict_grados:
                    if G.degree(i) == grado_check:
                        nodos_de_interes.append(i)
                        
                # me falta descartar los esenciales
                set_nodos_de_interes = set(nodos_de_interes)
                set_esenciales = set(lista_esenciales)
                nodos_de_interes = list(set_nodos_de_interes - set_esenciales)
                
            #print("grado check = %d" % (grado_check))
            
    size_comp_gigante = len(max(nx.connected_components(G), key=len))
    fraccion = size_comp_gigante/len(dict_grados)
    
    #G = nx.Graph()
    #G.add_edges_from(lista_bin)
    
    return fraccion

# 3.Funcion removiendo los esenciales:
def remover_esenciales(file):
    lista = ldata(file)
    G = nx.Graph()
    G.add_edges_from(lista)
    
    cant_nodos = len(G.nodes())
    
    set_esenciales = set(lista_esenciales)
    set_nodos_red = set(G.nodes())
    nodos_de_interes = list(set_nodos_red.intersection(set_esenciales))
    
    for n in nodos_de_interes:
        G.remove_node(n)
        
    size_comp_gigante = len(max(nx.connected_components(G), key=len))
    fraccion = size_comp_gigante/cant_nodos
    
    return fraccion

# Fin funciones ##########################################


#lista_files = [file_bin, file_mul, file_lit, file_reg]

rem_esenciales = []
rem_no_esenciales_bin = []
rem_no_esenciales_mul = []
rem_no_esenciales_lit = []
rem_no_esenciales_reg = []

iteraciones = 1000

#### BIN ####
print("RED BIN")
print("=======")
f = file_bin
rem_esenciales.append([f, remover_esenciales(f)])

for i in range(iteraciones):
    print("Iteracion %d" % (i))
    rem_no_esenciales_bin.append(remover_no_esenciales(f))


#### MUL ####
print("RED MUL")
print("=======")
f = file_mul
rem_esenciales.append([f, remover_esenciales(f)])

for i in range(iteraciones):
    print("Iteracion %d" % (i))
    rem_no_esenciales_mul.append(remover_no_esenciales(f))


#### LIT ####
print("RED LIT")
print("=======")
f = file_lit
rem_esenciales.append([f, remover_esenciales(f)])

for i in range(iteraciones):
    print("Iteracion %d" % (i))
    rem_no_esenciales_lit.append(remover_no_esenciales(f))


#### REG ####
print("RED REG")
print("=======")
f = file_reg
rem_esenciales.append([f, remover_esenciales(f)])

for i in range(iteraciones):
    print("Iteracion %d" % (i))
    rem_no_esenciales_reg.append(remover_no_esenciales(f))


print("RESULTADOS")
print(rem_esenciales)
print(rem_no_esenciales_bin)
print(rem_no_esenciales_mul)
print(rem_no_esenciales_lit)
print(rem_no_esenciales_reg)

results_bin = (np.mean(rem_no_esenciales_bin), np.std(rem_no_esenciales_bin))
results_mul = (np.mean(rem_no_esenciales_mul), np.std(rem_no_esenciales_mul))
results_lit = (np.mean(rem_no_esenciales_lit), np.std(rem_no_esenciales_lit))
results_reg = (np.mean(rem_no_esenciales_reg), np.std(rem_no_esenciales_reg))

results = open("./tabla_3_zotenko_resultados.txt", "w")
results.write("\tRemEs\tRemNoEs\tStd\n")
results.write("BIN\t%f\t%f\t%f\n" % (rem_esenciales[0][1], results_bin[0], results_bin[1]))
results.write("MUL\t%f\t%f\t%f\n" % (rem_esenciales[1][1], results_mul[0], results_mul[1]))
results.write("LIT\t%f\t%f\t%f\n" % (rem_esenciales[2][1], results_lit[0], results_lit[1]))
results.write("REG\t%f\t%f\t%f\n" % (rem_esenciales[3][1], results_reg[0], results_reg[1]))
results.close()
