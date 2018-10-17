
dulos ###########

from __future__ import division
import networkx as nx
import matplotlib.pylab as plt
from numpy import linspace

ftsize = 16

# Input files #################
# Funcion para leer input files
def ldata(archive):
    f=open(archive)
    data=[]
    for line in f:
        line=line.strip()
        col=line.split()
        data.append(col)
    return data

# Input file paths
archive_bin="./dataset/yeast_Y2H_curado.txt"
archive_mul="./dataset/yeast_AP-MS_curado.txt"
archive_lit="./dataset/yeast_LIT_copy.txt"
archive_reg="./dataset/yeast_LIT_Reguly_curado.txt"

# Nodos esenciales
# Abrir archivo
l_essential = ldata("./dataset/Essential_ORFs_paperHe_curado.txt")

###### Grafo bin
lista_bin = ldata(archive_bin)
Gbin = nx.Graph()
Gbin.add_edges_from(lista_bin)

###### Grafo mul
lista_mul = ldata(archive_mul)
Gmul = nx.Graph()
Gmul.add_edges_from(lista_mul)

###### Grafo lit
lista_lit = ldata(archive_lit)
Glit = nx.Graph()
Glit.add_edges_from(lista_lit)

###### Grafo reg
lista_reg = ldata(archive_reg)
Greg = nx.Graph()
Greg.add_edges_from(lista_reg)

####################################
def overlap_ess(lista):
	fraction = 0
	longitud = len(lista)
	num = 0
	for x in lista:
		for y in l_essential:
			if (x[0] == y[0]):
				num = num + 1
				fraction = num/longitud 
	return fraction

######################################
def fraccion_hubs_esenciales(grafo, cant_puntos):
    D_grados = grafo.degree()
    items = D_grados.items()
    sorted_items = sorted(items, key=lambda tup: tup[1])

    fracciones = linspace(0, 1, cant_puntos)
    overlap_hubs_ess = []
    for f in fracciones:
        print(f)
        cant_hubs = int(f*(len(sorted_items)))
        hubs = sorted_items[-cant_hubs:]
        overlap_hubs_ess.append(overlap_ess(hubs))
    return fracciones, overlap_hubs_ess
    
resultados_bin = fraccion_hubs_esenciales(Gbin, 1000)
resultados_mul = fraccion_hubs_esenciales(Gmul, 1000)
resultados_lit = fraccion_hubs_esenciales(Glit, 1000)
resultados_reg = fraccion_hubs_esenciales(Greg, 1000)
plt.plot(resultados_bin[0], resultados_bin[1], label='BIN')
plt.plot(resultados_mul[0], resultados_mul[1], label='MUL')
plt.plot(resultados_lit[0], resultados_lit[1], label='LIT')
plt.plot(resultados_reg[0], resultados_reg[1], label='REG')

plt.grid()
plt.legend()
plt.ylim([0,1])
plt.xlabel('fraccion de nodos que son hubs', fontsize = ftsize)
plt.ylabel('fraccion de hubs que son esenciales', fontsize = ftsize)
#plt.savefig('./figura_1_zotenko.png')
plt.show()
