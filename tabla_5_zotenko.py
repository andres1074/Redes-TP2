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

file_y2h = "./dataset/yeast_Y2H_curado.txt"
file_apms = "./dataset/yeast_AP-MS_curado.txt"
file_lit = "./dataset/yeast_LIT_copy.txt"
file_reg = "./dataset/yeast_LIT_Reguly_curado.txt"

####################################
# elijo con que red trabajo:
lista = ldata(file_y2h)

# elijo el nombre del achivo output:
f = "./pares_totales_y2h_1_vecino.txt"

####################################

G = nx.Graph()
G.add_edges_from(lista)

A = nx.adjacency_matrix(G)
A_2 = A*A

size = G.number_of_nodes()

print(size)

num_pares = 0
id_positivos = []

for i in range(size):
    for j in range(size):
        # quiero mirar solo media matriz, sin incluir la diagonal
        if j>i:
            print("elemento (%d, %d)" % (i, j))
            # quiero solo pares que no sean adyacentes:
            if A[i, j] == 0:
                # quiero que ambos nodos tengan grado mayor o igual a 3:
                # (1 en el caso de la Y2H)
                if (A_2[i, i] >= 1) and (A_2[j, j] >= 1):
                    # me fijo que haya al menos 3 vecinos en comun:
                    
                    # identifico el par que paso todas las pruebas hasta ahora:
                    id1 = (list(G.nodes()))[i]
                    id2 = (list(G.nodes()))[j]
                    
                    # guardo los vecinos de cada uno como sets:
                    vecinos1 = set(list(G.neighbors(id1)))
                    vecinos2 = set(list(G.neighbors(id2)))
                    
                    # ahora si me fijo que haya al menos 3 vecinos en comun:
                    cant_comunes = len(vecinos1.intersection(vecinos2))
                    
                    if cant_comunes >= 1:
                        # caso exitoso
                        num_pares = num_pares + 1
                        print("POSITIVO")
                        id_positivos.append(((list(G.nodes()))[i], (list(G.nodes()))[j]))
                    

print(num_pares)

results = open(f, "w")
for i in id_positivos:
    results.write("%s\t%s\n" % (i[0], i[1]))
results.close()

# obs: la cantidad de pares totales para Y2H me dio 522. Con un vecino en comun dio 22k
# para apms dio 11612
# para lit, 707
# para reg, 10793

###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################

# A[0,:] devuelve la FILA 0. Ojo, como es una fila, solo tiene una fila, y el indice de fila siempre sera 0, aunque yo le haya pedido la fila N

# cada elemento de la matriz es un PAR de nodos. Agarro un par que no sean adyacentes. Me fijo que ambos tengan grado mayor o igual a 3. Tengo que calcular:

# 0) alpha y beta. De donde salen? --> me los da Marco
# 2) No. total de pares: es N**2? --> NO. Son los pares no adyacentes con 3 (o 1, para Y2H) vecinos en comun
# 4) No. ESPERADO de pares del mismo tipo usando alpha y beta de la simulacion. Cómo se llega de alpha y beta al numero esperado? --> tendria mucho sentido que este "no. esperado" sea P_E
# 5) No. ESPERADO de pares del mismo tipo usando alpha y beta del ajuste lineal

# 1) probabilidad de que sean "del mismo tipo": Pe1*Pe2 + (1 - Pe1)*(1 - Pe2) <-- esto a que parte de la tabla corresponde? Las Pe se calculan conociendo alpha y beta
# 3) No. de pares del mismo tipo. La esencialidad se saca de la lista de esenciales, o con otro criterio? Tiene algo que ver aca la P_E?

# segun zotenko:
# P_E = 1 - (1 - alpha)**k*(1 - beta)

# Obs: Hay que excluir la AP-MS de este estudio
# Obs: a la red Y2H le piden que los pares tengan uno o mas vecinos en comun, porque aparentemente es "muy" rala