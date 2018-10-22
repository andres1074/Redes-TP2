from __future__ import division
import networkx as nx
import numpy as np
from operator import itemgetter
import matplotlib.pylab as plt
from scipy.stats import linregress

ftsize = 16

def ldata(archive):
    f=open(archive)
    data=[]
    for line in f:
        line=line.strip()
        col=line.split()
        data.append(col)
    return data
    
def despeje(param):
    desp = 1 - np.exp(param)
    return desp

# titulos figuras
tit_y2h = "Y2H"
tit_apms = "AP-MS"
tit_lit = "LIT"
tit_reg = "REG"

titulos = [tit_y2h, tit_apms, tit_lit, tit_reg]

# nombres figuras
path_y2h = "./fig_2_He_y2h.png"
path_apms = "./fig_2_He_apms.png"
path_lit = "./fig_2_He_lit.png"
path_reg = "./fig_2_He_reg.png"

paths = [path_y2h, path_apms, path_lit, path_reg]

# archivos de redes
file_y2h = "./dataset/yeast_Y2H_curado.txt"
file_apms = "./dataset/yeast_AP-MS_curado.txt"
file_lit = "./dataset/yeast_LIT_curado.txt"
file_reg = "./dataset/yeast_LIT_Reguly_curado.txt"

files = [file_y2h, file_apms, file_lit, file_reg]

# lista de esenciales:
lista_esenciales_raw = ldata("./dataset/Essential_ORFs_paperHe_curado.txt")
lista_esenciales = []

for l in lista_esenciales_raw:
    lista_esenciales.append(l[0])

set_esenciales = set(lista_esenciales)

# hasta que grado grafico antes de perder la linealidad

# para y2h: 14
# para apms: 36 o 20
# para lit: 12
# para reg: 19
gs = [14, 20, 12, 19]

results = open("./alphas_y_betas_regresiones.txt", "w")
results.write("\tALPHA\tBETA\tstd\tR**2\n")

for j in range(len(titulos)):
    print("Red " + titulos[j])
    print("=======")
    ############ PARAMETROS ############
    g = gs[j]

    # elijo la red
    lista = ldata(files[j])

    # titulo del grafico
    titulo = "Red " + titulos[j]

    # output path de la figura
    path = paths[j]

    ####################################

    G = nx.Graph()
    G.add_edges_from(lista)

    grados = dict(G.degree())
    por_grado = sorted(grados.items(), key=lambda kv: kv[1])
    #por_grado = sorted(grados,key=itemgetter(1),reverse=True)

    max_grados = max([i for i in grados.values()]) 

    # separo los nodos de la red por grado.

    # coordenadas de los puntos:
    coord_x = []
    coord_y = []
    
    for i in range(1, max_grados + 1):
        grado = i
        
        nodos_del_grado = []
        for n in por_grado:
            if n[1] == grado:
                nodos_del_grado.append(n[0])
        
        # habra grados que no esten en la red:
        if len(nodos_del_grado) == 0:
            continue
        
        # ahora interseco los nodos de grado 1 con los esenciales
        nodos_del_grado = set(nodos_del_grado)
        interseccion = nodos_del_grado.intersection(set_esenciales)

        # ahora calculo la proporcion del total de ese grado que representan:
        P_E = len(interseccion)/len(nodos_del_grado)

        # la coordenada y del punto es:
        # OJO! cuando P_E = 1, el log da -inf (?)
        y = np.log(1 - P_E)
        
        coord_x.append(grado)
        coord_y.append(y)

    # regresion lineal
    # recorto la cantidad de puntos que ajusto:

    coord_x = coord_x[:g]
    coord_y = coord_y[:g]

    slope, intercept, r_value, p_value, std_err = linregress(coord_x, coord_y)
    print(slope, intercept)

    t = np.linspace(0, g, 100)
    f = slope*t + intercept
    """
    plt.figure()
    plt.plot(coord_x, coord_y, 'o')
    plt.plot(t, f, '--r')
    plt.grid()
    plt.xlabel("grado", fontsize = ftsize)
    plt.ylabel("$ln(1 - P_E)$", fontsize = ftsize)
    plt.suptitle(titulo, fontsize = ftsize)
    plt.savefig(path)
    #plt.show()
    plt.close()
    """
    
    alpha = despeje(slope)
    beta = despeje(intercept)
    
    results.write("%s\t\t%.4f\t%.4f\t%.4f\t%.4f\n" % (titulos[j], alpha, beta, std_err, r_value**2))

results.close()