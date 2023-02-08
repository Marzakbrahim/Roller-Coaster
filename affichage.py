# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 22:43:41 2022

@author: adrie
"""
"""______________________________________IMPORTS_________________________________________"""
import numpy as np
import math as m
from spline import *
from cabine import *
import matplotlib.animation as animation

import mpl_toolkits.mplot3d.axes3d as p3



"""____________________________________FONCTIONS_________________________________________"""
def afficherCircuit(S, N, figure) :
    '''
    Fonction traçant le circuit à partir de la spline
    ---------
    Paramètres : 
        S : objet de la classe Spline
        N : nombres de points que l'on souhaite pour tracer le circuit
    --------
    Résultat : trace le graphe et renvoie 
        ax : graphique
        pointAnime : Animation du graphe
    
    '''
    lst_points = S.genererPoints(N)
    x = np.hstack((lst_points[:,0],S.pointsControle[0,0]))
    y = np.hstack((lst_points[:,1],S.pointsControle[0,1]))
    z = np.hstack((lst_points[:,2],S.pointsControle[0,2]))
    
    ax = figure.add_subplot(projection='3d')
    pointAnime, = ax.plot([], [], [], marker="o", color = 'yellow')
    
    figure.set_facecolor('black')

    ax.set_facecolor('black')
    ax.grid(False)
    ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
    ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
    ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
        
    ax.plot(x, y, z, color = 'white')
    ax.scatter(S.pointsControle[:,0], S.pointsControle[:,1], S.pointsControle[:,2], color = 'red')
    
    return ax, pointAnime
    
    

def afficherTrajectoire(C, ax, pointAnime, camera = 1, move = 0.5) :
    '''
    Fonction définissant la fonction d'animation pour le graphe et la cabine 
    en paramètres.
    ---------
    Paramètres : 
        C : objet de la classe Cabine
        ax : graphique
        pointAnime : variable d'animation
        camera : {1,2 ou 3} définit quel type de caméra on souhaite
        move : largeur de l'écran pour la caméra 2
    --------
    Résultat : trace le graphe et renvoie 
        maFonctionDanimation : fonction responsable de l'animation
    
    '''
    
    def maFonctionDanimation(t):
        
        C.update(t)
        pos_x = C.position[0]
        pos_y = C.position[1]
        pos_z = C.position[2]
        
        if camera == 2:
            
            ax.set_xlim((pos_x - move, pos_x + move))
            ax.set_ylim((pos_y - move, pos_y + move))
            ax.set_zlim((pos_z - move, pos_z + move))
            
            ax.view_init(C.direction[0], C.direction[1])
        
        elif camera == 3:
            mid_axis_x = sum(ax.get_xbound())/2
            mid_axis_y = sum(ax.get_ybound())/2
            mid_axis_z = sum(ax.get_zbound())/2

            cam_x1 = (ax.get_xbound()[0], mid_axis_x)
            cam_x2 = (mid_axis_x, ax.get_xbound()[1])

            cam_y1 = (ax.get_ybound()[0], mid_axis_y)
            cam_y2 = (mid_axis_y, ax.get_ybound()[1])

            cam_z1 = (ax.get_ybound()[0], mid_axis_y)
            cam_z2 = (mid_axis_y, ax.get_ybound()[1])

                
            if pos_x <= mid_axis_x :
                ax.set_xlim(cam_x1)
            else : 
                ax.set_xlim(cam_x2)
                
            if pos_y <= mid_axis_y :
                ax.set_ylim(cam_y1)
            else : 
                ax.set_ylim(cam_y2)
                    
            if pos_z <= mid_axis_z :
                ax.set_zlim(cam_z1)
            else : 
                ax.set_zlim(cam_z2)
            
        pointAnime.set_data_3d(pos_x,pos_y,pos_z)
        return pointAnime,
        
    
    return  maFonctionDanimation

def lire_fichier(file):
    '''
    Fonction lisant un fichier txt et renvoyant la liste des points de controles
    ---------
    Paramètres : 
        file : fichier à lire
    --------
    Résultat : trace le graphe et renvoie 
        pts_controles : liste de tuples
    
    '''
    f = open(file, 'r')
    
    if float(f.readline()) < 4 :
        print("Il n'y a pas assez de points de contrôles. Réessayez avec au moins 4 points.")
        return
    
    pts_controles = []
    point_courant = []
    for ligne in f:
        
        if ligne == '\n':
            if len(point_courant) == 3 :
                pts_controles.append(point_courant)
                point_courant = []
            continue
        
        #print('ligne : ', ligne)
        for nb in ligne.split() :
            point_courant.append(float(ligne))
            
    return pts_controles
        
    

def creer_fichier(filename):
    '''
    Fonction plaçant les variables saisies en input sous la forme d'un fichier txt
    ---------
    Paramètres : 
        filename : nom du fichier txt
    --------
    Résultat : trace le graphe et renvoie 
        None
    
    '''
    f = open(filename +'.txt', 'w')
    nb_pts = int(input('Entrez le nombre de points de contrôles : '))
    f.write(str(nb_pts) + '\n')
    f.write('\n')
    print('Entrez les points :')
    for i in range(nb_pts) :
        x = input('x' + str(i+1) + ' = ')
        y = input('y' + str(i+1) + ' = ')
        z = input('z'+ str(i+1) + ' = ')
        f.write(str(x) + '\n')
        f.write(str(y) + '\n')
        f.write(str(z) + '\n')
        f.write('\n')
    f.close()

def creer_spline_cabine(file, g = 9.8) :
    '''
    Fonction une spline et une cabine
    ---------
    Paramètres : 
        file : fichier txt contenant les points de controles
    --------
    Résultat : trace le graphe et renvoie 
        S : objet de classe Spline
        C : objet de classe Cabine
    
    '''
    S = Spline(lire_fichier(file))
    m = float(input('Entrer la masse de la cabine: '))
    k = float(input('Entrer le coefficient de frottements: '))
    v0 = float(input('Entrer la vitesse initiale de la cabine: '))
    N = float(input('Entrer le nombre de points du circuit: '))
    
    return S, Cabine(m,k,v0,S,N,g)

"""______________________________________TESTS___________________________________________"""
if __name__ == '__main__' :
    
    
    #creer_fichier('test2')
    #pts_controles = lire_fichier('test2.txt')
    
    #print(pts_controles)
    
    S,C = creer_spline_cabine('test')
    
    plt.close('all')
    figure = plt.figure()
    ax, ptAnime = afficherCircuit(S, 100, figure)
    
    
    
    maFonctionDanimation = afficherTrajectoire(C, ax, ptAnime, 0.5, camera = 1)
    
    dt = 0.01
    T = dt*np.ones(2)
    monAnimation = animation.FuncAnimation(fig=figure, func=maFonctionDanimation, frames=T, interval=50, blit=False)
    plt.show()

    
    