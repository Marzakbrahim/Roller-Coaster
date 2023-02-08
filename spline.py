# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 10:27:47 2022

@author: adrie
"""
"""______________________________________IMPORTS_________________________________________"""
import numpy as np
import matplotlib.pyplot as plt 
from numpy.linalg import *
import math as m
import matplotlib.animation as animation

import mpl_toolkits.mplot3d.axes3d as p3
import doctest

"""____________________________________FONCTIONS_________________________________________"""
class Spline :
    '''
    Classe modélisant une spline, permettant de creer notre circuit
    ---------
    Attributs : 
        pts : les points de contrôles de la spline
    --------
    
    '''
    
    def __init__(self, pts):
        """
        Fonction initialisant un spline par ses points de contÃ´les
        Parametre : 
            self : l'objet courant de type spline
            pts : Listes de points (liste de coordonnées en dim 3)
        Resultat : Cette fonction ne retourne rien

        """
        self.pointsControle = np.asarray(pts, dtype=float)
    
    def __str__(self) :
        """
        Fonction afichant un spline (ie : elle affiche ses points de controles)
        Parametre : 
        self : l'objet courant de type spline
        Resultat : une chaine de caracteres commence par "Liste des points de controles :" suivi par la listes des points de controles 
        
        """
        return 'Liste des points de contrôles : ' + str(self.pointsControle)
    
    
    def gradientPoints(self, axe):
        """
        Cette fonction nous permet de résoudre le système matriciel AX=B(ou AY=B OU AZ=B) où X contient les éléments qui nou aide à générer les points de la spline entres les points de contrôles  
        Paramètre : 
        self : l'objet courant de type spline
        axe : entier (=0 ou 1 ou 2) ; on travaille sur des points en 3D alors chaque point sera présenté par 3 composants et chaque composant sera généré par une polynome, "axe" corresponds à ces composants 
        Résultat : tableau de type array (A^-1*B)
        
        """
        n = len(self.pointsControle)
        A = np.diag(4*np.ones(n)) + np.diag(np.ones(n-1), k = -1) + np.diag(np.ones(n-1), k = 1)
        A[n-1,0] = 1
        A[0,n-1] = 1
        
        B = np.zeros((n,1))
        B[0] = 3*(self.pointsControle[1,axe] - self.pointsControle[n-1 ,axe])
        B[n - 1] = 3*(self.pointsControle[0,axe] - self.pointsControle[n-2,axe])
        
        for i in range(1, n-1) :
            B[i] = 3*(self.pointsControle[i+1,axe] - self.pointsControle[i-1,axe])
        
        return solve(A,B)
        
        
    def genereCoeffs(self, k, axe) :
        """
        Cette fonction génère les coefficients de polynôme qui represente les points de la spline entre le point de contrôle k et le point k+1  
        Paramètre : 
        self : l'objet courant de type spline
        k : intier ; c'est l'indice du 1er point de contrôle qui précise l'intervale qu'on veut définir par une polynôme
        axe : entier (=0 ou 1 ou 2) ; on travaille sur des points en 3D alors chaque point sera présenté par 3 composants et chaque composant sera généré par une polynome, "axe" corresponds à ces composants   
        Résultat : Liste de 4 floats correspent aux coefficients du polynôme qu'on veut construire
        
        """
        D = self.gradientPoints(axe);
        
        a0 = self.pointsControle[k,axe]
        a1 = float(D[k])
        a2 = float(3*(self.pointsControle[k+1,axe] - self.pointsControle[k,axe]) - 2*D[k] - D[k+1])
        a3 = float(2*(self.pointsControle[k,axe] - self.pointsControle[k+1,axe]) + D[k] + D[k+1])
    
        return [a0, a1, a2, a3]       
    
    def generePolynomes(self):
        """
        Cette fonction génère un dictionnaire dont chaque clé est l'indice du 1er point de l'intervalle [k,k+1] est sa valeur sera une liste de coefficients du polynôme qui évalue les points de cet intervalle.  
        Paramètre : 
        self : l'objet courant de type spline
        Résultat : un dictionnaire dont chaque clé est l'indice du 1er point de l'intervalle [k,k+1] est sa valeur sera une liste de coefficients du polynôme qui évalue les points de cet intervalle.
        """
        pt0 = self.pointsControle[0]
        spline_copy = Spline(np.vstack((self.pointsControle, pt0)))
        n = len(self.pointsControle)
        dico_polyn = dict()
        
        for k in range(n):
            
            dico_polyn[k] = (spline_copy.genereCoeffs(k,0), spline_copy.genereCoeffs(k, 1), spline_copy.genereCoeffs(k, 2))
        return dico_polyn
    
    def estimerLongueur(self, pas):
        """
        Cette fonction estime la longueur d'un spline
            self : l'objet courant de type spline
            pas: un flottant qui représente la distance souhaitée entre les points fournits par ParcoursUnit
        Résultat :
            spline: c'est un tableau np qui représente notre objet spline sous forme d'un tableau contenant les points de contrôles d'une manière verticale
            longueur : c'est un float qui représente la longueur de la spline
        """
        n = len(self.pointsControle)
        longueur = 0
        lst_coeffs = self.generePolynomes()
        spline = np.array([[],[],[]]).T
        
        for k in range(n):
            points = ParcoursUnit(lst_coeffs[k], pas)
            for i in range(len(points)-1):
                A = points[i]
                B = points[i+1];
                longueur += distance3D(A, B)
            spline = np.vstack((spline, points[:-1]))
        pt_initial = [self.pointsControle[0,0], self.pointsControle[0,1], self.pointsControle[0,2]]
        spline = np.vstack((spline, pt_initial))
        return spline, longueur
    
    
    def genererPoints(self, N) :
        """
        Cette fonction génére plusieurs points entre les points de contrôles qu'on a.
        Paramètres :
            self : l'objet courant de type spline
            N : on se sert de cet entier pour construir le pas des abscisses(pas=1/N**2)
        Résultat :
            points_finaux : c'est un tableau qui manifeste les nouveaux points qu'on a créé d'une manière verticale
        """
        pas = 1/N**4
        A = self.pointsControle[0]
        points_finaux = np.array([A])
        
        spline, longueur = self.estimerLongueur(pas)
        step = longueur/N
        n = spline.shape[0]
        for i in range(1,n):
            B = spline[i]
            d = distance3D(A, B)
            
            if d > step :
                A = B
                points_finaux = np.vstack((points_finaux, A))
                
        return points_finaux
    
    
    
def ParcoursUnit(lst_coeffs, pas):
    """
    Cette fonctionne génère plusieurs points entre deux points de contrôles

    Paramètres :
        lst_coeffs : c'est un dictionnaire dont chaque clé sera l'indice d'un axe et la valeur correspondante est les coefficients du polynôme qui représente les paramètres dans cet axe.
        pas : c'est le pas entre chaques deux calcules. 
    Returns :
        La fonctionne retourne un tableau de points entre deux points de contrôles
    """
    points = np.array([[],[],[]]).T
    t = 0
    epsilon = 10**(-6)
    while t <= 1 + epsilon:
        x = calculPolynome(lst_coeffs[0], t)
        y = calculPolynome(lst_coeffs[1], t)
        z = calculPolynome(lst_coeffs[2], t)
        
        points = np.vstack((points, [x,y,z]))
        t += pas
    return points

            
            
def calculPolynome(lst_coeffs, t) :
    """
    Cette fonctionne nous fait calculer le résultat d'un polynôme appliqué sur un réel

    Paramètres :
        lst_coeffs : c'est une liste de floats qui représente les coefficients du polynôme.
        t : le point sur lequel on veux appliquer le polynôme
    Returns :
        La fonctionne retourne un réel qui répresente le résultat de l'évaluation du polynôme sur t"polynôme(t)"
    >>> print(calculPolynome([0,1,2,3], 1))
    6
    """
    n = len(lst_coeffs)
    somme = 0
    for i in range(n) :
        somme += lst_coeffs[i]*t**i
    return somme
            
            
            
def distance3D(A,B) :
    """
    Cette fonctionne nous permet de savoir la distance entre deux points dans la D3

    Paramètres :
    A,B : deux listes dont chaqun contient les paramètres d'un point

    Returns :
    La fonctionne retourne un rélle qui répresente la distance entre A et B.
    >>> print(distance3D([0,0,0],[1,1,1]))
    1.7320508075688772
    """
    #Renvoie la distance entre deux points 3D A et B
    return m.sqrt((A[0] - B[0])**2 + (A[1] - B[1])**2 + (A[2] - B[2])**2)



"""______________________________________TESTS___________________________________________"""
if __name__ == '__main__' :

    S = Spline([(1,2,3), (3,2,1), (1,1,1), (0,0,0)])
    
    Dx = S.gradientPoints(0)
    Dy = S.gradientPoints(1)
    Dz = S.gradientPoints(2)
    
    fx = S.genereCoeffs(0,0)
    fy = S.genereCoeffs(0,1)
    fz = S.genereCoeffs(0,2)
    
    print(Dx)
    print(Dy)
    print(Dz)
    
    print(fx)
    print(fy)
    print(fz)
    
    print(S)
    
    dico = S.generePolynomes()
    pts = ParcoursUnit(dico[3], 0.01)
    
    spline, longueur = S.estimerLongueur(0.1)
    print(longueur)
    
    
    plt.close('all')
    
    fig = plt.figure()
    pts_finaux = S.genererPoints(5)
    
    
    
    x = np.hstack((pts_finaux[:,0], S.pointsControle[0,0]))
    y = np.hstack((pts_finaux[:,1], S.pointsControle[0,1]))
    z = np.hstack((pts_finaux[:,2], S.pointsControle[0,2]))
    
    
    ax = fig.add_subplot(projection='3d')
    pointAnime, = ax.plot([], [], [], marker="o", color = 'yellow')
    fig.set_facecolor('black')
    
    ax.set_facecolor('black')
    ax.grid(False)
    '''ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
    ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
    ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))'''
    
    
    
    ax.plot(x, y, z, color = 'white')
    ax.scatter(S.pointsControle[:,0], S.pointsControle[:,1], S.pointsControle[:,2], color = 'red')
    ax.scatter(x,y,z, marker = '.', color = 'blue')
    
    move = 0.5
    def maFonctionDanimation1(t):
        '''ax.set_xlim((x[t] - move, x[t] + move))
        ax.set_ylim((y[t] - move, y[t] + move))
        ax.set_zlim((z[t] - move, z[t] + move))'''
        pointAnime.set_data_3d(x[t],y[t],z[t])
        
        #ax.view_init(30, m.cos(t))
        
        return pointAnime,
    
    mid_axis_x = sum(ax.get_xbound())/2
    mid_axis_y = sum(ax.get_ybound())/2
    mid_axis_z = sum(ax.get_zbound())/2
    
    cam_x1 = (ax.get_xbound()[0], mid_axis_x)
    cam_x2 = (mid_axis_x, ax.get_xbound()[1])
    
    cam_y1 = (ax.get_ybound()[0], mid_axis_y)
    cam_y2 = (mid_axis_y, ax.get_ybound()[1])
    
    cam_z1 = (ax.get_ybound()[0], mid_axis_y)
    cam_z2 = (mid_axis_y, ax.get_ybound()[1])
    
    def maFonctionDanimation2(t):  
        
        if x[t] <= mid_axis_x :
            ax.set_xlim(cam_x1)
        else : 
            ax.set_xlim(cam_x2)
        
        if y[t] <= mid_axis_y :
            ax.set_ylim(cam_y1)
        else : 
            ax.set_ylim(cam_y2)
            
        '''if z[t] <= mid_axis_z :
            ax.set_zlim(cam_z1)
        else : 
            ax.set_zlim(cam_z2)'''
            
        pointAnime.set_data_3d(x[t],y[t],z[t])
        
        #ax.view_init(30, m.cos(t))
        
        return pointAnime,
    
    #monAnimation = animation.FuncAnimation(fig=fig, func=maFonctionDanimation1, frames=range(x.shape[0]), interval=100, blit=False)
    
    plt.show()
    
    doctest.testmod()
    
    
    
    







