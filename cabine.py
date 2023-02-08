# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 21:35:09 2022

@author: adrie
"""
"""______________________________________IMPORTS_________________________________________"""
import numpy as np
import math as m
from spline import *


"""____________________________________FONCTIONS_________________________________________"""
class Cabine:
    
    '''
    Classe modélisant la cabine de notre circuit
    ---------
    Attributs : 
        indice_circuit : l'indice du circuit auquel se trouve la cabine 
                        exemple : si la cabine est entre 1 et 2, l'indice est 1.
        masse : un flottant représente la masse de la cabine
        coeff_frottement : un flottant représente le coefficient de frottement de la cabine
        vitesse : un flottant qui représente la vitesse de la cabine
        circuit : tableau np définissant le circuit sur lequel est la cabine
        direction : array définissant la direction de la cabine
        position : array définissant la position de la cabine
        G : un flottant qui représente la constance de gravité (9.8 par défaut)
    --------
        
    '''
    
    def __init__(self, m, k, v0, S, N = 100, g = 9.8):
        """
        Fonction initialisant une cabine par sa masse et son coefficient de frottememnt et le cercuit oÃ¹ elle et sa vitesse initoale 
        Parametres : 
            self : l'objet courant de type cabine
            m : un flottant représente la masse de la cabine
            k : un flottant représente le coefficient de frottement de la cabine
            v0 : un flottant qui représente la vitesse initiale de la cabine
            S : l'objet de la classe spline définissant le circuit sur lequel est la cabine
            N : le nombre de points du circuit
            g : un flottant qui représente la constance de gravité (9.8 par défaut)
        Resultat : Cette fonction ne retourne rien
        """
        points = S.genererPoints(N)
        u0 = points[1] - points[0]
        dir0 = u0/np.linalg.norm(u0)
        
        self.indice_circuit = 0
        self.masse = m
        self.coeff_frottement = k
        self.vitesse = v0
        self.circuit = points
        self.direction = dir0
        self.position = points[0]
        self.G = g
    
    def forcesFrottements(self):
        return -self.coeff_frottement*self.vitesse
    
    def get_m(self):
        """
        Fonction membre retournant la masse de la cabine
        Paramètre : 
            self : l'objet courant de type cabine    
        Résultat : 
            La fonction retourne un flottent qui représente la masse de la cabine.

        """
        return self.masse
    
    def get_k(self):
        """
        Fonction membre retournant le coefficient de frottement de l'objet.
        Paramètre : 
            self : l'objet courant de type cabine    
        Résultat : 
            La fonction retourne un flottent qui représente le coefficient de frottement de la cabine.

        """
        return self.coeff_frottement
    
    def get_v(self):
        """
        Fonction membre retournant la vitesse de la cabine.
        Paramètre : 
            self : l'objet courant de type cabine    
        Résultat : 
            La fonction retourne un flottent qui représente la vitesse de la cabine.

        """
        return self.vitesse
    
    def get_u(self):
        """
        Fonction membre retournant la direction de la cabine.
        Paramètre : 
            self : l'objet courant de type cabine    
        Résultat : 
            La fonction retourne une tuple de 3 flottants qui représente vecteur de la direction de la cabine.

        """
        return self.direction
    
    def get_pos(self):
        """
        Fonction membre retournant la position de la cabine dans le circuit (spline).
        Paramètre : 
            self : l'objet courant de type cabine    
        Résultat : 
            La fonction retourne une tuple de 3 flottants qui représente la position de la cabine dans le cercuit.

        """
        return self.position
    
    def get_circuit(self):
        """
        Fonction membre retournant la position de la cabine dans le circuit (spline).
        Paramètre : 
            self : l'objet courant de type cabine    
        Résultat : 
            La fonction retourne le tableau np des pts equidistants du circuit.

        """
        return self.circuit
    
    
    def set_v(self, new_v):
        """
        Fonction membre qui mise à jour la vitesse la cabine.
        Paramètre : 
            self : l'objet courant de type cabine 
            new_v : un flottant qui représent la nouvelle vitesse de la cabine
        Résultat : 
            La fonction ne retourne rien, elle remplace l'ancienne vitesse da la cabine par la nouvelle.

        """
        self.vitesse = new_v 
        
    def set_u(self, new_dir):
        """
        Fonction membre qui met à jour la direction la cabine.
        Paramètre : 
            self : l'objet courant de type cabine 
            new_dir : une tuple de trois flottants qui représent la nouvelle direction de la cabine
        Résultat : 
        La fonction ne retourne rien, elle remplace l'ancienne direction de la cabine par la nouvelle.

        """
        self.direction = new_dir
    
    def set_pos(self, new_pos):
        """
        Fonction membre qui met à jour la position la cabine dans le cercuit(spline).
        Paramètre : 
            self : l'objet courant de type cabine 
            new_pos : une tuple de trois flottants qui représent la nouvelle position de la cabine dans le cercuit
        Résultat : 
            La fonction ne retourne rien, elle remplace l'ancienne position da la cabine par la nouvelle.

        """
        self.position = new_pos
    
    def get_acceleration(self) :
        """
        Fonction membre qui calcule l'accéleration la cabine.
        Paramètre : 
            self : l'objet courant de type cabine 
        Résultat : 
            La fonction retourne un flottant qui représente l'accéleration actuelle.

        """
        return -self.G*self.direction[2] - (self.coeff_frottement/self.masse)*self.vitesse 
    
    
    
    def update(self, dt):
        '''
        Fonction mettant à jour les attributs de la cabine (position, vitesse..),
        au temps t + dt.
        ---------
        Paramètres : l'intervalle de temps de dt pendant lequel la cabine se déplace
        --------
        Résultat : Ne renvoie rien
        
        '''
        n = len(self.circuit)
        k = self.indice_circuit
        
        ###################################
        
        
        if self.vitesse < 0:  
            #Si la vitesse est négative on va vérifier si le point suivant est
            #avant ou après le point k du circuit
            suivant = k
            suivant2 = k - 1
            
            #gère les cas où l'on est en bout de circuit
            if k == 0:
                suivant2 = n - 1
        
        else:
            #Si la vitesse est positive on vérifie par rapport au point k+1 du
            #circuit
            suivant = k+1
            suivant2 = k+2
            
            if k == n - 1:
                suivant = 0
                suivant2 = 1
                
            elif k == n - 2 :
                suivant2 = 0
        ####################################
        
        p_t = self.position
        u_t = self.direction
        v_t = self.vitesse
        acceleration = self.get_acceleration()
        
        d = v_t*dt
        
        dks = distance3D(p_t, self.circuit[suivant])
        distance_eval = distance3D(p_t, p_t + d*u_t) - dks 
        #Si distance_eval > 0 on a dépassé le point suivant dans le circuit, on
        #doit donc calculer le temps que l'on met à arriver au point suivant et le 
        #temps qu'il reste (dt2)
       
        #Condition d'arrêt
        if distance_eval < 0 :
            v_t = v_t + acceleration*dt
            self.set_v(v_t)
            self.set_pos(p_t + d*u_t)
            
        
        else : 
            dt1 = dks/v_t
            dt2 = dt - dt1
            
            #ARRIVEE AU PT K+1 : on met à jour vitesse et position
            v_t = v_t + acceleration*dt1
            self.set_v(v_t)
            self.set_pos(self.circuit[suivant])
            
            #APRES LE PT K+1 : on met à jour la direction 
            u_new = self.circuit[suivant2] - self.circuit[suivant]
            
            if v_t >= 0 :
                self.indice_circuit = suivant
            else :
                self.indice_circuit = suivant2
                u_new = -u_new 
            u_t = u_new/np.linalg.norm(u_new)
            self.set_u(u_t)
            
            self.update(dt2)



"""______________________________________TESTS___________________________________________"""
if __name__ == '__main__' :
    
    plt.close('all')
    #S = Spline([(0,0,0), (3,0,-6), (2,2,2), (1.5,2,4)])
    S = Spline([(0,0,0), (3,0,-6), (2,2,2), (1.5,2,4)])
    pts_finaux = S.genererPoints(100)
    
    x = np.hstack((pts_finaux[:,0],S.pointsControle[0,0]))
    y = np.hstack((pts_finaux[:,1],S.pointsControle[0,1]))
    z = np.hstack((pts_finaux[:,2],S.pointsControle[0,2]))
    C = Cabine(100000, 200, 0.1, S)
    
    
    '''pos_cabine = np.array([S.pointsControle[0]])
    for t in T :
        print('\n')
        print("t = ", t)
        
        C.update(dt)
        pos_cabine = np.vstack((pos_cabine, C.position))'''
        
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    pointAnime, = ax.plot([], [], [], marker="o", color = 'blue')
    '''
    fig.set_facecolor('black')
    
    ax.set_facecolor('black')
    ax.grid(False)
    ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
    ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
    ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))'''
    
    ax.plot(x, y, z, color = 'black')
    ax.scatter(S.pointsControle[:,0], S.pointsControle[:,1], S.pointsControle[:,2], color = 'red')
    #ax.scatter(x,y,z, marker = '.', color = 'blue')
    
    move = 0.5
    
    def maFonctionDanimation1(t):
        #print('frame = ', t)
        C.update(t)
        #print('\n')
        #print("t = ", t)
        pos_x = C.position[0]
        pos_y = C.position[1]
        pos_z = C.position[2]
        
        #ax.set_xlim((pos_x - move, pos_x + move))
        #ax.set_ylim((pos_y - move, pos_y + move))
        #ax.set_zlim((pos_z - move, pos_z + move))
        
        pointAnime.set_data_3d(pos_x,pos_y,pos_z)
        #n = [C.direction(), 
        #ax.view_init(90,0)
        
        return pointAnime,
        
    dt = 0.01
    T = dt*np.ones(100)
    monAnimation = animation.FuncAnimation(fig=fig, func=maFonctionDanimation1, frames= T, interval=40, blit=False)

    plt.show()
        