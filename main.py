# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 23:57:31 2022

@author: adrie
"""

import matplotlib.pyplot as plt
from affichage import *
import matplotlib.animation as animation

import mpl_toolkits.mplot3d.axes3d as p3

#creer_fichier('test')

S,C = creer_spline_cabine('test.txt')

plt.close('all')
figure = plt.figure()
ax, ptAnime = afficherCircuit(S, 100, figure)


maFonctionDanimation = afficherTrajectoire(C, ax, ptAnime, 0.5, camera = 2)

dt = 0.01
T = dt*np.ones(2)
monAnimation = animation.FuncAnimation(fig=figure, func=maFonctionDanimation, frames=T, interval=50, blit=False)
plt.show()