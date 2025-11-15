# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 13:36:22 2024

@author: Victor HELIES
"""
"""
This file was modified on 07/11/2024 and the result does not correspond to the orignial BMW120D.
The gearbox ratio have been optimized with geometrical step to improve 0-100kph performance, thus the speed curve doesn't correspond to the report anymore.
"""

import sympy as sy
import numpy as np
import matplotlib.pyplot as plt




# Constantes de la voiture :

L = 4.33 # Longueur de la voiture (m) -> Empattement : 2,69m ; Porte-à-faux avant : 0,77m ; Porte-à-faux arrière : 0,87m
l = 1.77 # Largeur de la voiture (m)
h = 1.42 # Hauteur de la voiture (m)
gas = 0.140 # Garde au sol de la voiture (m)
m_vide = 1558 # Masse à vide de la voiture (Kg) -> 1370kg
delta = 0.0157 # Coefficient de résistance au roulement (N/Kg) -> 0.01N/kg
P_mot = 65000 # Puissance moteur (W) -> 177ch


# Partie aérodynamique :

rho = 1.164 # Densité de l'air (kg/m^3)
S = 2.1 # Surface projetée de la voiture (m^2)
Cx = 0.3392 # Coefficient de traînée de la voiture (sans unité)


# Partie résistance au roulement :

# delta = 10 (déjà défini)
M = m_vide # Masse de la voiture (kg)


# Partie pente :

# M = m_vide (déjà défini)
g = 9.81 # Gravité (m/s^2)
alpha = 0 # Angle de la pente (rad)


# Partie moteur :

eta_trans = 0.92 # Rendement de transmission (sans unité)
# P_mot = 132 kW (déjà défini)


# Calcul de la vitesse maximale du véhicule :

x = sy.Symbol('x')
eq = sy.Eq(1/2 * rho * S * Cx * x**3 + delta * M * x + M * g * np.sin(alpha) * x, eta_trans * P_mot)
sol = sy.solve(eq) # Vitesse de la voiture (km/h)
V = sol[0] * 3.6
print("Vitesse maximale de la BMW Série 1 (F20) 120d :", int(V), "km/h")



    
#constantes du moteur

w_max = 4500 # Régime moteur maximal (tr/min)
w_ralenti = 1000 # Régime moteur au ralenti (tr/min)
liste_w = [w_ralenti + (w_max - w_ralenti) * i / 100 for i in range(101)] # Régime moteur (tr/min)




R_roue = 0.308  # Rayon de la roue (m) : 18 pouces + épaisseur du pneu = 31,4cm
k_diff = 1 # Réducteur au différentiel (sans unité) -> Supposé égal à 1 en ligne droite
k_bdv = 1 # Réducteur à la boîte de vitesse (sans unité) -> Inconnu
C_max = 208 # Couple maximal (Nm)






points = {0:0, 1000:160, 1500:200, 1750:208, 2000:205, 2500:195, 3000:187, 3500:173, 4000:150, 4500:105} #courbe de couple moteur
abscisses = list(points.keys())
pentes = [(points[abscisses[i + 1]] - points[abscisses[i]]) / (abscisses[i + 1] - abscisses[i]) for i in range(len(abscisses) - 1)]

def C(x):
    """Renvoie le couple en fonction des RPM."""
    for i in range(len(abscisses) - 1):
        if abscisses[i] <= x < abscisses[i + 1]:
            return pentes[i] * (x - abscisses[i]) + points[abscisses[i]]
    return 0

liste_w = np.linspace(1000, 4500, 100)  # Régime moteur (tr/min)
plt.plot(liste_w, [C(w) for w in liste_w])
plt.xlabel('Régime moteur (tr/min)')
plt.ylabel('Couple (Nm)')
plt.title('Couple en fonction du régime moteur')
plt.grid()
plt.show()





# Constantes du moteur et du véhicule
#rapports
n_rapports = 4
i_1 = 28
i_n = 3.1
i_tot = i_1/i_n
q = i_tot**(1/(n_rapports-1))
rapports = []
for i in range(n_rapports-1):
    rapports.append(1/(i_n*q**(n_rapports-i)))  # Rapports de vitesse
rapports.append(1/i_n)
w_ralenti = 1000*2*np.pi/60  # Régime moteur au ralenti (rad/s)
regime = [w_ralenti]
dt = 0.0001
t = [0]
t0 = 0
i=0
w = w_ralenti#tr/min->rad/s
tmax = 12# secondes
v = [w_ralenti*R_roue*k_diff*rapports[0]]
t_pilote = 0.3#le temps que le pilote met pour passer ses vitesses

for k_bdv in rapports:
    i+=1
    while w*60/(2*np.pi)<w_max:#on change de vitesse quand le moteur atteind 4300 tr/min (déterminé par lecture graphique)
        if t0>tmax:#sécurité mise en place pour éviter que l'algorithme diverge
            break
        w = w + 15*dt*(eta_trans*C(w*60/(2*np.pi))/(R_roue*k_bdv*k_diff) - 0.5*rho*S*Cx*(R_roue*k_bdv*k_diff*w)**2 - M*delta - M*g*np.sin(alpha))/(M*R_roue*k_bdv*k_diff)#méthode d'euler explicite
        regime.append(w*60/(2*np.pi))#                                       ^ remplacer le 0.2 par un 0.5
        t.append(dt+t0)
        v.append(w*k_bdv*k_diff*R_roue*3.6)#w en rad/s -> v en km/h (remplacer le *2 par un *3.6)
        t0 += dt
        t1 = t0
    while abs(t0-t1)<=t_pilote:#le pilote met un certain temps à changer de vitesse
        if t0>tmax:
            break
        w = w + dt*( - 0.5*rho*S*Cx*(R_roue*k_bdv*k_diff*w)**2 - M*delta - M*g*np.sin(alpha))/(M*R_roue*k_bdv*k_diff)#pas de couple moteur ici, car le pilote a débrayé
        regime.append(w)
        t.append(dt+t0)
        t0 += dt
        v.append(w*k_bdv*k_diff*R_roue*3.6)#(remplacer le *2 par un *3.6)
    if i<len(rapports):
        w = w*k_bdv/rapports[i]#pour avoir la vitesse du bolide
    else:
        w = w_ralenti

        
# Tracé du régime moteur en fonction du temps
plt.plot(t, regime)
plt.xlabel('Temps (s)')
plt.ylabel('Régime moteur (tr/min)')
plt.title('Régime moteur en entrée de boite en fonction du temps')
plt.grid()
plt.show()

plt.figure()
plt.plot(t, v)
plt.xlabel('Temps (s)')
plt.ylabel('Vitesse (km/h)')
plt.title('Vitesse en fonction du temps')
plt.grid()
plt.show()
 
