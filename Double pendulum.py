# -*- coding: utf-8 -*-
"""
Created on Sat Jul  7 18:17:09 2018

@author: Rishi Acharya
"""

import math
import matplotlib

G = 9.80665
DURATION = 20
STEP = 0.0001

L_1 = 1
L_2 = 1
M_1 = 1
M_2 = 1
THETA_IN = 0.5
PHI_IN = 0.5
THETA_DOT_IN = 0
PHI_DOT_IN = 0


def theta_2dot(l_1:float, l_2:float,
               m_1:float, m_2:float,
               theta:float, phi:float,
               theta_dot:float, phi_dot:float
               )->float:
    n_1 = m_2 * G * math.sin(theta)
    n_2 = m_2 * l_1 * (theta_dot**2) * math.sin(theta-phi)
    n_3 = m_2 * l_2 * phi_dot *\
          (theta_dot-phi_dot) * math.sin(theta-phi)
    n_4 = m_2 * l_2 * theta_dot * phi_dot * math.sin(theta - phi)
    n_5 = m_1 * G * math.sin(theta)
    n_6 = m_2 * G * math.sin(theta)
    d_1 = m_1 * L_1
    d_2 = m_2 * L_1
    d_3 = m_2 * l_1 * math.cos(theta - phi)
    out = (n_1 - n_2 + n_3 - n_4 - n_5 - n_6)/(d_1 + d_2 - d_3)
    return out


def phi_2dot(l_1:float, l_2:float,
            m_1:float, m_2:float,
            theta:float, phi:float,
            theta_dot:float, phi_dot:float
            )->float:
    n_1 = l_1 * theta_dot * phi_dot * math.sin(theta - phi)
    n_2 = G * math.sin(phi)
    n_3 = l_1 * theta_2dot(l_1, l_2, m_1, m_2,
                           theta, phi, theta_dot, phi_dot)\
                           * math.cos(theta - phi)
    n_4 = l_1 * theta_dot * (theta_dot-phi_dot) * math.sin(theta - phi)
    out = (n_1 - n_2 - n_3 + n_4)/l_2
    return out

def get_energy(l_1:float, l_2:float,
               m_1:float, m_2:float,
               theta:float, phi:float,
               theta_dot:float, phi_dot:float
               )->float:
    e_k = (0.5 * m_1 * (l_1**2) * (theta_dot**2))\
        + (0.5 * m_2 * ((theta_dot*l_1*math.cos(theta))
           + (phi_dot*l_2*math.cos(phi))\
           - (theta_dot*l_1*math.sin(theta))\
           - (phi_dot*l_2*math.sin(phi))))
    e_p = (m_1*G*l_1)\
         -(m_1*G*l_1*math.cos(theta))\
         +(m_2*G*l_1)\
         +(m_2*G*l_2)\
         -(m_2*G*l_1*math.cos(theta))\
         -(m_2*G*l_2*math.cos(phi))
    e_t = e_k + e_p
    return e_t

E_T = get_energy(L_1, L_2, M_1, M_2,
                 THETA_IN, PHI_IN, THETA_DOT_IN, PHI_DOT_IN)

t = 0

thetas = []
theta_dots = []
phis = []
phi_dots = []
times = []

thetas.append(THETA_IN)
phis.append(PHI_IN)
theta_dots.append(THETA_DOT_IN)
phi_dots.append(PHI_DOT_IN)
times.append(t)

print("Start loop")

while t <= DURATION:
    theta_dot = theta_dots[-1] + (theta_2dot(L_1, L_2, M_1, M_2,
                           thetas[-1], phis[-1], theta_dots[-1],
                           phi_dots[-1])\
                           * STEP)
    theta = thetas[-1] + theta_dot * STEP
    phi_dot = phi_dots[-1] + (phi_2dot(L_1, L_2, M_1, M_2,
                           thetas[-1], phis[-1], theta_dots[-1],
                           phi_dots[-1])\
                           * STEP)
    phi = phis[-1] + phi_dot * STEP
    theta_dots.append(theta_dot)
    thetas.append(theta)
    phi_dots.append(phi_dot)
    phis.append(phi)
    times.append(t)
    t += STEP

energy_change = get_energy(L_1, L_2, M_1, M_2,
                           thetas[-1], phis[-1],
                           theta_dots[-1], phi_dots[-1])\
               - get_energy(L_1, L_2, M_1, M_2,
                           thetas[0], phis[0],
                           theta_dots[0], phi_dots[0])

print(energy_change)

fig_theta = matplotlib.pyplot.figure()
axes_theta = fig_theta.add_subplot(111)
axes_theta.plot(thetas, theta_dots, c='tab:blue')
axes_theta.grid()
axes_theta.set(
        xlabel='Theta', ylabel='Theta dot',
        title='Phase portrait- theta'
        )
matplotlib.pyplot.savefig('double pendulum phase portrait theta.jpg')

fig_phi = matplotlib.pyplot.figure()
axes_phi = fig_phi.add_subplot(111)
axes_phi.plot(phis, phi_dots, c='tab:red')
axes_phi.grid()
axes_phi.set(
        xlabel='Phi', ylabel='phi dot',
        title='Phase portrait- phi'
        )
matplotlib.pyplot.savefig('double pendulum phase portrait phi.jpg')

fig_angles = matplotlib.pyplot.figure()
axes_angles = fig_angles.add_subplot(111)
axes_angles.plot(thetas, phis, c='tab:green')
axes_angles.grid()
axes_angles.set(
        xlabel='Theta', ylabel='Phi',
        title='Angles'
        )
matplotlib.pyplot.savefig('double pendulum angles.jpg')

fig_theta_time = matplotlib.pyplot.figure()
axes_theta_time = fig_theta_time.add_subplot(111)
axes_theta_time.plot(times, thetas, c='tab:orange')
axes_theta_time.grid()
axes_theta_time.set(
        xlabel='Time', ylabel='Theta',
        title='Time theta'
        )
matplotlib.pyplot.savefig('double pendulum time theta.jpg')

fig_phi_time = matplotlib.pyplot.figure()
axes_phi_time = fig_phi_time.add_subplot(111)
axes_phi_time.plot(times, phis, c='tab:cyan')
axes_phi_time.grid()
axes_phi_time.set(
        xlabel='Time', ylabel='Phi',
        title='Time phi'
        )
matplotlib.pyplot.savefig('double pendulum time phi.jpg')