#%%
from cProfile import label
from re import L
from IPython.display import HTML
from curses import meta
from enum import unique
from fileinput import filename
import itertools
from time import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import vertex_model as model
import vertex_model.initialisation as init
from vertex_model.forces import TargetArea, Tension, Perimeter, Pressure
from matplotlib import animation, rc
import random
import dill
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import os
from vertex_model import simulation_parser
from vertex_model.simulation_parser import *
import matplotlib
import gc
from ast import literal_eval
import scipy
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
from vertex_model import Gobal_Constant
from vertex_model.Gobal_Constant import dt, viscosity, t_G1, t_G2, t_S, A_c, J, pos_d, T1_eps, P, microns, time_hours, expansion_constant #file wit
from vertex_model import script_for_memory_issues
from vertex_model.script_for_memory_issues import *
from vertex_model.run_select import *
matplotlib.rc('font', family='sansserif', serif='cm10')
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
random.seed(123)
def save_multi_image(filename):
    pp = PdfPages(filename)
    fig_nums = plt.get_fignums()
    figs = [plt.figure(n) for n in fig_nums]
    for fig in figs:
        fig.savefig(pp, format='pdf')
    pp.close()
metadata_dic=simulation_parser.make_metadata_dic()



#first figure is effect of area on E
#%%
A0=0.5 #target
K=1#elasticity
A=np.linspace(0,1,100)
E=K*(A-A0)**2

#plot
fig = plt.figure(dpi=500)
ax = fig.add_subplot(1, 1, 1)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.plot(A,E,label=r"$A_{a}$")
plt.xlabel(r"$A_{a}$ (AU)")
plt.ylabel(r"$E$  (AU)")
plt.axvline(x=A0,color="r",label="Target Area")
plt.legend(loc="lower right")
plt.title("Deviations from the target area raise the energy of the system.")
plt.text(0.9, 0.25,s=f'K = {K}\n$A_{0}$ = {A0}',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.savefig(os.path.join("output","E~A_figure.png"))
# %% E~l
Lambda=0.04
l=np.linspace(0,100,100)
E=Lambda*l
fig = plt.figure(dpi=500)
ax = fig.add_subplot(1, 1, 1)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.plot(l,E)
plt.xlabel(r"$l_{i,j}$ (AU)")
plt.ylabel(r"$E$  (AU)")

plt.title("Increases in the length of an edge slowly but steadily\n increase the energy of the system,")
plt.text(0.9, 0.2,s=rf'$\Lambda$ = {Lambda}',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.savefig(os.path.join("output","E~l_figure.png"))

# %% E~L
Gamma=0.075
L=np.linspace(0,10,100)
E=Gamma/2*L**2
fig = plt.figure(dpi=500)
ax = fig.add_subplot(1, 1, 1)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.plot(L,E)
plt.xlabel(r"$L_{a}$ (AU)")
plt.ylabel(r"$E$  (AU)")

plt.title("The perimeter of the cell increases the energy of the system quadratically.")
plt.text(0.9, 0.2,s=rf'$\Gamma$ = {Gamma}',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)

plt.savefig(os.path.join("output","E~L_capital_figure.png"))


# %% A0~z
z=np.linspace(-0.2,1,100)
A0=1/2*(1+np.maximum(z,0)**2)
Gamma=0.075


fig = plt.figure(dpi=500)
ax = fig.add_subplot(1, 1, 1)

ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.plot(z,A0,label=r"$z_{a}$")

plt.axvline(x=0.75,color="r",label=r"$z_{c}$")
plt.axhline(y=1.3,color="y",label=r"$A_{c}$")
plt.xlabel(r"$z_{a}$ (AU)")
plt.ylabel(r"$A_{0}$  (AU)")
plt.legend(loc="upper left")
plt.title("As the nucleus moves towards the apex,\n apical area increases quadratically.")

plt.savefig(os.path.join("output","A0_z_figure.png"))


# %% A0~age
t=np.linspace(0,2,100)
g_a=1
A0=1/2*(t*g_a+1)


fig = plt.figure(dpi=500)
ax = fig.add_subplot(1, 1, 1)

ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.plot(t,A0,label=r"$t$")

plt.axvline(x=0.89,color="r",label=r"$t_{c}$")
plt.axhline(y=1.3,color="y",label=r"$A_{c}$")
plt.xlabel(r"$t$ (AU)")
plt.ylabel(r"$A_{0}$  (AU)")
plt.legend(loc="lower right")
plt.title("Area increases as a function of a cell's age")
plt.text(0.9, 0.3,s=r'$g_{a}$'+f' = {g_a}',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.savefig(os.path.join("output","A0_t_figure.png"))

# %% z~t
z_0_t=[0,0,1,1]
k_t=[12,0,24,24]
dt=1000
z=[1]*dt

t=np.linspace(0,1,dt)
t_G1=0.4
t_S=0.33
t_G2=0.16
t_M=0.11

for i in range(0,len(t)-1):
    if(0<=t[i]<t_G1):
        k=k_t[0]
        z_0=z_0_t[0]
    elif(t_G1<=t[i]<t_G1+t_S):
        k=k_t[1]
        z_0=z_0_t[1]
    elif(t_G1+t_S<=t[i]<t_G1+t_S+t_G2):
        k=k_t[2]
        z_0=z_0_t[2]
    elif(t_G1+t_S+t_G2<=t[i]):
        k=k_t[3]
        z_0=z_0_t[3]
    
    z[i+1]=z[i]+k*(z_0-z[i])*1/dt


fig = plt.figure(dpi=500)
ax = fig.add_subplot(1, 1,1)
plt.plot(t,z,label=r"$z$")
plt.axhline(y=0.75,color='red',label=r"$z_{c}$")
plt.axvline(x=t_G1,color="black")
plt.axvline(x=t_G1+t_S,color="black")
plt.axvline(x=t_G1+t_S+t_G2,color="black")
plt.axvline(x=t_G1+t_S+t_G2+t_M,color="black")
plt.ylim(0,1)

plt.xlabel(r"$t$ (AU)")
plt.ylabel(r"$z$  (AU)")
plt.legend(loc="lower left")
plt.title("Nuclear dynamics in the purely active model")
plt.text(t_G1/2, 0.5,s='G1',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.text(np.mean([t_G1,t_G1+t_S]), 0.5,s='S',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.text(np.mean([t_G1+t_S,t_G1+t_S+t_G2])-0.02, 0.5,s='G2',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.text(np.mean([t_G1+t_S+t_G2,t_G1+t_S+t_G2+t_M])-0.04, 0.5,s='M',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
          
plt.savefig(os.path.join("output","z_t_figure.png"))

# %% f_c~delta_z
delta_z=np.linspace(-.75,0.75,1000)
a=1
s=0.1
f_c=a*delta_z/s*np.e**(-(delta_z)**2/(2*s**2))



fig = plt.figure(dpi=500)
ax = fig.add_subplot(1, 1, 1)

ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.plot(delta_z,f_c,label=r"$\Delta z$")
plt.xlabel(r"$\Delta z$ (AU)")
plt.ylabel(r"$F_{c}$  (AU)")
plt.legend(loc="lower right")
plt.title("Crowding force is maximum at $\Delta z = s$")
plt.text(0.9, 0.2,s=f"$s$ = {s}\n$a=1$",
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)

plt.savefig(os.path.join("output","f_c~delta_z_figure.png"))
# %%
