#%%
import ipdb
import itertools
# get_ipython().magic(u'matplotlib inline')
import numpy as np
import matplotlib
# %matplotlib qt5
import matplotlib.pyplot as plt
from matplotlib import animation, rc
# from IPython.display import HTML
import vertex_model as model
from vertex_model.run_select import run_simulation_INM, definecolors
import vertex_model.initialisation as init
import sys
import vertex_model.timer as timer
import time
import dill



# parameters of the nucleus A-B stochastic dynamics
# sys.argv to pass from command line
# call as, e.g.: python test_M1.py 100 0.15

# parameters of the vertex model

G=0.075 #Cortical contractility
L=0.04 # Line tension
K=1.0 #elastic energy coefficient

#%%
# adding the nuclear movement parameters to 'params'

#sys.argv=['x',0, 0, 29, 29, 0.0, 0.1, 1, 5, 10]

try:
    k=sys.argv[1:5]
    k = [ int(x) for x in k ]
     # inertia array 
    D=float(sys.argv[5]) # Browian diffusion constant #0.15
    a=float(sys.argv[6]) # crowding force strength absolute
    no_simulations=float(sys.argv[7]) #number of simulations to run #1
    duration=float(sys.argv[8])#duration of simulation
    type_=int(sys.argv[9]) 
except:
    print("Missing command line parameters. Execute as\n")
    print("   python test_M1.py  <k_g1> <k_s> <k_g2> <k_m>  <D>  <a> <no_of_simulations> <duration_of_simulation> <sim_type> \n")
    exit()

#%%
# k=[12,0,24,24]#g1,s,g2,m
# D=0.00001
# a=0
# no_simulations=1
# duration=306
s=0.1
params = [K,G,L,k,D,s,a]
# np.random.seed(1999)
rand = np.random.RandomState(1999)
x=params
timend=duration
rand=rand
sim_type=type_
# lifespan=100.0
#%%

def run_sim (run_id):
    timestart = time.time()
    history = run_simulation_INM(params,duration,rand,type_)  # timend = 1 for test
    #params=x
    #timeend=duration
    #rand=rand
    #type=sim_type
    timeend = time.time()
    timer.timer(timestart,timeend)
    with open (f"simulations/model_active_{sim_type}_k_{k}_D_{D}_a_{a}_run_{run_id}_sim-duration_{duration}.pkl", "wb") as file:
        dill.dump(history, file)

    return timeend - timestart


for run in range(1,int(no_simulations+1)):
	run_sim(run)


# %%
