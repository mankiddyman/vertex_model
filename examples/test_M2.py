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
import vertex_model.timer as timer
import sys

import time
timestart = time.time()

type_=5

# parameters of the vertex model
G=0.075 #Cortical Contractility
L=0.04 # Line tension
K=1.0# Elastic energy coefficient

# parameters of the nucleus A-B stochastic dynamics
# sys.argv to pass from command line
# call as, e.g.: python test_M1.py 100 0.15
k=float(sys.argv[1]) # 100
D=float(sys.argv[2]) # 0.15
r=float(sys.argv[3]) # size of crowding force relative to inertia constant #1
no_simulations=float(sys.argv[4]) #number of simulations to run #1
duration=float(sys.argv[5]) #duration of simulation
# parameters of the crowding force
s=0.2
a=r*k

# adding the nuclear movement parameters to 'params'
params = [K,G,L,k,D,s,a]

np.random.seed(1999)
rand = np.random.RandomState(1999)

#Cell_cycle_lengths_progenitors = []
#Cell_cycle_lengths_pois = []
#final_history = []
for i in range(no_simulations):
    history = run_simulation_INM(params,duration,rand,type_)  # timend = 1 for test
    # Rows 32-49 -> copy in ANALYSIS
    """
    last = history[-1] # get state of the cells object from the previous step -> this + following rows -> analysis
    gone = last.empty()
    divided = gone
    diffntd = gone
    if 'poisoned' in last.properties:
	    pois = last.properties['poisoned']==1
	    healty = last.properties['poisoned']==0
	    divided = gone & healty
	    diffntd = gone & pois
    progenitor_deadcells_last = np.where(divided)[0]
    differentiated_last = np.where(diffntd)[0]
    correct_cell_cycle_lengths = last.properties['age'][progenitor_deadcells_last]
    pois_cell_cycle_lengths = last.properties['age'][differentiated_last]
    
    #new_cell_cycle_lengths = last.properties['age'][deadcells_last] # cell cycle lengths = age of deadcells
    Cell_cycle_lengths_progenitors = np.append(Cell_cycle_lengths_progenitors, correct_cell_cycle_lengths)
    Cell_cycle_lengths_pois = np.append(Cell_cycle_lengths_pois, pois_cell_cycle_lengths)
    final_history = np.append(final_history, history) # get a final history array

    print(last.properties['nucl_pos'])
    print(min(last.properties['nucl_pos']))
    print(max(last.properties['nucl_pos']))
    print(last.properties['A0_final'])
    print(len(last.properties['final_nucl_pos']))
    """

timeend = time.time()
timer.timer(timestart,timeend)

import dill
with open (f"simulations/model_2_k{k}_D{D}_run_{run_id}.pkl", "wb") as file:
    dill.dump(history, file)



