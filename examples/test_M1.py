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


type_=4

# parameters of the vertex model
G=0.075 #Cortical contractility
L=0.04 # Line tension
K=1.0 #elastic energy coefficient

# parameters of the nucleus A-B stochastic dynamics
# sys.argv to pass from command line
# call as, e.g.: python test_M1.py 100 0.15
try:
	k=float(sys.argv[1]) # inertia contribution to IKNM #example 100
	D=float(sys.argv[2]) # Browian diffusion constant #0.15
	r=float(sys.argv[3]) # size of crowding force relative to inertia constant #1
	no_simulations=float(sys.argv[4]) #number of simulations to run #1
	duration=float(sys.argv[5]) #duration of simulation 
except:
	print("Missing command line parameters. Execute as\n")
	print("   python test_M1.py  <k>  <D>  <r = a/k> <no_of_simulations> <duration_of_simulation> \n")
	exit()

# parameters of the crowding force
s=0.2
a=r*k

# adding the nuclear movement parameters to 'params'
params = [K,G,L,k,D,s,a]

def run_sim (run_id):
    timestart = time.time()
    history = run_simulation_INM(params,duration,rand,type_)  # timend = 1 for test
    timeend = time.time()
    timer.timer(timestart,timeend)
    with open (f"simulations/model_1_k{k}_D{D}_run_{run_id}.pkl", "wb") as file:
        dill.dump(history, file)

    return timeend - timestart

np.random.seed(1999)
rand = np.random.RandomState(1999)

for run in range(1,int(no_simulations+1)):
	run_sim(run)




