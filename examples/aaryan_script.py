#setup
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation, rc
#from IPython.display import HTML
import vertex_model as model
from vertex_model.run_select import run_simulation_INM, definecolors
import vertex_model.initialisation as init

#now
#### Parameter used in the simulation are picked up form Gobal_Constant.py file

### Choose simulation to run
# def run_simulation_INM(x, timend,rand, type):
#     sim_type 0 simulation_with_division_clone (no differentiation rate)
#     sim_type 1 simulation_with_division_clone_differentiation (all differentiation rate)
#     sim_type 2 simulation_with_division_clone_differenciation_3stripes (2 population with and without diffentiation rate)
#     sim_type 3 simulation_with_division_clone_whole_tissue_differenciation (differentiation rate everywhere)

#use type 2 and parameter i=4 for floor plate
#run simulation with the choosen parameters ; line 521 in run_select.py

rand =  np.random.RandomState() #random number to choose Lambda
params_dict = ["K":1,"G":0.075,"L":0.05,]  # K=x[0],G=x[1],L=x[2]
history= run_simulation_INM(params,300,rand,type_) #return hist; simulation time = 300?
