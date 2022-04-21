from time import time
import dill
import os
import numpy as np
import pandas as pd
from pathlib import Path
#need to list all files in simulations folder
def make_simulations_dic(filter):
    vertex_model_path=Path(os.getcwd()).parent.absolute()
    simulations_path=os.path.join(vertex_model_path,"examples","simulations")
    all_simulations=os.listdir(simulations_path)
    simulations_model_1=[]
#filter list for only model 1
    if filter==1:
        for simulation in all_simulations:
            if simulation.startswith("model_1"):
                simulations_model_1.append(simulation)
#filter list for only model 2
    if filter==2:
        for simulation in all_simulations:
            if simulation.startswith("model_2"):
                simulations_model_1.append(simulation)
#now opening files and appending to dataframe
    simulations_dict={}
    for simulation in simulations_model_1:
        with open(os.path.join("simulations",simulation),"rb") as file:
            history=dill.load(file)
            for i in range(0,len(history),1):
                simulations_dict[str(simulation)+"_timepoint_"+str(i+1)+"_of_"+str(len(history))]=history[i]
    return simulations_dict

def make_metadata_dic(dictionary):
    keys=dictionary.keys()
    result=[]
    for i in keys:
        temp = i.split(".pkl")
        result.append(temp[0] + ".pkl")
    metadata_dic=np.unique(result)   
    return metadata_dic
def slicedict(d, s):
    return {k:v for k,v in d.items() if k.startswith(s)}