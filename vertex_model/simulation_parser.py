from time import time
import dill
import os
import numpy as np
import pandas as pd
from pathlib import Path
import re
import collections
#need to list all files in simulations folder
def remove_values_from_list(the_list, val):
   return [value for value in the_list if value != val]
def flatten(t):
    return [item for sublist in t for item in sublist]
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
    if filter==0:
        simulations_model_1=all_simulations
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

#new function for requesting file indexes based on files in directory
#eg order the indexes of files with model 1 k =100 d=0.01 etc
#%%
def index_of_simulation(model_type="none",k="none",d="none",r="none",run="none",duration="none"): #acceptable values are "all" or a float/integer 
    query={"model_type":str(model_type),"k":str(k),"d":str(d),"r":str(r),"run":str(run),"duration":str(duration)}
    for ind in list(query):
        if query[ind]=="none":
            query.pop(ind)
            continue
        if((query[ind] != "all") and # %%
isinstance(int(float(query[ind])), int)!=True and # %%
isinstance(float(int(query[ind])), float)!=True):
            print(f"please specify {ind}  as either \"all\" or an integer")
            return 
        
    #now that parameters have passed check we first list all files
    simulations_dict=make_simulations_dic(filter=0)
    metadata_dic=make_metadata_dic(simulations_dict)
    metadictionary_list=[]
    for i in range(0,len(metadata_dic)):
        test_file=metadata_dic[i]
        file_model=test_file[0:7][-1]
        file_k=re.search(r'_k_(.*)_D_',test_file).group(1)
        file_D=re.search(r'_D_(.*)_r_',test_file).group(1)   
        file_r=re.search(r'_r_(.*)_run_',test_file).group(1)
        file_run=re.search(r'run_(.*)_sim-duration_',test_file).group(1)
        file_duration=re.search(r"_sim-duration_(.*).pkl",test_file).group(1)
        simulation_metadictionary={}
        simulation_metadictionary={"name":test_file}
        simulation_metadictionary["model_type"]=file_model
        simulation_metadictionary["index"]=i
        simulation_metadictionary["k"]=file_k
        simulation_metadictionary["d"]=file_D
        simulation_metadictionary["r"]=file_r
        simulation_metadictionary["run"]=file_run
        simulation_metadictionary["duration"]=file_duration
        metadictionary_list.append(simulation_metadictionary)
    matches_index=[]
    matches_filename=[]
    for key_query in list(query):
       
        for dictionary in metadictionary_list:
            
            for key_dict in list(dictionary):
                
                if key_query==key_dict:
                    if str(float(query[key_query]))==str(float(dictionary[key_dict])):
                        matches_index.append(dictionary["index"])
                        matches_filename.append(dictionary["name"])
    threshold_for_count=len(list(query))
    b=[]
    print("look here aaryan. 4")
    for item,count in collections.Counter(matches_index).items():
        if count > threshold_for_count:
            b.append(item)
        if count < threshold_for_count:
            matches_index=remove_values_from_list(matches_index,item)
    if b!=[]:
        matches_index=b
    b=[]
    for item,count in collections.Counter(matches_filename).items():
        if count > threshold_for_count:
            b.append(item)
        if count < threshold_for_count:
            matches_filename=remove_values_from_list(matches_filename,item)
    if b!=[]:
        matches_filename=b
    matches_index=list(dict.fromkeys(matches_index))
    matches_filename=list(dict.fromkeys(matches_filename))        
    return [matches_index,matches_filename]
    #need to take a file and parse its k and stuff
    # %%

#search time using query as the input
#first make a dictionary of all of metadata_dic


