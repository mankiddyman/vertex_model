from time import time
import dill
import os
import numpy as np
import pandas as pd
from pathlib import Path
import re
import collections
import numpy
import sys
from ast import literal_eval
#need to list all files in simulations folder
def remove_values_from_list(the_list, val):
   return [value for value in the_list if value != val]
def flatten(t):
    return [item for sublist in t for item in sublist]
#%%
#get ids of neighbours for each cell
def neighbour_finder(sim_object,index_of_interest):
    a=sim_object
    cell_ids = a.mesh.face_id_by_edge
    neigh_ids=cell_ids[a.mesh.edges.reverse]
    n_neigbours=len(np.where(cell_ids==index_of_interest)[0])
    
    if n_neigbours!=0:
        cell_edge_indices=np.where(cell_ids==index_of_interest)[0]
        cell_reverse_edge_indices=a.mesh.edges.reverse[cell_edge_indices]
        neighbour_indices=cell_ids[cell_reverse_edge_indices]
    
    else:
         neighbour_indices=[]
    return neighbour_indices






#%%
def make_metadata_dic():
    vertex_model_path=Path(os.getcwd()).parent.absolute()
    simulations_path=os.path.join(vertex_model_path,"examples","simulations")
    all_simulations=os.listdir(simulations_path)
    a=[]
    for i in range(0,len(all_simulations)):
        if all_simulations[i].endswith(".pkl"):
            a.append(all_simulations[i])
    return a
#%%
def slicedict(d, s):
    return {k:v for k,v in d.items() if k.startswith(s)}
#%%
def make_simulations_dic(metadata_index):
    simulation=make_metadata_dic()[metadata_index]
#now opening files and appending to dataframe
    simulations_dict={}
    with open(os.path.join("simulations",simulation),"rb") as file:
        history=dill.load(file)
        for i in range(0,len(history),1):
            simulations_dict[str(simulation)+"_timepoint_"+str(i+1)+"_of_"+str(len(history))]=history[i]

    return simulations_dict
#simulations_dict=make_simulations_dic(filter=0)
#metadata_dic=make_metadata_dic(simulations_dict)
#new function for requesting file indexes based on files in directory
#eg order the indexes of files with model 1 k =100 d=0.01 etc
#%%
def index_of_simulation(model_type="none",k_g1="none",k_g2="none",k_s="none",k_m="none",d="none",a="none",run="none",duration="none"): 
    #acceptable values are "all" or a float/integer 
    query={"model_type":str(model_type),"k_g1":str(k_g1),"k_s":str(k_s),"k_g2":str(k_g2),"k_m":str(k_m),"d":str(d),"a":str(a),"run":str(run),"duration":str(duration)}
    for ind in list(query):
        if query[ind]=="none" or query[ind]=="all":
            query.pop(ind)
            continue
        if query[ind]=="active":
            query[ind]=3
            continue
        
    #now that parameters have passed check we first list all files
   
    metadata_dic=make_metadata_dic()
    metadictionary_list=[]
    for i in range(0,len(metadata_dic)):
        test_file=metadata_dic[i]
        file_model=re.search(r"model_(.*)_k",test_file).group(1)
        if file_model=="active":
            file_model=3
        file_k=re.search(r'_k_(.*)_D_',test_file).group(1)
        file_k=file_k[1:-1].split(",")
        file_k_g1=file_k[0]
        file_k_g2=file_k[2]
        file_k_s=file_k[1]
        file_k_m=file_k[3]
        file_D=re.search(r'_D_(.*)_a_',test_file).group(1)   
        file_a=re.search(r'_a_(.*)_run_',test_file).group(1)
        file_run=re.search(r'run_(.*)_sim-duration_',test_file).group(1)
        file_duration=re.search(r"_sim-duration_(.*).pkl",test_file).group(1)
        simulation_metadictionary={}
        simulation_metadictionary={"name":test_file}
        simulation_metadictionary["model_type"]=file_model
        simulation_metadictionary["index"]=i
        simulation_metadictionary["k"]=file_k
        simulation_metadictionary['k_g1']=file_k_g1
        simulation_metadictionary['k_s']=file_k_s
        simulation_metadictionary['k_g2']=file_k_g2
        simulation_metadictionary['k_m']=file_k_m
        simulation_metadictionary["d"]=file_D
        simulation_metadictionary["a"]=file_a
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
                        #continue
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
def get_time(input):
    output=re.search(r"timepoint.*$",str(input))[0]
    return output
#now turning a simulation into a csv

#test index
#%%
def make_df(metadata_index):#NB pass only one experiment at a time
    experiment_metadata=make_metadata_dic()[metadata_index]
    new_dict=make_simulations_dic(metadata_index)
    df_master=pd.DataFrame()
    for i in range(0,len(new_dict.values())):
        a=list(new_dict.values())[i] #the dictionary object for the i'th timepoint
#will first try getting a df
        numRows=a.__len__()
        df=pd.DataFrame(index=range(numRows))
#first collumn should be file name 
        df['file_index']=[metadata_index]*numRows
        df['file_name']=list(new_dict)[i]
#need function to add a collumn for timepoints
        df['time_point']=df["file_name"].apply(get_time)
        df['time']=[int(re.search(r"timepoint_(.*)_of_",i).group(1)) for i in df['time_point']]
        df['cell_index']=range(0,numRows)
        df['parent']=a.properties['parent']
        df["K"]=[a.properties['K']]*numRows
        if df['file_name'][0].startswith("model_active"):
            df['k(t)']=a.properties['k(t)']
        df['Gamma']=[a.properties['Gamma']]*numRows
        df['Lambda']=[a.properties['Lambda']]*numRows
        df['age']=a.properties['age']
        df['dead']=a.empty()
        df['parent']=a.properties['parent']
        #df['parent_group']=a.properties['parent_group'] #have to figure out what this means
        df['ageingrate']=a.properties['ageingrate']

        #df['ids_division']=a.properties['ids_division']
        #df['force_x']=a.properties['force_x']
        #df['force_y']=a.properties['force_y']
        df['zposn']=a.properties['zposn']
        df['A0']=a.properties['A0']
        df['k']=[a.properties['k']]*numRows
        df['k_g1']=a.properties['k'][0]
        df['D']=a.properties['D']
        df['nucl_pos']=a.properties['nucl_pos']
        df["apical_area"]=a.mesh.area
        if experiment_metadata.startswith("model_2"):
            df['r']=[a.properties['a']/a.properties['k']]*numRows
        elif experiment_metadata.startswith("model_1"):
            df['r']=0
        elif experiment_metadata.startswith("model_active"):
            df['a']=[a.properties['a']]*numRows
        df['x_centre']=list(a.mesh.centres()[0])
        df['y_centre']=list(a.mesh.centres()[1])
        df['neighbour_ids']=0*numRows 
        temp_list=[]
        for j in range(0,numRows):
            temp_list.append(neighbour_finder(a,j))
        df['neighbour_ids']=temp_list
        for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'K', 'Gamma','Lambda', 'dead', 'parent', 'A0','D','x_centre','y_centre']:
            
            df[col]=df[col].astype('category')
        df['neighbour_ids']=df['neighbour_ids'].apply(list)
        df['n_neighbours']=df['neighbour_ids'].apply(len)
        df_master=pd.concat([df_master,df])
   
        del df
    filepath=os.path.join("simulations",experiment_metadata[:-3]+"csv")
    df_master.to_csv(filepath,index=False)
    return df_master

def retrieve_sim_metadata(sim_name): #function to be used in make_simulation_metadataframe() where u pass a simulation name and it returns each parameter
    #retrieving each parameter from string parsing of .pkl file
    file_name=sim_name
    model_type=re.search(r"model_(.*)_k",sim_name).group(1)
    k_list=literal_eval(re.search(r"k_(.*)_D",sim_name).group(1))
    k_g1=k_list[0]
    k_s=k_list[1]
    k_g2=k_list[2]
    k_m=k_list[3]
    D=float(re.search(r"D_(.*)_a",sim_name).group(1))
    a=float(re.search(r"a_(.*)_run",sim_name).group(1))
    run=int(re.search(r"run_(.*)_sim",sim_name).group(1))
    duration=float(re.search(r"_sim-duration_(.*).pkl",sim_name).group(1))
    return [file_name,model_type,k_g1,k_s,k_g2,k_m,D,a,run,duration]

def make_simulation_metadataframe(): #will make a dataframe containing simulation metadata to retrieve particular simulations with ease.\
    #first list every simulation
    all_sims=make_metadata_dic()
    #now running retrieve_sim_metadata on each file
    parameters=[]
    for i in range(0,len(all_sims)):
        parameters.append(retrieve_sim_metadata(all_sims[i]))
    df=pd.DataFrame(parameters)
    df.columns=['file_name','model_type','k_g1','k_s','k_g2','k_m','D','a','run','duration']
    
    
    return(df)
#%%

