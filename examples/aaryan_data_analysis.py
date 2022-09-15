#%%
from pickle import FALSE
from collections import Counter
from nis import cat
import math
from re import A
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
from textwrap import wrap
from itertools import chain
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
sns.set(rc={"figure.dpi":500, 'savefig.dpi':500})
t_G1=0.4
t_S=0.33
t_G2=0.16
t_M=0.11
metadataframe=make_simulation_metadataframe()
palette=['b','orange','green','red','purple']
#%%
#function for 3 histograms
def plot_3_histograms(experiment_names,bins_list,seperator,constant):
    #will make a massive dataframe for all 3? sure
    df_master=pd.DataFrame()#contains unfiltered data from all 3 experiments
    df=pd.DataFrame()#empty datafrmae
    collumns=['file_index','file_name',"time_point","cell_index",'time',"dead","k_g1","D",'a',"nucl_pos","apical_area","age","neighbour_ids","n_neighbours"]

    for i in experiment_names:
        a=os.path.join("simulations",i[:-3]+"csv")
        df=pd.read_csv(a,usecols=collumns)
        df_master=pd.concat([df_master,df])
    del df
    for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead', 'k_g1','D','time','a',"n_neighbours"]:
            df_master[col]=df_master[col].astype('category')
    df_master['dead']=df_master['dead'].astype('bool')    

    palette=['b','orange','green','red','purple']
    #applying thermalised cells filter
    max_thermal_cell_index=max(df_master.loc[(df_master.time==1)].cell_index)
    df_master=df_master.query(f'cell_index > {max_thermal_cell_index}')

#now histogram of Area segregated by seperator, additionally want labels of critical area.
    
    #first need to get only alive cells at final timepoint
    df_alive_final_timepoint=df_master.loc[(df_master.dead==False)&(df_master.time==306)]
    #now i want to collect A seperated by seperator
    seperator_list=list(df_alive_final_timepoint[seperator].cat.categories)
    x=[[],[],[]]
    for i in range(0,len(seperator_list)):
        x[i]=(df_alive_final_timepoint.loc[(df_alive_final_timepoint[seperator]==seperator_list[i])]['apical_area'])
    
    kwargs = dict(bins=bins_list[0],histtype='step',density=True)
    fig1, ax=plt.subplots(dpi=500)
    plt.subplots_adjust(left=0.2,bottom=0.2, top = 0.9, right = 0.9)
    plt.hist(x[0], **kwargs, color='g', label=f'${seperator}$={seperator_list[0]}, $n={len(x[0])}$')
    plt.hist(x[1], **kwargs, color='b', label=f'${seperator}$={seperator_list[1]}, $n={len(x[1])}$')
    plt.hist(x[2], **kwargs, color='r', label=f'${seperator}$={seperator_list[2]}, $n={len(x[2])}$')
    plt.xlim(0,2)
    plt.axvline(x=1.3,label="$A_{c}$",color='black')
    plt.xlabel("Apical Area")
    plt.ylabel("Density")
    plt.text(0.85,0.5,s=f"D = {df_alive_final_timepoint['D'].iloc[0]}\n"+f"{constant}={df_alive_final_timepoint[constant].iloc[0]}"+f"\n Bins={bins_list[0]}",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
    plt.legend(loc="upper right",prop={'size': 8})
    plt.show()


    #now need histogram of age
    df_dead_final_timepoint=df_master.loc[(df_master.dead==True)&(df_master.time==306)]
    #now collect ages seperated by seperator
    seperator_list=list(df_dead_final_timepoint[seperator].cat.categories)
    x=[[],[],[]]
    for i in range(0,len(seperator_list)):
        x[i]=(df_dead_final_timepoint.loc[(df_dead_final_timepoint[seperator]==seperator_list[i])]['age'])
    
    kwargs = dict(bins=bins_list[1],histtype='step',density=True)
    fig2, ax=plt.subplots(dpi=500)
    plt.subplots_adjust(left=0.2,bottom=0.2, top = 0.9, right = 0.9)
    plt.hist(x[0], **kwargs, color='g', label=f'${seperator}$={seperator_list[0]}, $n={len(x[0])}$')
    plt.hist(x[1], **kwargs, color='b', label=f'${seperator}$={seperator_list[1]}, $n={len(x[1])}$')
    plt.hist(x[2], **kwargs, color='r', label=f'${seperator}$={seperator_list[2]}, $n={len(x[2])}$')
    plt.xlim(0,3)
    plt.axvline(x=1-t_M,label="$t_{c}$",color="black")
    plt.xlabel("Cell Cycle Length")
    plt.ylabel("Density")
    plt.text(0.85,0.5,s=f"D = {df_dead_final_timepoint['D'].iloc[0]}\n"+f"{constant}={df_dead_final_timepoint[constant].iloc[0]}"+f"\n Bins={bins_list[1]}",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
    plt.legend(loc="upper right")
    plt.show()



    #now need histogram of nucl_pos
    #need alive cells at final timepoint
    df_alive_final_timepoint=df_master.loc[(df_master.dead==False)&(df_master.time==306)]
    
    seperator_list=list(df_dead_final_timepoint[seperator].cat.categories)
    x=[[],[],[]]
    for i in range(0,len(seperator_list)):
        x[i]=(df_alive_final_timepoint.loc[(df_alive_final_timepoint[seperator]==seperator_list[i])]['nucl_pos'])

    kwargs = dict(bins=bins_list[2],histtype='step',density=True)
    fig3, ax=plt.subplots(dpi=500)
    plt.subplots_adjust(left=0.2,bottom=0.2, top = 0.9, right = 0.9)
    plt.hist(x[0], **kwargs, color='g', label=f'${seperator}$={seperator_list[0]}, $n={len(x[0])}$')
    plt.hist(x[1], **kwargs, color='b', label=f'${seperator}$={seperator_list[1]}, $n={len(x[1])}$')
    plt.hist(x[2], **kwargs, color='r', label=f'${seperator}$={seperator_list[2]}, $n={len(x[2])}$')
    plt.xlabel("Nuclear Position")
    plt.axvline(x=0.75,color='black',label="$z_{c}$")
    plt.ylabel("Density")
    plt.text(0.41,0.8,s=f"D = {df_alive_final_timepoint['D'].iloc[0]}\n"+f"{constant}={df_alive_final_timepoint[constant].iloc[0]}"+f"\n Bins={bins_list[2]}",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
    plt.legend(loc="upper left")
    plt.show()

    #experimental plot of z against age
   
    df_alive=df_master.loc[(df_master.dead==False)]
    nuclear_positions=[[],[],[]]
    ages=[[],[],[]]
    for i in range(0,len(seperator_list)):
        nuclear_positions[i]=df_alive.loc[(df_alive[seperator]==seperator_list[i])]['nucl_pos']
        ages[i]=df_alive.loc[(df_alive[seperator]==seperator_list[i])]['age']
    fig4,ax=plt.subplots(dpi=500,)
    plt.hist2d(ages[0],nuclear_positions[0],bins=bins_list[3],norm=mpl.colors.LogNorm(),cmap='inferno')
    cb=plt.colorbar(label="$Log_{10}($Count in bin$)$")
    plt.xlabel("Age (AU)")
    plt.ylabel("Nuclear_position (AU)")
    plt.axhline(y=0.75,color='red',label=r"$z_{c}$")
    plt.axvline(x=t_G1,color="black")
    plt.axvline(x=t_G1+t_S,color="black")
    plt.axvline(x=t_G1+t_S+t_G2,color="black")
    plt.axvline(x=t_G1+t_S+t_G2+t_M,color="black")
    title = ax.set_title("\n".join(wrap(f"{seperator}={seperator_list[0]}, "+f"{constant}={df_alive[constant].iloc[0]}, "+f"$D$={df_alive['D'].iloc[0]}, "+f"Bins={bins_list[3]}", 60)))
    title.set_y(1.05)
    
    
    plt.text(t_G1/2, 0.5,s='G1',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5))
    plt.text(np.mean([t_G1,t_G1+t_S]), 0.5,s='S',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5))
    plt.text(np.mean([t_G1+t_S,t_G1+t_S+t_G2]), 0.5,s='G2',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5))
    plt.text(np.mean([t_G1+t_S+t_G2,t_G1+t_S+t_G2+t_M]), 0.5,s='M',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5))

    fig5,ax=plt.subplots(dpi=500)
    plt.hist2d(ages[1],nuclear_positions[1],bins=bins_list[4],norm=mpl.colors.LogNorm(),cmap='inferno')
    cb=plt.colorbar(label="$Log_{10}($Count in bin$)$")
    plt.xlabel("Age (AU)")
    plt.ylabel("Nuclear_position (AU)")
    plt.axhline(y=0.75,color='red',label=r"$z_{c}$")
    plt.axvline(x=t_G1,color="black")
    plt.axvline(x=t_G1+t_S,color="black")
    plt.axvline(x=t_G1+t_S+t_G2,color="black")
    plt.axvline(x=t_G1+t_S+t_G2+t_M,color="black")
    title = ax.set_title("\n".join(wrap(f"{seperator}={seperator_list[1]}, "+f"{constant}={df_alive[constant].iloc[0]}, "+f"$D$={df_alive['D'].iloc[0]}, "+f"Bins={bins_list[4]}", 60)))
    title.set_y(1.05)
    
    plt.text(t_G1/2, 0.5,s='G1',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5))
    plt.text(np.mean([t_G1,t_G1+t_S]), 0.5,s='S',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5))
    plt.text(np.mean([t_G1+t_S,t_G1+t_S+t_G2]), 0.5,s='G2',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5))
    plt.text(np.mean([t_G1+t_S+t_G2,t_G1+t_S+t_G2+t_M]), 0.5,s='M',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5))


    fig6,ax=plt.subplots(dpi=500)
    plt.hist2d(ages[2],nuclear_positions[2],bins=bins_list[5],norm=mpl.colors.LogNorm(),cmap='inferno')
    cb=plt.colorbar(label="$Log_{10}($Count in bin$)$")
    plt.xlabel("Age (AU)")
    plt.ylabel("Nuclear_position (AU)")
    plt.axhline(y=0.75,color='red',label=r"$z_{c}$")
    plt.axvline(x=t_G1,color="black")
    plt.axvline(x=t_G1+t_S,color="black")
    plt.axvline(x=t_G1+t_S+t_G2,color="black")
    plt.axvline(x=t_G1+t_S+t_G2+t_M,color="black")
    title = ax.set_title("\n".join(wrap(f"{seperator}={seperator_list[2]}, "+f"{constant}="+f"{df_alive[constant].iloc[0]}, "+f"$D$={df_alive['D'].iloc[0]}, "+f"Bins={bins_list[5]}", 60)))
    title.set_y(1.05)
    
    
    plt.text(t_G1/2, 0.5,s='G1',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5))
    plt.text(np.mean([t_G1,t_G1+t_S]), 0.5,s='S',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5))
    plt.text(np.mean([t_G1+t_S,t_G1+t_S+t_G2]), 0.5,s='G2',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5))
    plt.text(np.mean([t_G1+t_S+t_G2,t_G1+t_S+t_G2+t_M]), 0.5,s='M',
     horizontalalignment='center',
     verticalalignment='center',
     bbox=dict(facecolor='none', alpha=0.5))

  
    #we need alive cells at final timepoint
    df_alive_final_timepoint['neighbour_ids']=[i[1:-1].split(",") for i in df_alive_final_timepoint['neighbour_ids']]
    
    a=[[],[],[]]
    b=[[],[],[]]
    n=[[],[],[]]
    for iii in range(0,len(seperator_list)):#iii is the seperator
        df=df_alive_final_timepoint.loc[df_alive_final_timepoint[seperator]==seperator_list[iii]]
        n[iii]=len(df)
        for i in range(0,len(df['neighbour_ids'])): # i is each cell in each seperated simulation
            neighbours=df['neighbour_ids'].iloc[i]
            neighbours=[int(j) for j in neighbours]
            for ii in range(0,len(neighbours)): #ii is the index of the i'th cell's neighbour
                if df.loc[df.cell_index==neighbours[ii]].empty==False:
                    a[iii].append(df['apical_area'].iloc[i])
                    b[iii].append(df.loc[df.cell_index==neighbours[ii]]['apical_area'].values[0])
    
    
    fig7, ax=plt.subplots(dpi=500)
    
    plt.hexbin(a[0], b[0],gridsize=bins_list[6],  norm=mpl.colors.LogNorm(),linewidths=0.2,cmap='inferno')

    z = np.linspace(0,max(a[0]))
    plt.plot(z,z, marker=",", alpha=1,color='r',)
    cb = plt.colorbar(label='$Log_{10}($count in bin$)$')

    plt.xlabel("Cell Apical Area")
    plt.ylabel("Adjacent cell Apical Area")
    title = ax.set_title("\n".join(wrap(f"{seperator}={seperator_list[0]}, "+f"{constant}="+f"{df_alive_final_timepoint[constant].iloc[0]}, "+f"$D$={df_alive_final_timepoint['D'].iloc[0]}, "+f"Bins={bins_list[6]}"+f", $n$={n[0]}", 60)))
    plt.axvline(1.3,color='blue')
    plt.axhline(1.3,color='blue')
    title.set_y(1.05)


    
    fig8, ax=plt.subplots(dpi=500)
    
    plt.hexbin(a[1], b[1],gridsize=bins_list[7],  norm=mpl.colors.LogNorm(),linewidths=0.2,cmap='inferno')

    z = np.linspace(0,max(a[1]))
    plt.plot(z,z, marker=",", alpha=1,color='r',)
    cb = plt.colorbar(label='$Log_{10}($count in bin$)$')

    plt.xlabel("Cell apical area")
    plt.ylabel("Adjacent cell apical area")
    title = ax.set_title("\n".join(wrap(f"{seperator}={seperator_list[1]}, "+f"{constant}="+f"{df_alive_final_timepoint[constant].iloc[0]}, "+f"$D$={df_alive_final_timepoint['D'].iloc[0]}, "+f"Bins={bins_list[7]}"+f", $n$={n[1]}", 60)))
    title.set_y(1.05)
    
    plt.axvline(1.3,color='blue')
    plt.axhline(1.3,color='blue')
    
    fig9, ax=plt.subplots(dpi=500)
    
    plt.hexbin(a[2], b[2],gridsize=bins_list[8],  norm=mpl.colors.LogNorm(),linewidths=0.2,cmap='inferno')

    z = np.linspace(0,max(a[2]))
    plt.plot(z,z, marker=",", alpha=1,color='r',)
    cb = plt.colorbar(label='$Log_{10}($count in bin$)$')

    plt.xlabel("Cell apical area")
    plt.ylabel("Adjacent cell apical area")
    title = ax.set_title("\n".join(wrap(f"{seperator}={seperator_list[2]}, "+f"{constant}="+f"{df_alive_final_timepoint[constant].iloc[0]}, "+f"$D$={df_alive_final_timepoint['D'].iloc[0]}, "+f"Bins={bins_list[8]}"+f", $n$={n[2]}", 60)))
    title.set_y(1.05)
    
    plt.axvline(1.3,color='blue')
    plt.axhline(1.3,color='blue')


    #now final plot is n_neighbours distribution by seperator
    #we want
    #df_alive_final_timepoint
    n_neighbours=[[],[],[]]
    for i in range(0,len(seperator_list)):
        n_neighbours[i]=df_alive_final_timepoint.loc[(df_alive_final_timepoint[seperator]==seperator_list[i])]['n_neighbours']
    df_alive_final_timepoint['n_neighbours']=pd.to_numeric(df_alive_final_timepoint['n_neighbours'])
    df_alive_final_timepoint['a']=pd.to_numeric(df_alive_final_timepoint['a'])

    seperator_vector=[]
    for i in range(0,len(n_neighbours)):
        seperator_vector.append([seperator_list[i]]*len(n_neighbours[i])) 
    #now going to make dataframe of n_neighbours and the associated seperator
    #unnesting these two lists
    n_neighbours=list(chain.from_iterable(n_neighbours))
    seperator_vector=list(chain.from_iterable(seperator_vector))
    df=pd.DataFrame({'n_neighbours': n_neighbours,f'{seperator}':seperator_vector})
 
        #The data
    a=df.query(f'{seperator}=={seperator_list[0]}')['n_neighbours']
    b=df.query(f'{seperator}=={seperator_list[1]}')['n_neighbours']
    c=df.query(f'{seperator}=={seperator_list[2]}')['n_neighbours']

    categories=list(set(df['n_neighbours']))
    
    x=[]
    y=[]
    z=[]
    for i in range(0,len(categories)):
        x.append(Counter(a)[categories[i]])
        y.append(Counter(b)[categories[i]])
        z.append(Counter(c)[categories[i]])
    x=np.array(x)
    y=np.array(y)
    z=np.array(z)
    counts=list(x)+list(y)+list(z)
    proportion=list(x/sum(x))+list(y/sum(y))+list(z/sum(z))

    seperator_legend_label=list(chain.from_iterable([[f'{seperator_list[0]}, n={len(a)}']*len(categories),[f'{seperator_list[1]}, n={len(b)}']*len(categories),[f'{seperator_list[2]}, n={len(c)}']*len(categories)]))

    df_counts=pd.DataFrame({f'{seperator}':seperator_legend_label,'n_neighbours':categories*3,'counts':counts,'proportion':proportion})

    
    fig10=plt.figure()
    p=sns.barplot(data=df_counts,x='n_neighbours',y='proportion',hue=f'{seperator}',palette=sns.color_palette('Blues',3))
    p.set_title(("\n".join(wrap(f"{constant} = "+f"{df_alive_final_timepoint[constant].iloc[0]}, "+f"$D$ = {df_alive_final_timepoint['D'].iloc[0]}, ", 60))))
    fig10=p.figure

    #now gonna be writing a trajectory plot function
    #first gonna take out useless collumn from df_master
    #I need time, cell_index, age, k_g1, D, a, neighbour_ids, n_neighbours, apical_area ,dead , nucl_pos

    df_trajectory=df_master[["time", "cell_index" ,"dead", "n_neighbours","neighbour_ids","k_g1", "D", "a", "age", "apical_area"  , "nucl_pos"]].copy()
    
    return[fig1,fig2,fig3,fig4,fig5,fig6,fig7,fig8,fig9,fig10]
#now run the analyses
#excluding constant =0
#%% the plots
#before
#experiment_names=index_of_simulation(model_type="active",k_g1=0,d=1e-5,duration=306)[1][3:6]
#after
metadataframe=make_simulation_metadataframe()
#now we are gonna conduct analysis on type 9 - new age equation
experiment_names=list(metadataframe.query('k_g1==0 and k_m==29 and model_type=="active_10"and D==1e-5 and duration==306.0 and 0.1<=a<=10').iloc[:,0])

seperator='a'
bins_list=[50,100,100,500,500,500,200,200,200]
constant='k_g1'
passive_model=plot_3_histograms(experiment_names=experiment_names,seperator=seperator,bins_list=bins_list,constant=constant)

for i in range(0,len(passive_model)):
    passive_model[i].savefig(os.path.join("output","passive_model",f"fig{i}.jpg"))


experiment_names=list(metadataframe.query('a==0 and D==1e-5 and duration==306 and model_type=="active_10"and (k_g1==1 or k_g1==6 or k_g1==12)').iloc[:,0])
seperator='k_g1'
bins_list=[50,150,80,500,500,500,200,200,200]
constant='a'
active_model=plot_3_histograms(experiment_names=experiment_names,seperator=seperator,bins_list=bins_list,constant=constant)


for i in range(0,len(passive_model)):
    active_model[i].savefig(os.path.join("output","active_model",f"fig{i}.jpg"))

#now mixed model
experiment_names=list(metadataframe.query('a==1 and D==1e-5 and model_type=="active_10" and duration==306 and (k_g1==1 or k_g1==6 or k_g1==12)').iloc[:,0])

seperator='k_g1'
bins_list=[50,150,80,500,500,500,200,200,200]
constant='a'
mixed_model_constant_a=plot_3_histograms(experiment_names=experiment_names,seperator=seperator,bins_list=bins_list,constant=constant)

for i in range(0,len(mixed_model_constant_a)):
    mixed_model_constant_a[i].savefig(os.path.join("output","mixed_model_constant_a",f"fig{i}.jpg"))


experiment_names=list(metadataframe.query('k_g1==6 and duration ==306 and model_type=="active_10" and  D==1e-5 and a!=0').iloc[:,0])
seperator='a'
bins_list=[50,150,80,500,500,500,200,200,200]
constant='k_g1'
mixed_model_constant_k_g1=plot_3_histograms(experiment_names=experiment_names,seperator=seperator,bins_list=bins_list,constant=constant)

for i in range(0,len(mixed_model_constant_k_g1)):
    mixed_model_constant_k_g1[i].savefig(os.path.join("output","mixed_model_constant_k_g1",f"fig{i}.jpg"))


#now plotting new experiments with a<0.1
#so far we have explored a=[1e1,1e0,1e-1]
#now gonna do a=[1e-2,1e-3,1e-4]
experiment_names=list(metadataframe.query('k_g1==0 and k_m==29 and duration ==306 and model_type=="active_10" and D==1e-5 and 1e-5<a<=1e-2').iloc[:,0])
seperator='a'
bins_list=[50,150,80,500,500,500,200,200,200]
constant='k_g1'
passive_model_micro_crowding=plot_3_histograms(experiment_names=experiment_names,seperator=seperator,bins_list=bins_list,constant=constant)
for i in range(0,len(passive_model_micro_crowding)):
    passive_model_micro_crowding[i].savefig(os.path.join("output","passive_model_micro_crowding",f"fig{i}.jpg"))
# %%
