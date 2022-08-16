#%%
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

#parse all _simulations using make_df_simulations.sh
#figure 8 panel b & c

#%%
#d=0
metadata_indices=[0,0]
metadata_indices[0]=index_of_simulation(model_type=1,d=0,k=100)[0]+index_of_simulation(model_type=1,d=0,k=50)[0]+index_of_simulation(model_type=1,d=0,k=5)[0]
metadata_indices[1]=index_of_simulation(model_type=1,d=0,k=100)[1]+index_of_simulation(model_type=1,d=0,k=50)[1]+index_of_simulation(model_type=1,d=0,k=5)[1]

#will only read certain collumns to save memory
collumns=['file_index','file_name',"time_point","cell_index","age","dead","zposn","A0","k","D","nucl_pos","apical_area"]
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead', 'k','D']:
            df[col]=df[col].astype('category')
    df_master=pd.concat([df_master,df])
    del df
df=df_master

mpl.rcParams.update(mpl.rcParamsDefault)
#filter df_master for cells that are dead
#%%
k=[100,50,5]
df_master=df
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306") &(df_master.dead==False)  ]
xlabel="A0"
x1=df_master.loc[df_master.k==float(k[0]),xlabel]
x2=df_master.loc[df_master.k==float(k[1]),xlabel]
x3=df_master.loc[df_master.k==float(k[2]),xlabel]
kwargs = dict(bins=200,histtype='step',density=True)


fig, ax=plt.subplots()
plt.hist(x1, **kwargs, color='g', label=f'K={k[0]}')
plt.hist(x2, **kwargs, color='b', label=f'K={k[1]}')
plt.hist(x3, **kwargs, color='r', label=f'K={k[2]}')
#plt.xlim(-0.1,3)
plt.xlabel(xlabel.capitalize())
plt.ylabel("Density")
plt.text(0.89,0.7,s=f"D={df_master['D'].iloc[0]}",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.legend()
plt.show()



#%%
df_master=df
df_master=df_master.loc[df_master.dead==True]

xlabel="age"
x1=df_master.loc[df_master.k==float(k[0]),xlabel]
x2=df_master.loc[df_master.k==float(k[1]),xlabel]
x3=df_master.loc[df_master.k==float(k[2]),xlabel]
kwargs = dict(bins=1000,histtype='step',density=True)

fig1, ax=plt.subplots()
plt.hist(x1, **kwargs, color='g', label=f'K={k[0]}')
plt.hist(x2, **kwargs, color='b', label=f'K={k[1]}')
plt.hist(x3, **kwargs, color='r', label=f'K={k[2]}')
#plt.xlim(0.8,1.5)
plt.xlabel(xlabel.capitalize())
plt.ylabel("Density")
plt.text(0.9,0.7,s=f"D={df_master['D'].iloc[0]}",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.legend();
plt.show()

#%%
df_master=df
df_master=df_master.loc[(df_master.dead==False) &(df_master.time_point=="timepoint_306_of_306")]
xlabel="apical_area"
x1=df_master.loc[df_master.k==float(k[0]),xlabel]
x2=df_master.loc[df_master.k==float(k[1]),xlabel]
x3=df_master.loc[df_master.k==float(k[2]),xlabel]
kwargs = dict(bins=80,histtype='step',density=True)


fig2, ax=plt.subplots()
plt.hist(x1, **kwargs, color='g', label=f'K={k[0]}')
plt.hist(x2, **kwargs, color='b', label=f'K={k[1]}')
plt.hist(x3, **kwargs, color='r', label=f'K={k[2]}')
# plt.xlim(-0.1,1)
plt.xlabel(xlabel.capitalize())
plt.ylabel("Density")
plt.text(0.9,0.7,s=f"D={df_master['D'].iloc[0]}",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.legend();
plt.show()

# %%

pp = PdfPages(os.path.join("output","Figure_8_panel_b_D=0.pdf"))
pp.savefig(fig)
pp.savefig(fig1)
pp.savefig(fig2)
pp.close()

#now d=0.01
#%%
#getting the right files
metadata_indices=[0,0]
d=0.01
metadata_indices[0]=index_of_simulation(model_type=1,d=d,k=100)[0]+index_of_simulation(model_type=1,d=d,k=50)[0]+index_of_simulation(model_type=1,d=d,k=5)[0]
metadata_indices[1]=index_of_simulation(model_type=1,d=d,k=100)[1]+index_of_simulation(model_type=1,d=d,k=50)[1]+index_of_simulation(model_type=1,d=d,k=5)[1]

#will only read certain collumns to save memory
collumns=['file_index','file_name',"time_point","cell_index","age","dead","zposn","A0","k","D","nucl_pos","apical_area"]
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    
    df_master=pd.concat([df_master,df])
    del df
for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead','D']:
    df_master[col]=df_master[col].astype('category')    
df=df_master

mpl.rcParams.update(mpl.rcParamsDefault)
#filter df_master for cells that are dead
#%%
df_master=df
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306") &(df_master.dead==False)  ]
xlabel="A0"
x1=df_master.loc[df_master.k==float(k[0]),xlabel]
x2=df_master.loc[df_master.k==float(k[1]),xlabel]
x3=df_master.loc[df_master.k==float(k[2]),xlabel]
kwargs = dict(bins=200,histtype='step',density=True)


fig, ax=plt.subplots()
plt.hist(x1, **kwargs, color='g', label=f'K={k[0]}')
plt.hist(x2, **kwargs, color='b', label=f'K={k[1]}')
plt.hist(x3, **kwargs, color='r', label=f'K={k[2]}')
#plt.xlim(-0.1,3)
plt.xlabel(xlabel.capitalize())
plt.ylabel("Density")
plt.text(0.89,0.7,s=f"D={df_master['D'].iloc[0]}",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.legend()
plt.show()



#%%
df_master=df
df_master=df_master.loc[df_master.dead==True]

xlabel="age"
x1=df_master.loc[df_master.k==float(k[0]),xlabel]
x2=df_master.loc[df_master.k==float(k[1]),xlabel]
x3=df_master.loc[df_master.k==float(k[2]),xlabel]
kwargs = dict(bins=1000,histtype='step',density=True)

fig1, ax=plt.subplots()
plt.hist(x1, **kwargs, color='g', label=f'K={k[0]}')
plt.hist(x2, **kwargs, color='b', label=f'K={k[1]}')
plt.hist(x3, **kwargs, color='r', label=f'K={k[2]}')
#plt.xlim(0.8,1.5)
plt.xlabel(xlabel.capitalize())
plt.ylabel("Density")
plt.text(0.9,0.7,s=f"D={df_master['D'].iloc[0]}",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.legend();
plt.show()

#%%
df_master=df
df_master=df_master.loc[(df_master.dead==False) &(df_master.time_point=="timepoint_306_of_306")]
xlabel="apical_area"
x1=df_master.loc[df_master.k==float(k[0]),xlabel]
x2=df_master.loc[df_master.k==float(k[1]),xlabel]
x3=df_master.loc[df_master.k==float(k[2]),xlabel]
kwargs = dict(bins=80,histtype='step',density=True)


fig2, ax=plt.subplots()
plt.hist(x1, **kwargs, color='g', label=f'K={k[0]}')
plt.hist(x2, **kwargs, color='b', label=f'K={k[1]}')
plt.hist(x3, **kwargs, color='r', label=f'K={k[2]}')
#plt.xlim(-0.1,1)
plt.xlabel(xlabel.capitalize())
plt.ylabel("Density")
plt.text(0.9,0.7,s=f"D={df_master['D'].iloc[0]}",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.legend();
plt.show()


# %%

pp = PdfPages(os.path.join("output","Figure_8_panel_b_D=0.01.pdf"))
pp.savefig(fig)
pp.savefig(fig1)
pp.savefig(fig2)
pp.close()

#%%
#figure 9 panel b
#cell cycle length histogram categorised by D

k=100
d=[0,0.01,0.05,0.15,0.75]
metadata_indices[0]=index_of_simulation(model_type=1,d=d[0],k=k)[0]+index_of_simulation(model_type=1,d=d[1],k=k)[0]+index_of_simulation(model_type=1,d=d[2],k=k)[0]+index_of_simulation(model_type=1,d=d[3],k=k)[0]+index_of_simulation(model_type=1,d=d[4],k=k)[0]
metadata_indices[1]=index_of_simulation(model_type=1,d=d[0],k=k)[1]+index_of_simulation(model_type=1,d=d[1],k=k)[1]+index_of_simulation(model_type=1,d=d[2],k=k)[1]+index_of_simulation(model_type=1,d=d[3],k=k)[1]+index_of_simulation(model_type=1,d=d[4],k=k)[1]

collumns=['file_index','file_name',"time_point","cell_index","age","dead","zposn","A0","k","D","nucl_pos","apical_area"]
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    df_master=pd.concat([df_master,df])
    del df
for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead', 'k','D']:
            df_master[col]=df_master[col].astype('category')
df=df_master
#%%
df_master=df_master.loc[df_master.dead==True]

xlabel="age"
x1=df_master.loc[df_master.D==float(d[0]),xlabel]
x2=df_master.loc[df_master.D==float(d[1]),xlabel]
x3=df_master.loc[df_master.D==float(d[2]),xlabel]
x4=df_master.loc[df_master.D==float(d[3]),xlabel]
x5=df_master.loc[df_master.D==float(d[4]),xlabel]

kwargs = dict(bins=750,histtype='step',density=True)

fig, ax=plt.subplots()
plt.hist(x1, **kwargs, color='b', label=f'D={d[0]}')
plt.hist(x2, **kwargs, color='orange', label=f'D={d[1]}')
plt.hist(x3, **kwargs, color='green', label=f'D={d[2]}')
plt.hist(x4,**kwargs,color='red',label=f"D={d[3]}")
plt.hist(x5,**kwargs,color='purple',label=f"D={d[4]}")
#plt.xlim(0.8,1.5)
plt.xlabel(xlabel.capitalize())
plt.ylabel("Density")
plt.text(0.85,0.55,s=f"K={df_master['k'].iloc[0]}\nr=0",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.legend();
plt.show()

#%%
df_master=df
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306")  &(df_master.dead==False) ]
xlabel="A0"
x1=df_master.loc[df_master.D==float(d[0]),xlabel]
x2=df_master.loc[df_master.D==float(d[1]),xlabel]
x3=df_master.loc[df_master.D==float(d[2]),xlabel]
x4=df_master.loc[df_master.D==float(d[3]),xlabel]
x5=df_master.loc[df_master.D==float(d[4]),xlabel]
kwargs = dict(bins=30,histtype='step',density=True)


fig1, ax=plt.subplots()
plt.hist(x1, **kwargs, color='b', label=f'D={d[0]}')
plt.hist(x2, **kwargs, color='orange', label=f'D={d[1]}')
plt.hist(x3, **kwargs, color='green', label=f'D={d[2]}')
plt.hist(x4,**kwargs,color='red',label=f"D={d[3]}")
plt.hist(x5,**kwargs,color='purple',label=f"D={d[4]}")
#plt.xlim(-0.1,3)
plt.xlabel(xlabel.capitalize())
plt.ylabel("Density")
plt.text(0.85,0.55,s=f"K={df_master['k'].iloc[0]}\nr=0",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.legend()
plt.show()


#%%
df_master=df
df_master=df_master.loc[(df_master.dead==False) &(df_master.time_point=="timepoint_306_of_306")]
xlabel="apical_area"
x1=df_master.loc[df_master.D==float(d[0]),xlabel]
x2=df_master.loc[df_master.D==float(d[1]),xlabel]
x3=df_master.loc[df_master.D==float(d[2]),xlabel]
x4=df_master.loc[df_master.D==float(d[3]),xlabel]
x5=df_master.loc[df_master.D==float(d[4]),xlabel]
kwargs = dict(bins=50,histtype='step',density=True)


fig2, ax=plt.subplots()
plt.hist(x1, **kwargs, color='b', label=f'D={d[0]}')
plt.hist(x2, **kwargs, color='orange', label=f'D={d[1]}')
plt.hist(x3, **kwargs, color='green', label=f'D={d[2]}')
plt.hist(x4,**kwargs,color='red',label=f"D={d[3]}")
plt.hist(x5,**kwargs,color='purple',label=f"D={d[4]}")
# plt.xlim(-0.1,1)
plt.xlabel(xlabel.capitalize())
plt.ylabel("Density")
plt.text(0.85,0.55,s=f"K={df_master['k'].iloc[0]}\nr=0",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.legend();
plt.show()

#%%
pp = PdfPages(os.path.join("output","Figure_9.pdf"))
pp.savefig(fig)
pp.savefig(fig1)
pp.savefig(fig2)
pp.close()

#%%
#figure 10
#same three histograms using r as the levels

d=0
k=100
r=[0,0.1,0.3,1,3]
metadata_indices=index_of_simulation(model_type=2,d=d,k=k,r=r[0])+index_of_simulation(model_type=2,d=d,k=k,r=r[1])+index_of_simulation(model_type=2,d=d,k=k,r=r[2])+index_of_simulation(model_type=2,d=d,k=k,r=r[3])+index_of_simulation(model_type=2,d=d,k=k,r=r[4])
a=metadata_indices
b=[[],[]]
for i in range(0,len(a),2):#even represents index
   
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
metadata_indices=b

collumns=['file_index','file_name',"time_point","cell_index","age","dead","zposn","A0","k","D","nucl_pos","apical_area","r"]
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    
    df_master=pd.concat([df_master,df])
    del df
for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead', 'k','D']:
            df_master[col]=df_master[col].astype('category')    
df=df_master
#%%
df_master=df
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306")  &(df_master.dead==False) ]
xlabel="A0"
x1=df_master.loc[df_master.r==float(r[0]),xlabel]
x2=df_master.loc[df_master.r==float(r[1]),xlabel]
x3=df_master.loc[df_master.r==float(r[2]),xlabel]
x4=df_master.loc[df_master.r==float(r[3]),xlabel]
x5=df_master.loc[df_master.r==float(r[4]),xlabel]
kwargs = dict(bins=30,histtype='step',density=True)


fig, ax=plt.subplots()
plt.hist(x1, **kwargs, color='b', label=f'r={r[0]}')
plt.hist(x2, **kwargs, color='orange', label=f'r={r[1]}')
plt.hist(x3, **kwargs, color='green', label=f'r={r[2]}')
plt.hist(x4,**kwargs,color='red',label=f"r={r[3]}")
plt.hist(x5,**kwargs,color='purple',label=f"r={r[4]}")
#plt.xlim(-0.1,3)
plt.xlabel(xlabel.capitalize())
plt.ylabel("Density")
plt.text(0.85,0.55,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.legend()
plt.show()


#%%
df_master=df
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306")  &(df_master.dead==True)]

xlabel="age"
x1=df_master.loc[df_master.r==float(r[0]),xlabel]
x2=df_master.loc[df_master.r==float(r[1]),xlabel]
x3=df_master.loc[df_master.r==float(r[2]),xlabel]
x4=df_master.loc[df_master.r==float(r[3]),xlabel]
x5=df_master.loc[df_master.r==float(r[4]),xlabel]

kwargs = dict(bins=170,histtype='step',density=True)

fig1, ax=plt.subplots()
plt.hist(x1, **kwargs, color='b', label=f'r={r[0]}')
plt.hist(x2, **kwargs, color='orange', label=f'r={r[1]}')
plt.hist(x3, **kwargs, color='green', label=f'r={r[2]}')
plt.hist(x4,**kwargs,color='red',label=f"r={r[3]}")
plt.hist(x5,**kwargs,color='purple',label=f"r={r[4]}")
# plt.xlim(0.8,1.5)
plt.xlabel(xlabel.capitalize())
plt.ylabel("Density")
plt.text(0.85,0.55,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.legend();
plt.show()

#%%
df_master=df
df_master=df_master.loc[(df_master.dead==False) &(df_master.time_point=="timepoint_306_of_306")]
xlabel="apical_area"
x1=df_master.loc[df_master.r==float(r[0]),xlabel]
x2=df_master.loc[df_master.r==float(r[1]),xlabel]
x3=df_master.loc[df_master.r==float(r[2]),xlabel]
x4=df_master.loc[df_master.r==float(r[3]),xlabel]
x5=df_master.loc[df_master.r==float(r[4]),xlabel]
kwargs = dict(bins=100,histtype='step',density=True)


fig2, ax=plt.subplots()
plt.hist(x1, **kwargs, color='b', label=f'r={r[0]}')
plt.hist(x2, **kwargs, color='orange', label=f'r={r[1]}')
plt.hist(x3, **kwargs, color='green', label=f'r={r[2]}')
plt.hist(x4,**kwargs,color='red',label=f"r={r[3]}")
plt.hist(x5,**kwargs,color='purple',label=f"r={r[4]}")
#plt.xlim(-0.1,1)
plt.xlabel(xlabel.capitalize())
plt.ylabel("Density")
plt.text(0.85,0.55,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.legend();
plt.show()

#%%
pp = PdfPages(os.path.join("output","Figure_10.pdf"))
pp.savefig(fig)
pp.savefig(fig1)
pp.savefig(fig2)
pp.close()


#figure 11 2d histograms
#%%
#k=100,d=0.0,r=0.0
D=0
k=100
r=0
metadata_indices=index_of_simulation(model_type=2,d=D,k=k,r=r)
a=metadata_indices
b=[[],[]]
for i in range(0,len(a),2):#even represents index
   
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
metadata_indices=b

collumns=['file_index','file_name',"time_point","cell_index","age","dead","zposn","A0","k","D","nucl_pos","apical_area","r","neighbour_ids"]
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    
    df['neighbour_ids']=df['neighbour_ids'].apply(literal_eval)
    df_master=pd.concat([df_master,df])
    del df
for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead', 'k','D','r']:
            df_master[col]=df_master[col].astype('category')
df=df_master

#%%
#want final timepoint
#NB neigbours method already filters for only alive cells so we dont need to do this ourselves
df=df_master
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306")  ]
df_master=df_master.reset_index()

#pairwise cell-neighbour indexes
x=[] #list of cell indexes having neigbours
y=[] #list of neighbour indexes for cell x[i] - 
for i in range(0,len(df_master)):
    a=df_master.cell_index[i]
    b=df_master.neighbour_ids[i]
    if b==[]:
        pass
    else:
        for j in range(0,len(b)):
            x.append(a)
            y.append(b[j])
xlabel='nucl_pos'
x=np.asarray(x)
y=np.asarray(y)
#now nuclear position 2d histogram
a=[]
b=[]
for i in range(0,len(x)):
    a.append(df_master[xlabel][x[i]])
    b.append(df_master[xlabel][y[i]])

pearson=scipy.stats.pearsonr(a,b)

#%%plot
#plt.hexbin(a, b, gridsize=30, cmap='Blues')
fig, ax=plt.subplots()
plt.hexbin(a, b,gridsize=30, norm=mpl.colors.LogNorm(), cmap="Blues")
z = np.linspace(0,1)
plt.plot(z,z, marker=",", alpha=1,color='r',)
cb = plt.colorbar(label='count in bin')

plt.xlabel(xlabel.capitalize())
plt.ylabel("Neigbour "+str(xlabel))
plt.text(0.1,0.6,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}\nr={df_master['r'].iloc[0]}\nPearson's r={pearson[0]:.2f}\np={pearson[1]:.2f} ",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)

plt.show()

#%%
#figure 11
#K=100,d=0.75,r=0
D=0.75
k=100
r=0
metadata_indices=index_of_simulation(model_type=2,d=D,k=k,r=r)
a=metadata_indices
b=[[],[]]
for i in range(0,len(a),2):#even represents index
   
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
metadata_indices=b

collumns=['file_index','file_name',"time_point","cell_index","age","dead","zposn","A0","k","D","nucl_pos","apical_area","r","neighbour_ids"]
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    df['neighbour_ids']=df['neighbour_ids'].apply(literal_eval)
    df_master=pd.concat([df_master,df])
    del df
for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead', 'k','D','r']:
            df_master[col]=df_master[col].astype('category')
df=df_master

#%%
#want final timepoint
df=df_master
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306")]
df_master=df_master.reset_index()

#pairwise cell-neighbour indexes
x=[] #list of cell indexes having neigbours
y=[] #list of neighbour indexes for cell x[i] - 
for i in range(0,len(df_master)):
    a=df_master.cell_index[i]
    b=df_master.neighbour_ids[i]
    if b==[]:
        pass
    else:
        for j in range(0,len(b)):
            x.append(a)
            y.append(b[j])
xlabel='nucl_pos'
x=np.asarray(x)
y=np.asarray(y)
#now nuclear position 2d histogram
a=[]
b=[]
for i in range(0,len(x)):
    a.append(df_master[xlabel][x[i]])
    b.append(df_master[xlabel][y[i]])

pearson=scipy.stats.pearsonr(a,b)

#%%plot
#plt.hexbin(a, b, gridsize=30, cmap='Blues')
fig1, ax=plt.subplots()
plt.hexbin(a, b,gridsize=30,  cmap="Blues")
z = np.linspace(0,1)
plt.plot(z,z, marker=",", alpha=1,color='r',)
cb = plt.colorbar(label='count in bin')

plt.xlabel(xlabel.capitalize())
plt.ylabel("Neigbour "+str(xlabel))
plt.text(0.1,0.75,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}\nr={df_master['r'].iloc[0]}\nPearson's r={pearson[0]:.2f}\np={pearson[1]:.2f} ",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)

plt.show()
#%% 
#k=100,d=0,r=3
k=100
D=0.0
r=3
metadata_indices=index_of_simulation(model_type=2,d=D,k=k,r=r)
a=metadata_indices
b=[[],[]]
for i in range(0,len(a),2):#even represents index
   
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
metadata_indices=b

collumns=['file_index','file_name',"time_point","cell_index","age","dead","zposn","A0","k","D","nucl_pos","apical_area","r","neighbour_ids"]
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    df['neighbour_ids']=df['neighbour_ids'].apply(literal_eval)
    df_master=pd.concat([df_master,df])
    del df
for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead' , 'k','D',"r"]:
            df_master[col]=df_master[col].astype('category')
    
df=df_master

#%%
#want final timepoint
df=df_master
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306")]
df_master=df_master.reset_index()

#pairwise cell-neighbour indexes
x=[] #list of cell indexes having neigbours
y=[] #list of neighbour indexes for cell x[i] - 
for i in range(0,len(df_master)):
    a=df_master.cell_index[i]
    b=df_master.neighbour_ids[i]
    if b==[]:
        pass
    else:
        for j in range(0,len(b)):
            x.append(a)
            y.append(b[j])
xlabel='nucl_pos'
x=np.asarray(x)
y=np.asarray(y)
#now nuclear position 2d histogram
a=[]
b=[]
for i in range(0,len(x)):
    a.append(df_master[xlabel][x[i]])
    b.append(df_master[xlabel][y[i]])

pearson=scipy.stats.pearsonr(a,b)
#%%plot
#plt.hexbin(a, b, gridsize=30, cmap='Blues')
fig2, ax=plt.subplots()
plt.hexbin(a, b,gridsize=30, norm=mpl.colors.LogNorm(),  cmap="Blues")
z = np.linspace(0,1)
plt.plot(z,z, marker=",", alpha=1,color='r',)
cb = plt.colorbar(label='count in bin')

plt.xlabel(xlabel.capitalize())
plt.ylabel("Neigbour "+str(xlabel))
plt.text(0.1,0.6,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}\nr={df_master['r'].iloc[0]}\nPearson's r={pearson[0]:.2f}\np={pearson[1]:.2f} ",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)

plt.show()
#%% 
#k=100,d=0.75,r=3
k=100
D=0.75
r=3
metadata_indices=index_of_simulation(model_type=2,d=D,k=k,r=r)
a=metadata_indices
b=[[],[]]
for i in range(0,len(a),2):#even represents index
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
metadata_indices=b

collumns=['file_index','file_name',"time_point","cell_index","age","dead","zposn","A0","k","D","nucl_pos","apical_area","r","neighbour_ids"]
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    
    df['neighbour_ids']=df['neighbour_ids'].apply(literal_eval)
    df_master=pd.concat([df_master,df])
    del df
for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead' , 'k','D',"r"]:
            df_master[col]=df_master[col].astype('category')
df=df_master

#%%
#want final timepoint
df=df_master
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306")]
df_master=df_master.reset_index()

#pairwise cell-neighbour indexes
x=[] #list of cell indexes having neigbours
y=[] #list of neighbour indexes for cell x[i] - 
for i in range(0,len(df_master)):
    a=df_master.cell_index[i]
    b=df_master.neighbour_ids[i]
    if b==[]:
        pass
    else:
        for j in range(0,len(b)):
            x.append(a)
            y.append(b[j])
xlabel='nucl_pos'
x=np.asarray(x)
y=np.asarray(y)
#now nuclear position 2d histogram
a=[]
b=[]
for i in range(0,len(x)):
    a.append(df_master[xlabel][x[i]])
    b.append(df_master[xlabel][y[i]])

pearson=scipy.stats.pearsonr(a,b)
#%%plot
#plt.hexbin(a, b, gridsize=30, cmap='Blues')
fig3, ax=plt.subplots()
plt.hexbin(a, b,gridsize=30, norm=mpl.colors.LogNorm(), cmap="Blues")
z = np.linspace(0,1)
plt.plot(z,z, marker=",", alpha=1,color='r',)
cb = plt.colorbar(label='count in bin')

plt.xlabel(xlabel.capitalize())
plt.ylabel("Neigbour "+str(xlabel))
plt.text(0.1,0.6,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}\nr={df_master['r'].iloc[0]}\nPearson's r={pearson[0]:.2f}\np={pearson[1]:.2f} ",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)

plt.show()
#%%
pp = PdfPages(os.path.join("output","Figure_11.pdf"))
pp.savefig(fig)
pp.savefig(fig1)
pp.savefig(fig2)
pp.savefig(fig3)
pp.close()
#%% figure12
#barchart by d
k=100
r=0
metadata_indices=index_of_simulation(model_type=1,k=k,r=r)
a=metadata_indices
b=[[],[]]
for i in range(0,len(a),2):#even represents index
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
metadata_indices=b

collumns=['file_index','file_name',"time_point","cell_index","dead",'n_neighbours','D']
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    df_master=pd.concat([df_master,df])
    del df
for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead','D']:
    df_master[col]=df_master[col].astype('category')
df=df_master

df_master=df
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306") & (df_master.dead==False)]
df_master=df_master.reset_index()
#%%
palette=['b','orange','green','red','purple']

fig1, ax=plt.subplots()
sns.set_style('whitegrid')
ax= sns.boxplot(x='D',y='n_neighbours',data=df_master)
plt.text(0.1,0.8,s=f"K={k}\nr={r}\n",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.show()
fig2,ax=plt.subplots()
ax= sns.violinplot(x='D',y='n_neighbours',data=df_master)
plt.text(0.1,0.8,s=f"K={k}\nr={r}\n",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.show()
# %%
#now similar plot for r
k=100
d=0
metadata_indices=index_of_simulation(model_type=2,k=k,d=d)
a=metadata_indices
b=[[],[]]
for i in range(0,len(a),2):#even represents index
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
metadata_indices=b

collumns=['file_index','file_name',"time_point","cell_index","dead",'n_neighbours','r']
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    df_master=pd.concat([df_master,df])
    del df
for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead','r']:
    df_master[col]=df_master[col].astype('category')
df=df_master
df_master=df
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306") & (df_master.dead==False)]
df_master=df_master.reset_index()
#%%
fig3, ax=plt.subplots()
sns.set_style('whitegrid')
ax= sns.boxplot(x='r',y='n_neighbours',data=df_master)
plt.text(0.1,0.8,s=f"K={k}\nd={d}\n",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.show()
fig4,ax=plt.subplots()
ax= sns.violinplot(x='r',y='n_neighbours',data=df_master)
plt.text(0.1,0.8,s=f"K={k}\nd={d}\n",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.show()

# %%
pp = PdfPages(os.path.join("output","Figure_12.pdf"))
pp.savefig(fig1)
pp.savefig(fig2)
pp.savefig(fig3)
pp.savefig(fig4)
pp.close()

# %%
#figure 13
#K=100, d=0 r=0.0
k=100
d=0
r=0
metadata_indices=index_of_simulation(model_type=2,k=k,d=d,r=r)
a=metadata_indices
b=[[],[]]
for i in range(0,len(a),2):#even represents index
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
metadata_indices=b
collumns=['file_index','file_name',"time_point","cell_index","dead",'n_neighbours','apical_area']
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    df_master=pd.concat([df_master,df])
    del df
for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead','n_neighbours']:
    df_master[col]=df_master[col].astype('category')
df=df_master

df_master=df
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306") & (df_master.dead==False)]
df_master=df_master.reset_index()
#can drop empty levels of n_neigbours category
df_master['n_neighbours']=df_master['n_neighbours'].astype('int')
#now back
df_master['n_neighbours']=df_master['n_neighbours'].astype('category')
#df_master.n_neighbours.value_counts()
#to check the levels
#%%
sns.set(rc = {'figure.figsize':(6,6)})

fig1, ax=plt.subplots()
sns.set_style('whitegrid')
ax= sns.violinplot(y='apical_area',x='n_neighbours',data=df_master,color="skyblue")
plt.text(0.1,0.8,s=f"K={k}\nd={d}\nr={r}",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.show()

# %%
#K=100, d=0 r=0.0
k=100
d=0.75
r=0
metadata_indices=index_of_simulation(model_type=2,k=k,d=d,r=r)
a=metadata_indices
b=[[],[]]
for i in range(0,len(a),2):#even represents index
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
metadata_indices=b
collumns=['file_index','file_name',"time_point","cell_index","dead",'n_neighbours','apical_area']
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    df_master=pd.concat([df_master,df])
    del df
for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead','n_neighbours']:
    df_master[col]=df_master[col].astype('category')
df=df_master

df_master=df
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306") & (df_master.dead==False)]
df_master=df_master.reset_index()
#can drop empty levels of n_neigbours category
df_master['n_neighbours']=df_master['n_neighbours'].astype('int')
#now back
df_master['n_neighbours']=df_master['n_neighbours'].astype('category')
#df_master.n_neighbours.value_counts()
#to check the levels
#%%
sns.set(rc = {'figure.figsize':(6,6)})

fig2, ax=plt.subplots()
sns.set_style('whitegrid')
ax= sns.violinplot(y='apical_area',x='n_neighbours',data=df_master,color="skyblue")
plt.text(0.1,0.8,s=f"K={k}\nd={d}\nr={r}",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.show()

# %%
#K=100, d=0 r=0.0
k=100
d=0
r=3
metadata_indices=index_of_simulation(model_type=2,k=k,d=d,r=r)
a=metadata_indices
b=[[],[]]
for i in range(0,len(a),2):#even represents index
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
metadata_indices=b
collumns=['file_index','file_name',"time_point","cell_index","dead",'n_neighbours','apical_area']
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    df_master=pd.concat([df_master,df])
    del df
for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead','n_neighbours']:
    df_master[col]=df_master[col].astype('category')
df=df_master

df_master=df
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306") & (df_master.dead==False)]
df_master=df_master.reset_index()
#can drop empty levels of n_neigbours category
df_master['n_neighbours']=df_master['n_neighbours'].astype('int')
#now back
df_master['n_neighbours']=df_master['n_neighbours'].astype('category')
#df_master.n_neighbours.value_counts()
#to check the levels
#%%
sns.set(rc = {'figure.figsize':(6,6)})

fig3, ax=plt.subplots()
sns.set_style('whitegrid')
ax= sns.violinplot(y='apical_area',x='n_neighbours',data=df_master,color="skyblue")
plt.text(0.1,0.8,s=f"K={k}\nd={d}\nr={r}",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.show()

# %%
#K=100, d=0 r=0.0
k=100
d=0.75
r=3
metadata_indices=index_of_simulation(model_type=2,k=k,d=d,r=r)
a=metadata_indices
b=[[],[]]
for i in range(0,len(a),2):#even represents index
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
metadata_indices=b
collumns=['file_index','file_name',"time_point","cell_index","dead",'n_neighbours','apical_area']
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    df_master=pd.concat([df_master,df])
    del df
for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead','n_neighbours']:
    df_master[col]=df_master[col].astype('category')
df=df_master

df_master=df
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306") & (df_master.dead==False)]
df_master=df_master.reset_index()
#can drop empty levels of n_neigbours category
df_master['n_neighbours']=df_master['n_neighbours'].astype('int')
#now back
df_master['n_neighbours']=df_master['n_neighbours'].astype('category')
#df_master.n_neighbours.value_counts()
#to check the levels
#%%
sns.set(rc = {'figure.figsize':(6,6)})

fig4, ax=plt.subplots()
sns.set_style('whitegrid')
ax= sns.violinplot(y='apical_area',x='n_neighbours',data=df_master,color="skyblue")
plt.text(0.1,0.8,s=f"K={k}\nd={d}\nr={r}",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)
plt.show()

#%%
pp = PdfPages(os.path.join("output","Figure_13.pdf"))
pp.savefig(fig1)
pp.savefig(fig2)
pp.savefig(fig3)
pp.savefig(fig4)
pp.close()

#%%
#2d histogram of conditions for division
k=100
d=0.0
r=0
metadata_indices=index_of_simulation(model_type=2,k=k,d=d,r=r)
a=metadata_indices
b=[[],[]]
for i in range(0,len(a),2):#even represents index
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
metadata_indices=b
collumns=['file_index','file_name',"time_point","cell_index","dead",'n_neighbours','apical_area','age','nucl_pos','k','D','r']
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    df_master=pd.concat([df_master,df])
    del df
for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead','n_neighbours','k','D','r']:
    df_master[col]=df_master[col].astype('category')
df=df_master

df_master=df
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306") & (df_master.dead==False)]
df_master=df_master.reset_index()

#%%

#line of critical age
critical_age=t_G1+t_G2+t_S
nuclear_position=0.75
critical_area=A_c
#%%

fig1,ax=plt.subplots()
plt.hexbin(df_master['age'], df_master['nucl_pos'], 
gridsize=30, cmap='Blues')
plt.axhline(y=nuclear_position, color='r', linestyle='-',label="Critical_nuclear_position")
plt.axvline(x=critical_age,color='g',linestyle='-',label="Critical Age")
#plt.plot(x, color="k")
cb = plt.colorbar(label='count in bin')
plt.xlabel('Age')
plt.ylabel('Nucl_pos')
plt.text(0.1,0.6,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}\nr={df_master['r'].iloc[0]}\n ",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)


#plot.ax_joint.axhline(y=nuclear_position,color='r',label="critical_nuclear_position")
ax.legend(bbox_to_anchor=(1.3, 1.05))

plt.show()
#%%

fig2,ax=plt.subplots()
plt.hexbin(df_master['age'], df_master['apical_area'], 
gridsize=30, cmap='Blues')
plt.axhline(y=critical_area, color='r', linestyle='-',label="critical area")
plt.axvline(x=critical_age,color='g',linestyle='-',label="Critical Age")
#plt.plot(x, color="k")
cb = plt.colorbar(label='count in bin')
plt.xlabel('Age')
plt.ylabel("Apical_area")
plt.text(0.1,0.6,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}\nr={df_master['r'].iloc[0]}\n ",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)


#plot.ax_joint.axhline(y=nuclear_position,color='r',label="critical_nuclear_position")
ax.legend(bbox_to_anchor=(1.3, 1.05))

plt.show()


#%%

fig3,ax=plt.subplots()
plt.hexbin(df_master['nucl_pos'], df_master['apical_area'], 
gridsize=30, cmap='Blues')
plt.axhline(y=critical_area, color='r', linestyle='-',label="critical area")
plt.axvline(x=nuclear_position,color='g',linestyle='-',label="critical_nuclear_position")
#plt.plot(x, color="k")
cb = plt.colorbar(label='count in bin')
plt.xlabel('nucl_pos')
plt.ylabel("Apical_area")
plt.text(0.1,0.6,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}\nr={df_master['r'].iloc[0]}\n ",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)


#plot.ax_joint.axhline(y=nuclear_position,color='r',label="critical_nuclear_position")
ax.legend(bbox_to_anchor=(1.3, 1.05))

plt.show()
#%%
#2nd conditions plot
k=100
d=0.75
r=0
metadata_indices=index_of_simulation(model_type=2,k=k,d=d,r=r)
a=metadata_indices
b=[[],[]]
for i in range(0,len(a),2):#even represents index
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
metadata_indices=b
collumns=['file_index','file_name',"time_point","cell_index","dead",'n_neighbours','apical_area','age','nucl_pos','k','D','r']
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    df_master=pd.concat([df_master,df])
    del df
for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead','n_neighbours','k','D','r']:
    df_master[col]=df_master[col].astype('category')
df=df_master

df_master=df
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306") & (df_master.dead==False)]
df_master=df_master.reset_index()

#%%

#line of critical age
critical_age=t_G1+t_G2+t_S
nuclear_position=0.75
critical_area=A_c
#%%

fig4,ax=plt.subplots()
plt.hexbin(df_master['age'], df_master['nucl_pos'], 
gridsize=30, cmap='Blues')
plt.axhline(y=nuclear_position, color='r', linestyle='-',label="Critical_nuclear_position")
plt.axvline(x=critical_age,color='g',linestyle='-',label="Critical Age")
#plt.plot(x, color="k")
cb = plt.colorbar(label='count in bin')
plt.xlabel('Age')
plt.ylabel('Nucl_pos')
plt.text(0.1,0.6,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}\nr={df_master['r'].iloc[0]}\n ",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)


#plot.ax_joint.axhline(y=nuclear_position,color='r',label="critical_nuclear_position")
ax.legend(bbox_to_anchor=(1.3, 1.05))

plt.show()
#%%

fig5,ax=plt.subplots()
plt.hexbin(df_master['age'], df_master['apical_area'], 
gridsize=30, cmap='Blues')
plt.axhline(y=critical_area, color='r', linestyle='-',label="critical area")
plt.axvline(x=critical_age,color='g',linestyle='-',label="Critical Age")
#plt.plot(x, color="k")
cb = plt.colorbar(label='count in bin')
plt.xlabel('Age')
plt.ylabel("Apical_area")
plt.text(0.1,0.6,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}\nr={df_master['r'].iloc[0]}\n ",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)


#plot.ax_joint.axhline(y=nuclear_position,color='r',label="critical_nuclear_position")
ax.legend(bbox_to_anchor=(1.3, 1.05))

plt.show()


#%%

fig6,ax=plt.subplots()
plt.hexbin(df_master['nucl_pos'], df_master['apical_area'], 
gridsize=30, cmap='Blues')
plt.axhline(y=critical_area, color='r', linestyle='-',label="critical area")
plt.axvline(x=nuclear_position,color='g',linestyle='-',label="critical_nuclear_position")
#plt.plot(x, color="k")
cb = plt.colorbar(label='count in bin')
plt.xlabel('nucl_pos')
plt.ylabel("Apical_area")
plt.text(0.1,0.6,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}\nr={df_master['r'].iloc[0]}\n ",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)


#plot.ax_joint.axhline(y=nuclear_position,color='r',label="critical_nuclear_position")
ax.legend(bbox_to_anchor=(1.3, 1.05))

plt.show()

#%%
#3rd conditions plot
k=100
d=0
r=3
metadata_indices=index_of_simulation(model_type=2,k=k,d=d,r=r)
a=metadata_indices
b=[[],[]]
for i in range(0,len(a),2):#even represents index
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
metadata_indices=b
collumns=['file_index','file_name',"time_point","cell_index","dead",'n_neighbours','apical_area','age','nucl_pos','k','D','r']
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    df_master=pd.concat([df_master,df])
    del df
for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead','n_neighbours','k','D','r']:
    df_master[col]=df_master[col].astype('category')
df=df_master

df_master=df
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306") & (df_master.dead==False)]
df_master=df_master.reset_index()

#%%

#line of critical age
critical_age=t_G1+t_G2+t_S
nuclear_position=0.75
critical_area=A_c
#%%

fig7,ax=plt.subplots()
plt.hexbin(df_master['age'], df_master['nucl_pos'], 
gridsize=30, cmap='Blues')
plt.axhline(y=nuclear_position, color='r', linestyle='-',label="Critical_nuclear_position")
plt.axvline(x=critical_age,color='g',linestyle='-',label="Critical Age")
#plt.plot(x, color="k")
cb = plt.colorbar(label='count in bin')
plt.xlabel('Age')
plt.ylabel('Nucl_pos')
plt.text(0.1,0.6,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}\nr={df_master['r'].iloc[0]}\n ",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)


#plot.ax_joint.axhline(y=nuclear_position,color='r',label="critical_nuclear_position")
ax.legend(bbox_to_anchor=(1.3, 1.05))

plt.show()
#%%

fig8,ax=plt.subplots()
plt.hexbin(df_master['age'], df_master['apical_area'], 
gridsize=30, cmap='Blues')
plt.axhline(y=critical_area, color='r', linestyle='-',label="critical area")
plt.axvline(x=critical_age,color='g',linestyle='-',label="Critical Age")
#plt.plot(x, color="k")
cb = plt.colorbar(label='count in bin')
plt.xlabel('Age')
plt.ylabel("Apical_area")
plt.text(0.1,0.6,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}\nr={df_master['r'].iloc[0]}\n ",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)


#plot.ax_joint.axhline(y=nuclear_position,color='r',label="critical_nuclear_position")
ax.legend(bbox_to_anchor=(1.3, 1.05))

plt.show()


#%%

fig9,ax=plt.subplots()
plt.hexbin(df_master['nucl_pos'], df_master['apical_area'], 
gridsize=30, cmap='Blues')
plt.axhline(y=critical_area, color='r', linestyle='-',label="critical area")
plt.axvline(x=nuclear_position,color='g',linestyle='-',label="critical_nuclear_position")
#plt.plot(x, color="k")
cb = plt.colorbar(label='count in bin')
plt.xlabel('nucl_pos')
plt.ylabel("Apical_area")
plt.text(0.1,0.6,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}\nr={df_master['r'].iloc[0]}\n ",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)


#plot.ax_joint.axhline(y=nuclear_position,color='r',label="critical_nuclear_position")
ax.legend(bbox_to_anchor=(1.3, 1.05))

plt.show()
#%%
# 4th conditions plot
k=100
d=0.75
r=3
metadata_indices=index_of_simulation(model_type=2,k=k,d=d,r=r)
a=metadata_indices
b=[[],[]]
for i in range(0,len(a),2):#even represents index
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
metadata_indices=b
collumns=['file_index','file_name',"time_point","cell_index","dead",'n_neighbours','apical_area','age','nucl_pos','k','D','r']
df_master=pd.DataFrame()
for i in metadata_indices[1]:
    a=os.path.join("simulations",str(i[:-3]+"csv"))
    df=pd.read_csv(a,usecols=collumns)
    df_master=pd.concat([df_master,df])
    del df
for col in ['file_index', 'file_name', 'time_point', 'cell_index', 'dead','n_neighbours','k','D','r']:
    df_master[col]=df_master[col].astype('category')
df=df_master

df_master=df
df_master=df_master.loc[(df_master.time_point=="timepoint_306_of_306") & (df_master.dead==False)]
df_master=df_master.reset_index()

#%%

#line of critical age
critical_age=t_G1+t_G2+t_S
nuclear_position=0.75
critical_area=A_c
#%%

fig10,ax=plt.subplots()
plt.hexbin(df_master['age'], df_master['nucl_pos'], 
gridsize=30, cmap='Blues')
plt.axhline(y=nuclear_position, color='r', linestyle='-',label="Critical_nuclear_position")
plt.axvline(x=critical_age,color='g',linestyle='-',label="Critical Age")
#plt.plot(x, color="k")
cb = plt.colorbar(label='count in bin')
plt.xlabel('Age')
plt.ylabel('Nucl_pos')
plt.text(0.1,0.6,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}\nr={df_master['r'].iloc[0]}\n ",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)


#plot.ax_joint.axhline(y=nuclear_position,color='r',label="critical_nuclear_position")
ax.legend(bbox_to_anchor=(1.3, 1.05))

plt.show()
#%%

fig11,ax=plt.subplots()
plt.hexbin(df_master['age'], df_master['apical_area'], 
gridsize=30, cmap='Blues')
plt.axhline(y=critical_area, color='r', linestyle='-',label="critical area")
plt.axvline(x=critical_age,color='g',linestyle='-',label="Critical Age")
#plt.plot(x, color="k")
cb = plt.colorbar(label='count in bin')
plt.xlabel('Age')
plt.ylabel("Apical_area")
plt.text(0.1,0.6,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}\nr={df_master['r'].iloc[0]}\n ",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)


#plot.ax_joint.axhline(y=nuclear_position,color='r',label="critical_nuclear_position")
ax.legend(bbox_to_anchor=(1.3, 1.05))

plt.show()


#%%

fig12,ax=plt.subplots()
plt.hexbin(df_master['nucl_pos'], df_master['apical_area'], 
gridsize=30, cmap='Blues')
plt.axhline(y=critical_area, color='r', linestyle='-',label="critical area")
plt.axvline(x=nuclear_position,color='g',linestyle='-',label="critical_nuclear_position")
#plt.plot(x, color="k")
cb = plt.colorbar(label='count in bin')
plt.xlabel('nucl_pos')
plt.ylabel("Apical_area")
plt.text(0.1,0.6,s=f"K={df_master['k'].iloc[0]}\nD={df_master['D'].iloc[0]}\nr={df_master['r'].iloc[0]}\n ",bbox=dict(facecolor='none', alpha=0.5),transform=ax.transAxes)


#plot.ax_joint.axhline(y=nuclear_position,color='r',label="critical_nuclear_position")
ax.legend(bbox_to_anchor=(1.3, 1.05))

plt.show()



#%%
#saving conditions plots
pp = PdfPages(os.path.join("output","Conditions_plot.pdf"))
pp.savefig(fig1)
pp.savefig(fig2)
pp.savefig(fig3)
pp.savefig(fig4)
pp.savefig(fig5)
pp.savefig(fig6)
pp.savefig(fig7)
pp.savefig(fig8)
pp.savefig(fig9)
pp.savefig(fig10)
pp.savefig(fig11)
pp.savefig(fig12)

pp.close()
#now making animations?
#%%
simulation=index_of_simulation(model_type=1,k=100,r=0,d=0)
metadata_index=simulation[0][0]
experiment_metadata=make_metadata_dic()[metadata_index]
history=make_simulations_dic(metadata_index)
history=list(history.values())
for i in range(0,len(history)):
    cells=history[i]


#%%
fig=plt.figure()
ax = fig.gca()
# initialization function: plot the background of each frame
def init_fig():
    ax = plt.figure()
    return (ax,)
# animation function. This is called sequentially
def animate_fig(i):
    cells_array=history
    v_max = np.max((np.max(cells_array[0].mesh.vertices), np.max(cells_array[-1].mesh.vertices)))
    size = 2.0*v_max
    cells= history[i]
    return model.draw(cells,ax,size)
# call the animator. blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig,animate_fig, init_func=init_fig,
                               frames=len(history))
HTML(anim.to_html5_video())
#FFMpegWriter = manimation.writers['ffmpeg']

# %%
