from curses import meta
from enum import unique
from fileinput import filename
import itertools
from time import time
import numpy as np
import matplotlib.pyplot as plt
import vertex_model as model
import vertex_model.initialisation as init
from vertex_model.forces import TargetArea, Tension, Perimeter, Pressure
import random
import dill
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import os
from vertex_model import simulation_parser
from vertex_model.simulation_parser import *
import matplotlib
import matplotlib.backends.backend_pdf
import plotnine

from plotnine import *
matplotlib.rc('font', family='sansserif', serif='cm10')

matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
random.seed(123)
simulations_dict=simulation_parser.make_simulations_dic(filter=0)
metadata_dic=simulation_parser.make_metadata_dic(simulations_dict)

#list of available simulations


#plot 1

#%%
def plot_1(metadata_index):
    experiment_metadata=metadata_dic[metadata_index]
    new_dict=simulation_parser.slicedict(simulations_dict,experiment_metadata)
    cells=new_dict[list(new_dict)[0]] #picking first key in dictionary and saving its value to cells; ie looking at first timepoint
    #getting parameters of simulation
    cells_k=cells.properties["k"]
    cells_D=cells.properties["D"]
    df=pd.DataFrame()
    age=[]
    target_area=[]
    nuclear_position=[]


    for timepoint in range(0,len(new_dict)):
        cells=new_dict[list(new_dict)[timepoint]]
    #need to find first cell with age 0
        if timepoint==0:
            minimal_age=min(cells.properties["age"])
            cell_index=np.where(cells.properties["age"]==minimal_age)
        age.append(cells.properties['age'][cell_index][0])
        target_area.append(cells.properties['A0'][cell_index][0])
        nuclear_position.append(cells.properties['nucl_pos'][cell_index][0])
    df['age']=age
    df['target_area']=target_area
    df['nuclear_position']=nuclear_position
#plot one: Nucleus position & Target area ~ time
    fig, ax1 = plt.subplots() 
    ax1.set_xlabel('age')
    ax1.set_ylabel('Nuclear_position')
    ax1.yaxis.label.set_color('blue')
    ax1.plot(df['age'], df['nuclear_position'], color = 'blue',label=r"$\mathbf{z^{\alpha}}$")
    ax2 = ax1.twinx() 
    ax2.set_ylabel("Target_area")
    ax2.yaxis.label.set_color("green")
    ax2.plot(df['age'], df['target_area'],label=r"$\mathbf{A^{\alpha}_{0}}$", color = 'green')
    ax1.legend(loc=3)
    ax2.legend(loc=4)
    ax1.text(0.5,0.95,f"k={cells_k}\n D={cells_D}",bbox=dict(facecolor='none', edgecolor='black'))
    #plt.text(0.6,0.9,f"")
    
    return([[fig,ax1,ax2],df])
#%%
 #verifying correct experiments selected
#making threeway vertically stacked plot

#saving plots
#we want d=0.01 & k=5,50,100
#future aaryan add this as a function maybe?
#%%
a=index_of_simulation(model_type=1,k=5,d=0.01,duration=150)+index_of_simulation(model_type=1,k=50,d=0.01,duration=150)+index_of_simulation(model_type=1,k=100,d=0.01,duration=150)
b=[[],[]]
for i in range(0,len(a),2):#even represents index
   
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
plot_1_duration_150_d_0_01=b

pdf = matplotlib.backends.backend_pdf.PdfPages(os.path.join("output","plot_1_duration_150_d_0_01.pdf"))
for i in plot_1_duration_150_d_0_01[0]: #need to manually specify order to make the plots in order of k=100, 50 ,5 with D=0.01
    pdf.savefig(plot_1(i)[0][0])

pdf.close()
#next with duration 100
a=index_of_simulation(model_type=1,k=5,d=0.01,duration=102)+index_of_simulation(model_type=1,k=50,d=0.01,duration=102)+index_of_simulation(model_type=1,k=100,d=0.01,duration=102)
b=[[],[]]
for i in range(0,len(a),2):#even represents index
    
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
plot_1_duration_150_d_0_01=b

pdf = matplotlib.backends.backend_pdf.PdfPages(os.path.join("output","plot_1_duration_102_d_0_01.pdf"))
for i in plot_1_duration_150_d_0_01[0]: #need to manually specify order to make the plots in order of k=100, 50 ,5 with D=0.01
    pdf.savefig(plot_1(i)[0][0])

pdf.close()
#d=0.00
a=index_of_simulation(model_type=1,k=5,d=0.00,duration=150)+index_of_simulation(model_type=1,k=50,d=0.00,duration=150)+index_of_simulation(model_type=1,k=100,d=0.00,duration=150)
b=[[],[]]
for i in range(0,len(a),2):#even represents index
   
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
plot_1_duration_150_d_0_01=b

pdf = matplotlib.backends.backend_pdf.PdfPages(os.path.join("output","plot_1_duration_150_d_0_00.pdf"))
for i in plot_1_duration_150_d_0_01[0]: #need to manually specify order to make the plots in order of k=100, 50 ,5 with D=0.01
    pdf.savefig(plot_1(i)[0][0])

pdf.close()
#next with duration 100
a=index_of_simulation(model_type=1,k=5,d=0.00,duration=102)+index_of_simulation(model_type=1,k=50,d=0.00,duration=102)+index_of_simulation(model_type=1,k=100,d=0.00,duration=102)
b=[[],[]]
for i in range(0,len(a),2):#even represents index
    
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
plot_1_duration_150_d_0_01=b

pdf = matplotlib.backends.backend_pdf.PdfPages(os.path.join("output","plot_1_duration_102_d_0_00.pdf"))
for i in plot_1_duration_150_d_0_01[0]: #need to manually specify order to make the plots in order of k=100, 50 ,5 with D=0.01
    pdf.savefig(plot_1(i)[0][0])

pdf.close()

#%%
#plotting histograms of target_area,cell_cycle length & apical area for treatments of K=100,50,5 with D=0 & r=0
print(metadata_dic) #verifying correct experiments selected
#first make ur dataframe
#%%
def plot_2(metadata_indices):
    df=pd.DataFrame()
    target_area=[]
    nuclear_position=[]
    cells_k=[]
    cells_D=[]
    cell_cycle_lengths=[]
    apical_area=[]
    for experiment in metadata_indices:
        metadata_index=experiment
        experiment_metadata=metadata_dic[metadata_index]
        new_dict=simulation_parser.slicedict(simulations_dict,
        experiment_metadata)
    #need final timepoint
        cells=new_dict[list(new_dict)[-1]]
    #getting parameters of simulation
        cells_k.append([cells.properties["k"]]*len(cells.properties['A0']))
        cells_D.append([cells.properties["D"]]*len(cells.properties['A0']))
        target_area.append(cells.properties['A0'])
        apical_area.append(cells.mesh.area)
        dummy_list=np.asarray([0.1]*len(cells.properties['A0']))
        deadcells = np.where(cells.empty())[0]
        dummy_list[deadcells]=cells.properties['age'][deadcells]
        cell_cycle_lengths.append(dummy_list)

    df['target_area']=list(itertools.chain(*target_area))
    df['cells_k']=list(itertools.chain(*cells_k))
    df['cells_D']=list(itertools.chain(*cells_D))    
    df['cell_cycle_lengths']=list(itertools.chain(*cell_cycle_lengths))
    df['apical_area']=list(itertools.chain(*apical_area))
    #plot
    df_temp=df[df['target_area'].between(0,3)]
    #filtering target area between 0 &3
    
    histogram=ggplot(data=df_temp,mapping=aes(x='target_area',color='factor(cells_k)',y=plotnine.after_stat('density')))
    a=(histogram+plotnine.geom_histogram(alpha=0.5,fill=None)+ plotnine.theme(text=plotnine.element_text(family="sans-serif"))+plotnine.annotate('text', x=2.9, y = 4, label=f'D={cells_D[0][0]}\n r={0}')+plotnine.labs(color="K")+xlim(0,3)   
)
    df_temp=df[df['cell_cycle_lengths'].between(0.8,1.5)]
    histogram=ggplot(data=df_temp,mapping=aes(x='cell_cycle_lengths',color='factor(cells_k)',y=plotnine.after_stat('density')))
    b=(histogram+plotnine.geom_histogram(alpha=0.5,fill =None)+plotnine.theme(text=plotnine.element_text(family='sans-serif'))
    +plotnine.annotate('text',x=1.3,y=4,label=f'D={cells_D[0][0]}\n r={0}')+plotnine.labs(color="K")+xlim(0.8,1.5)
    )
    df_temp=df[df['apical_area'].between(0.0000001,1)]
    histogram=ggplot(data=df_temp,mapping=aes(x='apical_area',color='factor(cells_k)',y=plotnine.after_stat('density')))
    c=(histogram+plotnine.geom_histogram(alpha=0.5,fill=None)+plotnine.theme(text=plotnine.element_text(family='sans-serif'))+plotnine.annotate('text',x=0.5,y=10,label=f'D={cells_D[0][0]} \n r={0}')+plotnine.labs(colorl="K")+xlim(0,1))
    print(a)
    print(b)
    print(c)
    yield a
    yield b
    yield c
    return [a,b,c,df]

#%%
#need k=100,50,5 for d=0 and d=0.01 and duration 102 and 150

#d=0.0 duration=150
a=index_of_simulation(model_type=1,k=5,d=0.00,duration=150)+index_of_simulation(model_type=1,k=50,d=0.00,duration=150)+index_of_simulation(model_type=1,k=100,d=0.00,duration=150)
b=[[],[]]
for i in range(0,len(a),2):#even represents index
   
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
plot_2_duration_150_d_0_00=b
del a
del b
save_as_pdf_pages(plot_2(plot_2_duration_150_d_0_00[0]),filename="plot_2_D=0.01_duration=150.pdf",path=os.path.join("output"))
#d=0.01 duration=150
a=index_of_simulation(model_type=1,k=5,d=0.01,duration=150)+index_of_simulation(model_type=1,k=50,d=0.01,duration=150)+index_of_simulation(model_type=1,k=100,d=0.01,duration=150)
b=[[],[]]
for i in range(0,len(a),2):#even represents index
   
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
plot_2_duration_150_d_0_01=b
del a
del b

save_as_pdf_pages(plot_2(plot_2_duration_150_d_0_01[0]),filename="plot_2_D=0.01_duration=150.pdf",path=os.path.join("output"))
#d=0.0 duration=102


a=index_of_simulation(model_type=1,k=5,d=0.0,duration=102)+index_of_simulation(model_type=1,k=50,d=0.0,duration=102)+index_of_simulation(model_type=1,k=100,d=0.0,duration=102)
b=[[],[]]
for i in range(0,len(a),2):#even represents index
   
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
plot_2_duration_102_d_0_0=b
del a
del b


save_as_pdf_pages(plot_2(plot_2_duration_102_d_0_0[0]),filename="plot_2_D=0.00_duration=102.pdf",path=os.path.join("output"))
#d=0.01 duration=102
a=index_of_simulation(model_type=1,k=5,d=0.01,duration=102)+index_of_simulation(model_type=1,k=50,d=0.01,duration=102)+index_of_simulation(model_type=1,k=100,d=0.01,duration=102)
b=[[],[]]
for i in range(0,len(a),2):#even represents index
   
    b[0].append(a[i])
    b[1].append(a[i+1])
b[0]=flatten(b[0])
b[1]=flatten(b[1])
plot_2_duration_102_d_0_01=b
del a
del b


save_as_pdf_pages(plot_2(plot_2_duration_102_d_0_01[0]),filename="plot_2_D=0.01_duration=102.pdf",path=os.path.join("output"))
#%%

#plot of crowding force 
#first extract all model 2 dicts
simulations_dict=simulation_parser.make_simulations_dic(filter=2)
metadata_dic=simulation_parser.make_metadata_dic(simulations_dict)
def plot_crowding(metadata_indices):
    df=pd.DataFrame()
    target_area=[]
    nuclear_position=[]
    cells_k=[]
    cells_D=[]
    cells_r=[]
    cell_cycle_lengths=[]
    apical_area=[]
    for experiment in metadata_indices:
        metadata_index=experiment
        experiment_metadata=metadata_dic[metadata_index]
        new_dict=simulation_parser.slicedict(simulations_dict,
        experiment_metadata)
    #need final timepoint
        cells=new_dict[list(new_dict)[-1]]
    #getting parameters of simulation
        cells_k.append([cells.properties["k"]]*len(cells.properties['A0']))
        cells_D.append([cells.properties["D"]]*len(cells.properties['A0']))
        target_area.append(cells.properties['A0'])
        apical_area.append(cells.mesh.area)
        dummy_list=np.asarray([0.1]*len(cells.properties['A0']))
        deadcells = np.where(cells.empty())[0]
        dummy_list[deadcells]=cells.properties['age'][deadcells]
        cell_cycle_lengths.append(dummy_list)

    df['target_area']=list(itertools.chain(*target_area))
    df['cells_k']=list(itertools.chain(*cells_k))
    df['cells_D']=list(itertools.chain(*cells_D))    
    df['cell_cycle_lengths']=list(itertools.chain(*cell_cycle_lengths))
    df['apical_area']=list(itertools.chain(*apical_area))
    df['cells_r']=([0.0]*2102)+([0.1]*2102)+([0.3]*2102)+([1.0]*2102)+([3.0]*2102)#CHANGE IMM EDIATLEY
    #plot
    histogram=ggplot(data=df,mapping=aes(x='cell_cycle_lengths',fill='factor(cells_r)',y=plotnine.after_stat('density')))
    b=(histogram+plotnine.geom_histogram(alpha=0.5)+plotnine.theme(text=plotnine.element_text(family='sans-serif'))
    +plotnine.annotate('text',x=1.3,y=4,label=f'D={cells_D[0][0]}\n k={cells_k[0][0]}')+plotnine.labs(fill="r")
    )



plt.hist(df['r']


)

metadata_indices=range(0,5)
df.pivot(collumns='cells_r',values='cell_cycle_lengths').plot.hist(bins=100)
plt.show()
    #here
    metadata_index=3
experiment_metadata=metadata_dic[metadata_index]
new_dict=simulation_parser.slicedict(simulations_dict,
        experiment_metadata)
    #need final timepoint
cells=new_dict[list(new_dict)[-1]]
print(len(cells.properties['age']))
# %%
