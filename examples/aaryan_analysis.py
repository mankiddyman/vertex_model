from curses import meta
from enum import unique
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
import matplotlib
import matplotlib.backends.backend_pdf
import plotnine
from plotnine import ggplot, aes
matplotlib.rc('font', family='sansserif', serif='cm10')

matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
random.seed(123)


#plot 1
simulations_dict=simulation_parser.make_simulations_dic(filter=1)
metadata_dic=simulation_parser.make_metadata_dic(simulations_dict)
# we want k=100 and D=0.01

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
print(metadata_dic) #verifying correct experiments selected
#making threeway vertically stacked plot

#saving plots
#we want d=0.01 & k=5,50,100
pdf = matplotlib.backends.backend_pdf.PdfPages(os.path.join("output","plot_1.pdf"))
for i in [0,2,4]:
    pdf.savefig(plot_1(i)[0][0])

pdf.close()

#plotting histograms of target_area,cell_cycle length & apical area for treatments of K=100,50,5 with D=0 & r=0
print(metadata_dic) #verifying correct experiments selected
#first make ur dataframe
def plot_2(metadata_indices)
    df=pd.DataFrame()
    target_area=[]
    nuclear_position=[]
    cells_k=[]
    cells_D=[]
    for experiment in metadata_indices:
        metadata_index=experiment
        experiment_metadata=metadata_dic[metadata_index]
        new_dict=simulation_parser.slicedict(simulations_dict,
        experiment_metadata)
    #need final timepoint
        cells=new_dict[list(new_dict)[-1]]
    #getting parameters of simulation
        cells_k.append([cells.properties["k"]]*len(cells.properties['A0']))
        cells_D.append([cells.properties["D"]*len(cells.properties['A0'])])
        target_area.append(cells.properties['A0'])
    df['target_area']=target_area
    
    #plot
    histogram=ggplot(data=df,mapping=aes(x='target_area'))
    (histogram+plotnine.geom_histogram()+ plotnine.theme(text=plotnine.element_text(family="sans-serif")))

plt.figure()
newplot=sns.histplot(data=df,x='A0',bins=200)
plt.text(2.5, 800,f"k={cells_k} \n D={cells_D}", fontsize=9) #add text)
#plotting cell cycle length histogram

deadcells = np.where(cells.empty())[0]
cell_cycle_lengths = cells.properties['age'][deadcells]
plt.figure()
plt.hist(cell_cycle_lengths, bins = 50, density=True)
plt.title('Cell cycle length distribution Model 1')
plt.text(2.8,1,f"k={cells_k} \n D={cells_D} \n r={cells_r}")
plt.show()
# plt.figure()
# newplot=sns.histplot(data=df,x='age',)
# plt.text(1,200 ,f"k={cells_k} \n D={cells_D}", fontsize=9) #add text)
# plt.xlim(0.8,1.5)



# %%
