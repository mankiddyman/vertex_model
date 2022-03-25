import itertools
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
with open("model_1_k100.0_D0.01_run_1.pkl", "rb") as file:
    history = dill.load(file)

cells = history[-1]
random.seed(123)

#for reference of properties of cells
#print(cells.properties.keys())
#plot of isolated nucleus position against

N_of_cells=len(cells.properties['age'])
#pick a random cell
cell_of_interest=random.randrange(0,N_of_cells-1)

#making dataframe from dictionary
for i in cells.properties:
    print(i)
df=pd.DataFrame(cells.properties['age'],columns=["age"])
df['zposn']=cells.properties['zposn']
df['nucl_pos']=cells.properties["nucl_pos"]
df['A0']=cells.properties['A0']
df['age']=cells.properties['age']
cells_k=cells.properties['k']
cells_D=cells.properties['D']
if file.name.startswith("model_1"):
    cells_r=0
#sns.relplot(data=df,x='age',y='nucl_pos',s=5)
#sns.relplot(data=df,x='age',y='A0',s=5)
sns.lineplot(data=df,x='age',y='nucl_pos', color="b",)
ax2 = plt.twinx()
sns.lineplot(data=df,x='age',y='A0', color="g", ax=ax2)
plt.xlim(0,1)
plt.ylim(0.6,2)
plt.text(0.5, 2,f"k={cells_k} \n D={cells_D}", fontsize=9) #add text)

#moving on

#plotting target area histogram
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
