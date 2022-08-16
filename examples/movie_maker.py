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
import sys
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


random.seed(123)
#%%
metadata_index=int(sys.argv[1])
#%%
test_file=make_metadata_dic()[metadata_index]
print(test_file)
history=make_simulations_dic(metadata_index)
metadata=history
file_name=list(metadata)[0]
file_k=re.search(r'_k_(.*)_D_',file_name).group(1)
file_D=re.search(r'_D_(.*)_a_',file_name).group(1)   
file_a=re.search(r'_a_(.*)_run_',file_name).group(1)

history=list(history.values())

   

#%%
fig=plt.figure(figsize=[15,15])
ax = fig.gca()
# initialization function: plot the background of each frame
def init_fig():
    ax = plt.figure(figsize=[15,15])
    return (ax,)
# animation function. This is called sequentially
def animate_fig(i):
    file_timepoint=f"timepoint {i} out of {len(history)}"
    string_title=f"k={file_k}\nD={file_D}\na={file_a}\n{file_timepoint}"
    cells_array=history
    v_max = np.max((np.max(cells_array[0].mesh.vertices), np.max(cells_array[-1].mesh.vertices)))
    size = 2*v_max
    cells= history[i]
  
    
    return model.draw(cells,ax,size,string_for_title=string_title)
    
    
# call the animator. blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig,animate_fig, init_func=init_fig,
                               frames=len(history))

writervideo = animation.FFMpegWriter(fps=60) 
anim.save(os.path.join("output",file_name[:-23]+".mp4"), writer=writervideo)

#FFMpegWriter = manimation.writers['ffmpeg']


# %%
