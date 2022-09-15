
# # All the function to run a simulation

#Import libraries

#########################################################################PARALELL################################################

#########################################################################PARALELL################################################
# %matplotlib tk 
#get_ipython().magic(u'matplotlib') #to use model.animate and see video alive
import itertools
import multiprocessing
import time 
import random
from itertools import repeat
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import vertex_model as model
import vertex_model.initialisation as init
from vertex_model.forces import TargetArea, Tension, Perimeter, Pressure
import os
import math
import seaborn as sns
import warnings
warnings.filterwarnings('ignore') #Don't show warnings
from vertex_model.Gobal_Constant import dt, viscosity, t_G1, t_G2, t_S, A_c, J, pos_d, T1_eps, P, microns, time_hours, expansion_constant #file with necessary constants
import ipdb
diff_rate_hours=0.05 


#differentiation rate (1/h) 



def check_estage():
    print("Running a %s hours"%J + " %s"%pos_d_v)
    print("dt=%s"%dt)            #time step
    print("viscosity=%s" %viscosity)  #viscosity*dv/dt = F
    print("A_c=%s"%A_c) #critical area
    print("T1_eps =%s"%T1_eps)
    

# run simulation - input is a generator = simulation => returns array
def run(simulation,N_step,skip):
    return [cells.copy() for cells in itertools.islice(simulation,0,int(N_step),int(skip))] # itertools.islice looks through the generator's elements. => list of cell objects = the history



def division_axis(mesh,face_id,rand):
    """Choose a random division axis (given as a pair of boundary edges to bisect) for the given cell.
    
    The first edge is chosen randomly from the bounding edges of the cell with probability proportional 
    to edge length. The second edge is then fixed to be n_edge/2 from the first. 
    """
    edges = mesh.boundary(face_id)
    if edges==[-1]:
        print('here')
        os._exit(1)
    p = np.cumsum(mesh.length[edges])
    e0 = p.searchsorted(rand.rand()*p[-1])
    return edges[e0],edges[e0-len(edges)//2]  

def bin_by_xpos(cells,percentiles):
    vx = cells.mesh.vertices[0]
    #simple 'midpoint' as mean of vertex positions
    # np.bincount = Count number of occurrences of each value in array of non-negative ints.
    # 1.0, 0.5; mid_x / width?
    mid_x = np.bincount(cells.mesh.face_id_by_edge,weights=vx)
    counts = np.maximum(np.bincount(cells.mesh.face_id_by_edge),1.0)
    mid_x = mid_x / counts 
    width = cells.mesh.geometry.width
    return np.searchsorted(percentiles,(mid_x/width + 0.5) % 1.0)   

#simulation without division - thermalisation without ages of cells and INM - evolves a mesh
def basic_simulation(cells,force,dt=dt,T1_eps=0.04): # a generator
    while True:
        cells.mesh , number_T1 = cells.mesh.transition(T1_eps)
        F = force(cells)/viscosity
        expansion = 0.05*np.average(F*cells.mesh.vertices,1)*dt
        dv = dt*model.sum_vertices(cells.mesh.edges,F) 
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+ expansion)
        yield cells # outputs the cells object but continues to run (unlike return)!

# simulation with division and INM (no differentiation rate domain) # type 100
def simulation_with_division(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None): #(cells,force,dt=0.001,T1_eps=0.04,lifespan=100.0,rand=Non
    properties = cells.properties
    properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    properties['ageingrate'] = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells)) #degradation rate per each cell
    properties['ids_division'] = [] #save ids of the cell os the division when it's ready per each time step
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pD'] = []
    expansion = np.array([0.0,0.0])
    while True:
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        ready = np.where(~cells.empty() & (cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)))[0]  
        if len(ready): #these are the cells ready to undergo division at the current timestep
            # print(ready)
            properties['ageingrate'] =np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2*len(ready))))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  # Daugthers and parent have the same ids
            properties['ids_division'] = ready
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready] #New edges after division 
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                properties['Division_angle_pD']= np.append(properties['Division_angle_pD'],cells.mesh.edge_angle[division_new_edge][0])
        properties['age'] = properties['age']+dt*properties['ageingrate'] #add time step depending of the degradation rate 
        
        #Calculate z nuclei position (Apical-Basal movement), depending of the cell cycle phase time and age of the cell
        N_G1=1-1.0/t_G1*properties['age'] 
        N_S=0 
        N_G2=1.0/(t_G2)*(properties['age']-(t_G1+t_S))
        properties['zposn'] = np.minimum(1.0,np.maximum(N_G1,np.maximum(N_S,N_G2)))

        
        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2)
        
        #############BAETTI VID
        
        #ath BAETTI D VID
        cells.mesh , number_T1, d= cells.mesh.transition(T1_eps)  #check edges verifing T1 transition
        if len(number_T1)>0:
            for ii in number_T1:
                index = cells.mesh.face_id_by_edge[ii]
                properties['T1_angle_pD'] =np.append(properties['T1_angle_pD'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 

        #EG BAETTI VID
        len_modified=np.matrix.copy(cells.mesh.length)
        #len_modified[np.where((properties['parent_group'][cells.mesh.face_id_by_edge]==1) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==0))]*=0.12
        
        tens = (0.5*cells.by_edge('Lambda', 'Lambda_boundary')/len_modified)*cells.mesh.edge_vect
        tens= tens - tens.take(cells.mesh.edges.prev, 1)
        F+=tens

        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        properties['force_x'] = F[0]*viscosity
        properties['force_y'] = F[1]*viscosity   
        
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        yield cells 


def pairwise_crowding_force (delta_z, a=1.0, s=0.1):
    '''
    Calculates the pairwise contributions for the crowding force as a 
    function of the difference in AB nuclear position between neighbouring
    cells.
    It accepts as input an array (and the function is applied element-wise)
    '''
    y = delta_z/s
    return a * np.exp(-y**2/2.)*y

def crowding_force (cells, a=1.0, s=0.1):
    '''
    For every cell calculate the total crowding force
    as the sum over its neighbours of the pairwise crowding force

    '''

    # Get the positions of all the nuclei.
    nucl_pos = cells.properties['nucl_pos']

    # Get the ids of cells by their edges
    cell_ids = cells.mesh.face_id_by_edge   # ids of faces/cells
    neig_ids = cell_ids[cells.mesh.edges.reverse] # ids of their neighbours

    # position of nuclei in cells and their neighbours.
    z  = np.take(nucl_pos, cell_ids) # cells of interest
    zn = np.take(nucl_pos, neig_ids) # their neighbours

    # calculate the difference in nucl_pos between neighbouring cells.
    delta_z = z - zn   # the difference in z for all pairs
    f_pairs = pairwise_crowding_force(delta_z, a=a, s=s)  # the corresponding pairwise forces

    # calculate the sum of the pairwise contributions for each cell of interest
    force = np.bincount(cell_ids, f_pairs, cells.mesh.n_face)

    return force

def crowding_force_2(cells,a=0.1,s=0.1):
    """
    Modification of old crowding force equation to include interactions between non adjacent neighbours by using their x and y coordinates to position them
    we will only consider interactions with cells within 3 distance units of eachother
    we will also scale the xy scale to the same units and z
    1 Area Unit in simulation is 23 micrometers^2 in real life
    it follows that one 1 length unit in the x-y plane is sqrt(23) micrometers in real life 

    for the z axis we have to consider more carefully
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4528075/
    the following paper states the nuclei as having diameter 8 micrometers

    the thickness of the mouse neuroepithelium (A-B axis length) is a bit difficult for me to find precisely but this paper
    https://www.frontiersin.org/articles/10.3389/fcell.2016.00139/full
    says 50-100 micrometers

    we will use a value of 80 micrometers as the length of the axis and therefore s (the size of nucleus relative to z-axis) is a neat 0.1

    so to scale x-y plane to z

    we will multiply x-y coordinates by sqrt(23) and divide by 80 to get values in z coordinates
    we return an array containing the aggregate crowding forces for all cells
    """
    #first going to make a vector containing all x,y and z positions

    [x_array,y_array]=cells.mesh.centres()
    #conversion ofx-y
    x_array=x_array*np.sqrt(23)/80
    y_array=y_array*np.sqrt(23)/80
    
    z_array=cells.properties['nucl_pos']

    #now going to create a distance matrix for x,y and z
    #one matrix x_1 is x_array as the rows of a square matrix
    #one matrix x_2 is x_array as the collumn of a square matrix
    x_1,x_2=np.meshgrid(x_array,x_array)
    y_1,y_2=np.meshgrid(y_array,y_array)
    z_1,z_2=np.meshgrid(z_array,z_array)
    #x_2=np.tile(np.array([x_array]).transpose(),(1,n_cells))
    #y_1=np.tile(y_array,(n_cells,1))
    #y_2=np.tile(np.array([y_array]).transpose(),(1,n_cells))
    #z_1=np.tile(z_array,(n_cells,1))
    #z_2=np.tile(np.array([z_array]).transpose(),(1,n_cells))
    
    #if we do x_1 - x_2 = x_matrix we get a matrix which contains the distance in the x-plane (with direction) from cell A to cell B by looking at coordinates [A,B] in x_matrix

    x_matrix=x_1-x_2
    y_matrix=y_1-y_2
    z_matrix=z_1-z_2
    #example of x componenet of vector connecting cell 1869 to cell 1879 
    #x_matrix[1869,1879]
    #block begins here
    #now we want to make a euclidian distance matrix #NOTE THIS IS JUST THE MAGNITUDE OF THE DISTANCE AND SIGNS ARE ALL POSITIVE
    
    pieces=matrix_to_pieces([x_matrix,y_matrix,z_matrix])
    iterable_for_multiprocessing=[None]*len(pieces[0])
    for i in range(0,len(pieces[0])):
        iterable_for_multiprocessing[i]=[pieces[0][i],pieces[1][i],pieces[2][i],a,s]

    forces=multiprocess_euc_to_force(iterable_for_multiprocessing)
    #block begins here
    force_matrix=forces[0]
    for i in range(1,len(pieces[0])):
        force_matrix=np.concatenate((force_matrix,forces[i]),axis=0)
    
    forces_array=np.nansum(force_matrix,axis=1)
    
    return forces_array

def matrix_to_pieces(x_y_z_matrix_list):
    x_matrix=x_y_z_matrix_list[0]
    y_matrix=x_y_z_matrix_list[1]
    z_matrix=x_y_z_matrix_list[2]
    n_cores=12
    #we split the matrix into length/n_cores pieces with remainder of rows being in the last process
    length_pieces=[int(len(x_matrix)/n_cores)]*n_cores
    length_pieces[-1]=length_pieces[-1]+len(x_matrix)%n_cores
    length_pieces=np.cumsum(length_pieces)
    length_pieces=np.insert(length_pieces,obj=0,values=0)
    x_pieces=np.empty((n_cores, 0)).tolist()
    y_pieces=np.empty((n_cores, 0)).tolist()
    z_pieces=np.empty((n_cores, 0)).tolist()
    for i in range(0,n_cores):
        x_pieces[i]=x_matrix[length_pieces[i]:length_pieces[i+1],0:len(x_matrix)]
        y_pieces[i]=y_matrix[length_pieces[i]:length_pieces[i+1],0:len(y_matrix)]
        z_pieces[i]=z_matrix[length_pieces[i]:length_pieces[i+1],0:len(z_matrix)]
    return[x_pieces,y_pieces,z_pieces]
    
def euc_to_force(x_y_z_pieces_a_s_list):
    x=x_y_z_pieces_a_s_list[0]
    y=x_y_z_pieces_a_s_list[1]
    z=x_y_z_pieces_a_s_list[2]
    a=x_y_z_pieces_a_s_list[3]
    s=x_y_z_pieces_a_s_list[4]
    euc_matrix=np.sqrt(x**2+y**2+z**2)
    #to reduce number of calculations we will reduce the number of cells that can interact by setting interactions to only occur below 0.3 z-units or 5.00 x-y units
    euc_matrix[euc_matrix < 0.15]
    #
    force_matrix=-a*np.exp(-(euc_matrix**2/s**2))*z/euc_matrix
    return force_matrix

def multiprocess_euc_to_force(iterable_for_multiprocessing):
    with multiprocessing.Pool() as pool:
        force_pieces=[]
        for result in pool.map(euc_to_force,iterable_for_multiprocessing):
            force_pieces.append(result)
        return force_pieces


# simulation with division and INM (no differentiation rate domain) # type 0 Rebeca model 1 -> type 4: delayed drift + noise
def simulation_with_division_model_1(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None):
    properties = cells.properties
    properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    properties['ageingrate'] = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells)) #degradation rate per each cell
    #dummy = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells))
    #properties['ageingrate'] = 1.0/lifespan*np.ones(len(cells)) #degradation rate per each cell
    properties['ids_division'] = [] #save ids of the cell os the division when it's ready per each time step
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pD'] = []
    properties['nucl_pos'] = properties['zposn'].copy()

    expansion = np.array([0.0,0.0])

    
    while True:
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        ready = np.where(~cells.empty() & (cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)))[0]  # divides if nucleus pos > 0.75
        properties['ageingrate'][ready] = 0
        if len(ready): #these are the cells ready to undergo division at the current timestep
            # print("cells dividing --> ", ready)
            properties['ageingrate'] = np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2*len(ready))))
            #dummy = np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2*len(ready)))
            #properties['ageingrate'] =np.append(properties['ageingrate'], 1.0/lifespan*np.ones(2*len(ready)))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  # Daugthers and parent have the same ids

            '''
            2 possibilities for nucl_pos of daughter cells:
            - nucleus at z = 1
            - nucleus at same z as the parent (more realistic)
            '''
            # properties['nucl_pos'] = np.append(properties['nucl_pos'], np.ones(2*len(ready)))
            properties['nucl_pos'] = np.append(properties['nucl_pos'], np.repeat(properties['nucl_pos'][ready],2))

            properties['ids_division'] = ready
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready] #New edges after division 
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                properties['Division_angle_pD']= np.append(properties['Division_angle_pD'],cells.mesh.edge_angle[division_new_edge][0])
        properties['age'] = properties['age']+dt*properties['ageingrate'] #add time step depending of the degradation rate 
        

        #Calculate z nuclei position (Apical-Basal movement), depending of the cell cycle phase time and age of the cell
        N_G1=1-1.0/t_G1*properties['age']
        N_S=0
        N_G2=1.0/(t_G2)*(properties['age']-(t_G1+t_S))
        properties['zposn'] = np.minimum(1.0,np.maximum(N_G1,np.maximum(N_S,N_G2)))

        properties['nucl_pos'] += properties['k']*(properties['zposn'] - properties['nucl_pos'])*dt \
                                  + np.sqrt(2*properties['D']*dt)*np.random.randn(len(properties['nucl_pos']))
        
        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+properties['nucl_pos']**2) # target area now depends on nucl_pos
        properties['A0_initial'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2) # initial
        
        
        #############BAETTI VID
        
        #ath BAETTI D VID
        cells.mesh , number_T1, d = cells.mesh.transition(T1_eps)  # check edges verifing T1 transition
        if len(number_T1)>0:
            for ii in number_T1:
                index = cells.mesh.face_id_by_edge[ii]
                properties['T1_angle_pD'] =np.append(properties['T1_angle_pD'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 

        #EG BAETTI VID
        len_modified=np.matrix.copy(cells.mesh.length)
        #len_modified[np.where((properties['parent_group'][cells.mesh.face_id_by_edge]==1) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==0))]*=0.12
        
        tens = (0.5*cells.by_edge('Lambda', 'Lambda_boundary')/len_modified)*cells.mesh.edge_vect
        tens= tens - tens.take(cells.mesh.edges.prev, 1)
        F+=tens

        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        properties['force_x'] = F[0]*viscosity
        properties['force_y'] = F[1]*viscosity   
        
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        yield cells
       
        
# simulation with division and INM (no differentiation rate domain) # type 0 Rebeca model 1 -> type 5: drift = v0 & noise
def simulation_with_division_model_2(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None): #(cells,force,dt=0.001,T1_eps=0.04,lifespan=100.0,rand=Non
    properties = cells.properties
    properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    properties['ageingrate'] = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells)) #degradation rate per each cell
    properties['ids_division'] = [] #save ids of the cell os the division when it's ready per each time step
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pD'] = []

    properties['nucl_pos'] = properties['zposn'].copy() # nucleus pos with drift = v0 and noise
    
    expansion = np.array([0.0,0.0])

    
    while True:
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        ready = np.where(~cells.empty() & (cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)))[0]  # divides if nucleus pos > 0.75
        properties['ageingrate'][ready] = 0
        if len(ready): #these are the cells ready to undergo division at the current timestep
            # print("cells dividing --> ", ready)
            properties['ageingrate'] =np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2*len(ready))))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  # Daugthers and parent have the same ids
            
            '''
            2 possibilities for nucl_pos of daughter cells:
            - nucleus at z = 1
            - nucleus at same z as the parent (more realistic)
            '''
            # properties['nucl_pos'] = np.append(properties['nucl_pos'], np.ones(2*len(ready)))
            properties['nucl_pos'] = np.append(properties['nucl_pos'], np.repeat(properties['nucl_pos'][ready],2))

            properties['ids_division'] = ready
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready] #New edges after division 
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                properties['Division_angle_pD']= np.append(properties['Division_angle_pD'],cells.mesh.edge_angle[division_new_edge][0])
        properties['age'] = properties['age']+dt*properties['ageingrate'] #add time step depending of the degradation rate 

        #Calculate z nuclei position (Apical-Basal movement), depending of the cell cycle phase time and age of the cell
        N_G1=1-1.0/t_G1*properties['age']
        N_S=0 
        N_G2=1.0/(t_G2)*(properties['age']-(t_G1+t_S))
        properties['zposn'] = np.minimum(1.0,np.maximum(N_G1,np.maximum(N_S,N_G2)))

        v0 = np.zeros_like(properties['age']) # array for velocities
        G1_cells = np.where((0<=properties['age']) & (properties['age']<=t_G1))[0]
        v0[G1_cells] = -1.0/t_G1
        #S_cells = np.where(t_G1 < properties['age'] <= t_G1 + t_S)[0]
        #v0[S_cells] = 0
        G2_cells = np.where((t_G1+t_S<properties['age']) & (properties['age']<=t_G1+t_S+t_G2))[0]
        v0[G2_cells] = 1.0/t_G2
        #M_cells = np.where(t_G1+t_S+t_G2 < properties['age'])
        #v0[M_cells] = 0

        properties['nucl_pos'] += v0*dt + np.sqrt(2*properties['D']*dt)*np.random.randn(len(properties['zposn']))

        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+properties['nucl_pos']**2) # target area now depends on nucl_pos_2
        properties['A0_initial'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2) # rerun simul with this
        
    
        #############BAETTI VID
        
        #ath BAETTI D VID
        cells.mesh , number_T1, d = cells.mesh.transition(T1_eps)   #  #check edges verifing T1 transition
        if len(number_T1)>0:
            for ii in number_T1: #
                index = cells.mesh.face_id_by_edge[ii]
                properties['T1_angle_pD'] =np.append(properties['T1_angle_pD'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 

        #EG BAETTI VID
        len_modified=np.matrix.copy(cells.mesh.length)
        #len_modified[np.where((properties['parent_group'][cells.mesh.face_id_by_edge]==1) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==0))]*=0.12
        
        tens = (0.5*cells.by_edge('Lambda', 'Lambda_boundary')/len_modified)*cells.mesh.edge_vect
        tens= tens - tens.take(cells.mesh.edges.prev, 1)
        F+=tens

        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        properties['force_x'] = F[0]*viscosity
        properties['force_y'] = F[1]*viscosity   
        
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        yield cells
               
# simulation with division and INM (no differentiation rate domain) # type 0 Rebeca model 3 -> type 6: delayed drift + noise + force
def simulation_with_division_model_3(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None): #(cells,force,dt=0.001,T1_eps=0.04,lifespan=100.0,rand=Non
    properties = cells.properties
    properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    properties['ageingrate'] = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells)) #degradation rate per each cell
    properties['ids_division'] = [] #save ids of the cell os the division when it's ready per each time step
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pD'] = []

    properties['nucl_pos'] = properties['zposn'].copy()
    
    properties['force_z'] = []
    
    expansion = np.array([0.0,0.0])
    
    while True:
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        ready = np.where(~cells.empty() & (cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)))[0]  # divides if nucleus pos > 0.75
        if len(ready): #these are the cells ready to undergo division at the current timestep
            properties['ageingrate'] =np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2*len(ready))))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  # Daugthers and parent have the same ids
            # daughter cells = same positions as parents - do this NEXT
            properties['nucl_pos'] = np.append(properties['nucl_pos'], np.ones(2*len(ready)))

            properties['ids_division'] = ready
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready] #New edges after division 
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                properties['Division_angle_pD']= np.append(properties['Division_angle_pD'],cells.mesh.edge_angle[division_new_edge][0])
        #properties['age'] = properties['age']+dt*properties['ageingrate'] #add time step depending of the degradation rate 

        # IMPORTANT: add age ingrate only for alive cells
        alivecells = np.where(~cells.empty())[0]
        properties['age'][alivecells] += dt*properties['ageingrate'][alivecells] #add time step depending of the degradation rate 

        #Calculate z nuclei position (Apical-Basal movement), depending of the cell cycle phase time and age of the cell
        N_G1=1-1.0/t_G1*properties['age']
        N_S=0 
        N_G2=1.0/(t_G2)*(properties['age']-(t_G1+t_S))
        properties['zposn'] = np.minimum(1.0,np.maximum(N_G1,np.maximum(N_S,N_G2)))

        properties['nucl_pos'] += properties['k']*(properties['zposn'] - properties['nucl_pos'])*dt \
                                  + crowding_force(cells, a=properties['a'], s=properties['s'])*dt \
                                  + np.sqrt(2*properties['D']*dt)*np.random.randn(len(properties['zposn']))

        
        
        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+properties['nucl_pos']**2) # target area now depends on nucl_pos
        properties['A0_initial'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2) # initial
        
    
        #############BAETTI VID
        
        #ath BAETTI D VID
        cells.mesh , number_T1, d= cells.mesh.transition(T1_eps)  #check edges verifing T1 transition
        if len(number_T1)>0:
            for ii in number_T1:
                index = cells.mesh.face_id_by_edge[ii]
                properties['T1_angle_pD'] =np.append(properties['T1_angle_pD'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 

        #EG BAETTI VID
        len_modified=np.matrix.copy(cells.mesh.length)
        #len_modified[np.where((properties['parent_group'][cells.mesh.face_id_by_edge]==1) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==0))]*=0.12
        
        tens = (0.5*cells.by_edge('Lambda', 'Lambda_boundary')/len_modified)*cells.mesh.edge_vect
        tens= tens - tens.take(cells.mesh.edges.prev, 1)
        F+=tens

        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        properties['force_x'] = F[0]*viscosity
        properties['force_y'] = F[1]*viscosity   
        
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        yield cells       

# simulation with division and INM (no differentiation rate domain) # type 0 Rebeca model 4 -> type 7: drift = v0 & noise & force
def simulation_with_division_model_4(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None): #(cells,force,dt=0.001,T1_eps=0.04,lifespan=100.0,rand=Non
    properties = cells.properties
    properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    properties['ageingrate'] = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells)) #degradation rate per each cell
    properties['ids_division'] = [] #save ids of the cell os the division when it's ready per each time step
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pD'] = []

    properties['nucl_pos'] = properties['zposn'].copy()
    
    properties['force_z'] = []
    
    expansion = np.array([0.0,0.0])
   
    while True:
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        ready = np.where(~cells.empty() & (cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)) & (cells.properties['nucl_pos']>=0.75))[0]  # divides if nucleus pos > 0.75
        if len(ready): #these are the cells ready to undergo division at the current timestep
            properties['ageingrate'] =np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2*len(ready))))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  # Daugthers and parent have the same ids
            
            # daughter cells = same positions as parents - do this NEXT
            properties['nucl_pos'] = np.append(properties['nucl_pos'], np.ones(2*len(ready)))

            properties['ids_division'] = ready
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready] #New edges after division 
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                properties['Division_angle_pD']= np.append(properties['Division_angle_pD'],cells.mesh.edge_angle[division_new_edge][0])
        #properties['age'] = properties['age']+dt*properties['ageingrate'] #add time step depending of the#properties['age_ready'] = properties['age'][ready] # get age of cells which are READY to divide

        # IMPORTANT: add age ingrate only for alive cells
        alivecells = np.where(~cells.empty())[0]
        properties['age'][alivecells] += dt*properties['ageingrate'][alivecells] #add time step depending of the degradation rate 
        
        #Calculate z nuclei position (Apical-Basal movement), depending of the cell cycle phase time and age of the cell
        N_G1=1-1.0/t_G1*properties['age']
        N_S=0 
        N_G2=1.0/(t_G2)*(properties['age']-(t_G1+t_S))
        properties['zposn'] = np.minimum(1.0,np.maximum(N_G1,np.maximum(N_S,N_G2)))

        v0 = np.zeros_like(properties['age']) # array for velocities
        G1_cells = np.where(0 <= properties['age'] <= t_G1)[0]
        v0[G1_cells] = -1.0/t_G1
        #S_cells = np.where(t_G1 < properties['age'] <= t_G1 + t_S)[0]
        #v0[S_cells] = 0
        G2_cells = np.where(t_G1+t_S < properties['age'] <= t_G1+t_S+t_G2)[0]
        v0[G2_cells] = 1.0/t_G2
        #M_cells = np.where(t_G1+t_S+t_G2 < properties['age'])
        #v0[M_cells] = 0

        properties['nucl_pos'] += v0*dt \
                                  + crowding_force(cells, a=properties['a'], s=properties['s'])*dt \
                                  + np.sqrt(2*properties['D']*dt)*np.random.randn(len(properties['zposn'])) 

        
        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+properties['nucl_pos']**2) # target area now depends on nucl_pos
        properties['A0_initial'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2) # initial
        
        #############BAETTI VID
        
        #ath BAETTI D VID
        cells.mesh , number_T1, d= cells.mesh.transition(T1_eps)  #check edges verifing T1 transition
        if len(number_T1)>0:
            for ii in number_T1:
                index = cells.mesh.face_id_by_edge[ii]
                properties['T1_angle_pD'] =np.append(properties['T1_angle_pD'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 

        #EG BAETTI VID
        len_modified=np.matrix.copy(cells.mesh.length)
        #len_modified[np.where((properties['parent_group'][cells.mesh.face_id_by_edge]==1) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==0))]*=0.12
        
        tens = (0.5*cells.by_edge('Lambda', 'Lambda_boundary')/len_modified)*cells.mesh.edge_vect
        tens= tens - tens.take(cells.mesh.edges.prev, 1)
        F+=tens

        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        properties['force_x'] = F[0]*viscosity
        properties['force_y'] = F[1]*viscosity   
        
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        yield cells       

# simulation with division, INM and no differentiation rate - type 0 with both pD and pMN
def simulation_with_division_clone(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None): #(cells,force,dt=0.001,T1_eps=0.04,lifespan=100.0,rand=Non
    properties = cells.properties
    properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    properties['ageingrate'] = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells)) #degradation rate per each cell
    properties['ids_division'] = [] #save ids of the cell os the division when it's ready per each time step
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pMN'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pMN'] = []
    properties['Division_angle_pD'] = []
    properties['ed']=[]
    properties['T1_edge']=[]
    properties['ids_removed']=[]

    expansion = np.array([0.0,0.0])
    while True:
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        ready = np.where(~cells.empty() & (cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)))[0]  
        if len(ready): #these are the cells ready to undergo division at the current timestep
            properties['ageingrate'] =np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2.0*len(ready))))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  # Daugthers and parent have the same ids
            properties['ids_division'] = ready
            properties['parent_group'] = np.append(properties['parent_group'],np.repeat(properties['parent_group'][ready],2)) #use to draw clones
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready] #New edges after division
            properties['ed'].append(edge_pairs)
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                properties['la']=[commun_edges, division_new_edge]
                if properties['parent_group'][cells.mesh.n_face-2*(i+1)]==1:
                    properties['Division_angle_pMN']= np.append(properties['Division_angle_pMN'],cells.mesh.edge_angle[division_new_edge][0])
                else:
                    properties['Division_angle_pD']= np.append(properties['Division_angle_pD'],cells.mesh.edge_angle[division_new_edge][0])
        properties['age'] = properties['age']+dt*properties['ageingrate'] #add time step depending of the degradation rate 
        
        """Calculate z nuclei position (Apical-Basal movement), depending of the cell cycle phase time and age of the cell"""
        N_G1=1-1.0/t_G1*properties['age'] #nuclei position in G1 phase
        N_S=0
        N_G2=1.0/(t_G2)*(properties['age']-(t_G1+t_S))  #nuclei position in G2 and S phase
        properties['zposn'] = np.minimum(1.0,np.maximum(N_G1,np.maximum(N_S,N_G2)))
        
        
        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2)
        

        cells.mesh , number_T1, ids_removed= cells.mesh.transition(T1_eps)
        if len(ids_removed)>0:
            properties['ids_removed'].append(ids_removed)

        if len(number_T1)>0:
            for ii in number_T1:
                properties['T1_edge']=np.append(properties['T1_edge'], ii)
                index = cells.mesh.face_id_by_edge[ii]
                if properties['parent_group'][index]==1:
                    properties['T1_angle_pMN'] =np.append(properties['T1_angle_pMN'],cells.mesh.edge_angle[ii])
                else:
                    properties['T1_angle_pD'] =np.append(properties['T1_angle_pD'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 
        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        properties['force_x'] = F[0]*viscosity
        properties['force_y'] = F[1]*viscosity   
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        yield cells 


# simulation with division and all the cells are differentated as pNM diff_rate_hours (1/h) - type 1
def simulation_with_division_clone_differentiation(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None): #(cells,force,dt=0.001,T1_eps=0.04,lifespan=100.0,rand=Non
    properties = cells.properties
    properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    properties['ageingrate'] = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells)) #degradation rate per each cell
    properties['ids_division'] = [] #save ids of the cell os the division when it's ready per each time step
    properties['poisoned'] = np.zeros(len(cells)) ### to add diferenciation rate in PMN

    expansion = np.array([0.0,0.0])
    while True:
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        ready = np.where(~cells.empty() & (cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)))[0]  
        if len(ready): #these are the cells ready to undergo division at the current timestep
            properties['ageingrate'] =np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2.0*len(ready))))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  # Daugthers and parent have the same ids
            properties['ids_division'] = ready
            properties['parent_group'] = np.append(properties['parent_group'],np.repeat(properties['parent_group'][ready],2)) #use to draw clones
            properties['poisoned'] = np.append(properties['poisoned'], np.zeros(2*len(ready))) ### to add diferenciation rate in PMN
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready] #New edges after division 
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
        ###### Defferenciation rate
        properties['poisoned'] = properties['poisoned'] - (properties['poisoned']-1) * (~(cells.empty()) & (rand.rand(len(cells)) < (dt*diff_rate_hours*time_hours)))
        properties['age'] = properties['age']+dt*properties['ageingrate'] #add time step depending of the degradation rate 
        
        """Calculate z nuclei position (Apical-Basal movement), depending of the cell cycle phase time and age of the cell"""
        N_G1=1-1.0/t_G1*properties['age'] #nuclei position in G1 phase
        N_S=0
        N_G2=1.0/(t_G2)*(properties['age']-(t_G1+t_S))  #nuclei position in G2 and S phase
        properties['zposn'] = np.minimum(1.0,np.maximum(N_G1,np.maximum(N_S,N_G2)))
        
        
        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2)*(1.0-cells.properties['poisoned'])
        
        cells.mesh , number_T1= cells.mesh.transition(T1_eps)  #check edges verifing T1 transition
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 
        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        yield cells 


# simulation with division with INM and 2 diferent populations (with and without differentiation rate) - type 2
def simulation_with_division_clone_differenciation_3stripes(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None): #(cells,force,dt=0.001,T1_eps=0.04,lifespan=100.0,rand=Non

    
    properties = cells.properties
    properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    properties['ageingrate'] = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells)) #degradation rate per each cell
    #floor plate cells
    no=len(properties['ageingrate'][np.where(properties['parent_group']==3)])
    properties['ageingrate'][np.where(properties['parent_group']==3)]=np.random.normal(0.5/lifespan,0.1/lifespan,no)

    properties['ids_division'] = [] #save ids of the cell os the division when it's ready per each time step
    properties['ids_division_1'] = [] #save ids of the cell os the division when it's ready per each time step and region
    properties['ids_division_02'] = [] #save ids of the cell os the division when it's ready per each time step and region
    properties['poisoned'] = np.zeros(len(cells)) ### to add diferenciation rate in PMN
    # properties['differentiation_rate']= np.zeros(len(cells),dtype=int)
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pMN'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pMN'] = []
    properties['Division_angle_pD'] = []
    properties['deleted_edges']=[]
    properties['edges_division']=[]
    properties['T1_edge']=[]
    properties['nucl_pos'] = []
    expansion = np.array([0.0,0.0])
    
    # added by Rebeca:
    D = 0.15
    k = 100 # a constant
    a=0.2 # a constant controlling the size of force
    s=0.2
    
    
    while True:
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        
        ready = np.where((~cells.empty()  &(cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)) & (cells.properties['zposn']>=0.75)) | (~cells.empty() & (cells.properties['parent_group']==3) &(cells.mesh.area>=0.4) & (cells.properties['age']>=(t_G1+t_S+t_G2)) & (cells.properties['zposn']>=0.75)))[0] 

        if len(ready): #these are the cells ready to undergo division at the current timestep
            #properties['ageingrate'] =np.append(properties['ageingrate'], np.repeat(properties['ageingrate'][ready],2))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  # Daugthers and parent have the same ids
            properties['ids_division'] = ready
            properties['parent_group'] = np.append(properties['parent_group'],np.repeat(properties['parent_group'][ready],2)) #use to draw clones
            
            properties['poisoned'] = np.append(properties['poisoned'], np.zeros(2*len(ready))) ### to add diferenciation rate in PMN
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready]  #New edges after division 
            properties['edges_division'].append(edge_pairs)
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh  
            

           
            
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                
                
                if properties['parent_group'][ready[i]]==3:
                    properties['ageingrate'] = np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2)))
                else:
                    properties['ageingrate'] = np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2)))
                   # print len(properties['ageingrate'])
                #if properties['parent_group'][cells.mesh.n_face-2*(i+1)]==1:
                 #   properties['Division_angle_pMN']= np.append(properties['Division_angle_pMN'],cells.mesh.edge_angle[division_new_edge][0])
                #else:
                 #   properties['Division_angle_pD']= np.append(properties['Division_angle_pD'],cells.mesh.edge_angle[division_new_edge][0])
            #for ids in ready:
             #   if properties['parent_group'][ids]==1:
              #      properties['ids_division_1'] = np.append(properties['ids_division_1'], ids)
               # else:
                #    properties['ids_division_02'] = np.append(properties['ids_division_02'], ids)
        ###### Differentiation rate
        properties['differentiation_rate'] = time_hours*dt*(np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))[properties['parent_group']] #Used 0.02, 0.0002 & 1/13
        properties['poisoned'] = properties['poisoned'] - (properties['poisoned']-1) * (~(cells.empty()) & (rand.rand(len(cells)) < properties['differentiation_rate']))
        
        # add age ingrate only for alive cells
        alivecells = np.where(~cells.empty())[0]
        properties['age'][alivecells] += dt*properties['ageingrate'][alivecells] #add time step depending of the degradation rate 
        
        
        N_G1=1-1.0/t_G1*properties['age'] #nuclei position in G1 phase
        N_S=0
        N_G2=1.0/(t_G2)*(properties['age']-(t_G1+t_S))  #nuclei position in G2 and S phase
        
        properties['zposn'] = np.minimum(1.0,np.maximum(N_G1,np.maximum(N_S,N_G2)))
       
        
        properties['zposn'][np.where(properties['parent_group']==3)]=0
       

        
        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2)*(1.0-cells.properties['poisoned'])

        
        cells.mesh , number_T1, del_edges= cells.mesh.transition(T1_eps)
        
        properties['deleted_edges'].append(del_edges)
        #if len(number_T1)>0:
         #   for ii in number_T1:
          #      properties['T1_edge']=np.append(properties['T1_edge'], ii)
               # index = cells.mesh.face_id_by_edge[ii]
                #if properties['parent_group'][index]==1:
                 #   properties['T1_angle_pMN'] =np.append(properties['T1_angle_pMN'],cells.mesh.edge_angle[ii])
                #else:
                 #   properties['T1_angle_pD'] =np.append(properties['T1_angle_pD'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 

        #ADDING THE FLOOR PLATE
        len_modified=np.matrix.copy(cells.mesh.length)
        #this next line can be use for modifying all edges at the boundary
        #len_modified[np.where( ((properties['parent_group'][cells.mesh.face_id_by_edge]==3) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==3)) | ((properties['parent_group'][cells.mesh.face_id_by_edge]==2) & (properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==3)) | ((properties['parent_group'][cells.mesh.face_id_by_edge]==4) & (properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==3)) )]*=0.004
        
        #modifying tension w.r.t. area
        for n in np.where( ((properties['parent_group'][cells.mesh.face_id_by_edge]==3) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==3)))[0]:
            if abs(cells.mesh.area[cells.mesh.face_id_by_edge[n]]-cells.mesh.area[cells.mesh.face_id_by_edge[cells.mesh.edges.reverse[n]]])>0.4:
                len_modified[n]*=0.002
   
        tens = (0.5*cells.by_edge('Lambda', 'Lambda_boundary')/len_modified)*cells.mesh.edge_vect
        tens= tens - tens.take(cells.mesh.edges.prev, 1)
        F+=tens

        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        #properties['force_x'] = F[0]*viscosity
        #properties['force_y'] = F[1]*viscosity
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        
        yield cells 

# simulation with division with INM and all the cells high without differentiation rate - type 3 ??
def simulation_with_division_clone_whole_tissue_differenciation(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None): #(cells,force,dt=0.001,T1_eps=0.04,lifespan=100.0,rand=None
    properties = cells.properties
    properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    properties['ageingrate'] = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells)) #degradation rate per each cell
    properties['ids_division'] = [] #save ids of the cell os the division when it's ready per each time step
    properties['ids_division_1'] = [] #save ids of the cell os the division when it's ready per each time step and region
    properties['ids_division_02'] = [] #save ids of the cell os the division when it's ready per each time step and region
    properties['poisoned'] = np.zeros(len(cells)) ### to add diferenciation rate in PMN
    # properties['differentiation_rate']= np.zeros(len(cells),dtype=int)
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pMN'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pMN'] = []
    properties['Division_angle_pD'] = []
    properties['ids_removed']=[]
    expansion = np.array([0.0,0.0])
    while True:
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        ready = np.where(~cells.empty() & (cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)))[0]  
        if len(ready): #these are the cells ready to undergo division at the current timestep
            properties['ageingrate'] =np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2.0*len(ready)))) # use only for alive cells at current timestep
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  # Daugthers and parent have the same ids
            properties['ids_division'] = ready
            properties['parent_group'] = np.append(properties['parent_group'],np.repeat(properties['parent_group'][ready],2)) #use to draw clones
            properties['poisoned'] = np.append(properties['poisoned'], np.zeros(2*len(ready))) ### to add diferenciation rate in PMN
            # properties['differentiation_rate'] = np.append(properties['differentiation_rate'], diff_rate_hours*np.ones(2*len(ready)))          
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready] #New edges after division 
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                if properties['parent_group'][cells.mesh.n_face-2*(i+1)]==1:
                    properties['Division_angle_pMN']= np.append(properties['Division_angle_pMN'],cells.mesh.edge_angle[division_new_edge][0])
                else: #antes estaba mal y ponia pD!!!! las simulaciones guardadas estan mal
                    properties['Division_angle_pMN']= np.append(properties['Division_angle_pMN'],cells.mesh.edge_angle[division_new_edge][0])
            for ids in ready:
                if properties['parent_group'][ids]==1:
                    properties['ids_division_1'] = np.append(properties['ids_division_1'], ids)
                else:#antes estaba mal y ponia 0_2!!!! las simulaciones guardadas estan mal
                    properties['ids_division_1'] = np.append(properties['ids_division_1'], ids)
        ###### Defferentiation rate
        properties['differentiation_rate'] = time_hours*dt*(np.array([0.0,diff_rate_hours,0.0,0.0,0.0,diff_rate_hours, 0.0]))[properties['parent_group']] #Used 0.02, 0.0002 & 1/13
        properties['poisoned'] = properties['poisoned'] - (properties['poisoned']-1) * (~(cells.empty()) & (rand.rand(len(cells)) < properties['differentiation_rate']))
        
        # add age ingrate only for alive cells
        alivecells = np.where(~cells.empty())[0]
        properties['age'][alivecells] += dt*properties['ageingrate'][alivecells] #add time step depending of the degradation rate 
        
        N_G1=1-1.0/t_G1*properties['age'] #nuclei position in G1 phase
        N_S=0
        N_G2=1.0/(t_G2)*(properties['age']-(t_G1+t_S))  #nuclei position in G2 and S phase
        properties['zposn'] = np.minimum(1.0,np.maximum(N_G1,np.maximum(N_S,N_G2)))
        
        
        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2)*(1.0-cells.properties['poisoned'])
        
        cells.mesh , number_T1, ids_removed= cells.mesh.transition(T1_eps)

        properties['ids_removed'].append(ids_removed) #####BAETTI VID
        if len(number_T1)>0:
            for ii in number_T1:
                index = cells.mesh.face_id_by_edge[ii]
                if properties['parent_group'][index]==1:
                    properties['T1_angle_pMN'] =np.append(properties['T1_angle_pMN'],cells.mesh.edge_angle[ii])
                else:#antes estaba mal y ponia pD!!!! las simulaciones guardadas estan mal
                    properties['T1_angle_pMN'] =np.append(properties['T1_angle_pMN'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 
        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        properties['force_x'] = F[0]*viscosity
        properties['force_y'] = F[1]*viscosity
    
        
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        yield cells 
        
    
# simulation with division with INM and 2 diferent populations (with and without differentiation rate) - type 4 PREVIOUS (type 2 modified by Rebeca)
def simulation_with_division_clone_differenciation_3stripes_model_1(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None): #(cells,force,dt=0.001,T1_eps=0.04,lifespan=100.0,rand=Non

    
    properties = cells.properties
    properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    properties['ageingrate'] = np.random.normal(1.0/lifespan,0.2/lifespan,len(cells)) #degradation rate per each cell
    #floor plate cells
    no=len(properties['ageingrate'][np.where(properties['parent_group']==3)])
    properties['ageingrate'][np.where(properties['parent_group']==3)]=np.random.normal(0.5/lifespan,0.1/lifespan,no)

    properties['ids_division'] = [] #save ids of the cell os the division when it's ready per each time step
    properties['ids_division_1'] = [] #save ids of the cell os the division when it's ready per each time step and region
    properties['ids_division_02'] = [] #save ids of the cell os the division when it's ready per each time step and region
    properties['poisoned'] = np.zeros(len(cells)) ### to add diferenciation rate in PMN
    # properties['differentiation_rate']= np.zeros(len(cells),dtype=int)
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pMN'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pMN'] = []
    properties['Division_angle_pD'] = []
    properties['deleted_edges']=[]
    properties['edges_division']=[]
    properties['T1_edge']=[]
    properties['nucl_pos']=[]
    expansion = np.array([0.0,0.0])
    
    # added by Rebeca:
    D = 0.15
    k = 100 # a constant
    a=0.2 # a constant controlling the size of force
    s=0.2
    

    
    while True:
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        
        # added by Rebeca: properties['zposn'] >= 0.75
        ready = np.where((~cells.empty()  &(cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)) & (cells.properties['nucl_pos']>=0.75)) | (~cells.empty() & (cells.properties['parent_group']==3) &(cells.mesh.area>=0.4) & (cells.properties['age']>=(t_G1+t_S+t_G2)) & (cells.properties['nucl_pos']>=0.75)))[0]  

        if len(ready): #these are the cells ready to undergo division at the current timestep
            #properties['ageingrate'] =np.append(properties['ageingrate'], np.repeat(properties['ageingrate'][ready],2))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  # Daugthers and parent have the same ids
            properties['ids_division'] = ready
            properties['parent_group'] = np.append(properties['parent_group'],np.repeat(properties['parent_group'][ready],2)) #use to draw clones
            
            properties['poisoned'] = np.append(properties['poisoned'], np.zeros(2*len(ready))) ### to add diferenciation rate in PMN
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready]  #New edges after division 
            properties['edges_division'].append(edge_pairs)
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
            
       
            
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                
                

              
                if properties['parent_group'][ready[i]]==3:
                    properties['ageingrate'] = np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2)))
                else:
                    properties['ageingrate'] = np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2)))
                   # print len(properties['ageingrate'])
                #if properties['parent_group'][cells.mesh.n_face-2*(i+1)]==1:
                 #   properties['Division_angle_pMN']= np.append(properties['Division_angle_pMN'],cells.mesh.edge_angle[division_new_edge][0])
                #else:
                 #   properties['Division_angle_pD']= np.append(properties['Division_angle_pD'],cells.mesh.edge_angle[division_new_edge][0])
            #for ids in ready:
             #   if properties['parent_group'][ids]==1:
              #      properties['ids_division_1'] = np.append(properties['ids_division_1'], ids)
               # else:
                #    properties['ids_division_02'] = np.append(properties['ids_division_02'], ids)
        ###### Differentiation rate
        properties['differentiation_rate'] = time_hours*dt*(np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))[properties['parent_group']] #Used 0.02, 0.0002 & 1/13
        properties['poisoned'] = properties['poisoned'] - (properties['poisoned']-1) * (~(cells.empty()) & (rand.rand(len(cells)) < properties['differentiation_rate']))
        
        # add age ingrate only for alive cells -> not necessarily needed as all cells are alive in type 2
        alivecells = np.where(~cells.empty())[0]
        properties['age'][alivecells] += dt*properties['ageingrate'][alivecells] #add time step depending of the degradation rate 
        
        
        N_G1_new = np.ones_like(properties['age']) # vector with ones for as many cells (ages) as they are at the current timestep
        N_S_new = np.zeros_like(properties['age'])
        N_G2_new = np.zeros_like(properties['age'])
        N_M_new = np.ones_like(properties['age'])
        
        N_G1=1-1.0/t_G1*properties['age']
        N_G1_new = N_G1_new + k*(N_G1 - N_G1_new)*dt + np.sqrt(2*D*dt)*np.random.randn(len(properties['age'])) #nuclei position in G1 phase
       
        N_S=0 
        N_S_new = N_S_new + k*(N_S - N_S_new)*dt + np.sqrt(2*D*dt)*np.random.randn(len(properties['age']))
        
        N_G2=1.0/(t_G2)*(properties['age']-(t_G1+t_S))  #nuclei position in G2 and S phase
        N_G2_new = N_G2_new + k*(N_G2 - N_G2_new)*dt + np.sqrt(2*D*dt)*np.random.randn(len(properties['age']))
        
        N_M = 1
        N_M_new = N_M_new + k*(N_M - N_M_new)*dt + np.sqrt(2*D*dt)*np.random.randn(len(properties['age']))
        properties['zposn'] = np.minimum(1.0,np.maximum(N_G1,np.maximum(N_S,N_G2)))
        properties['nucl_pos'] = np.minimum(1.0,np.maximum(N_G1_new,np.maximum(N_S_new,N_G2_new)))
        
        
        properties['zposn'][np.where(properties['parent_group']==3)]=0

        
        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2)*(1.0-cells.properties['poisoned'])

        
        cells.mesh , number_T1, del_edges= cells.mesh.transition(T1_eps)
        
        properties['deleted_edges'].append(del_edges)
        #if len(number_T1)>0:
         #   for ii in number_T1:
          #      properties['T1_edge']=np.append(properties['T1_edge'], ii)
               # index = cells.mesh.face_id_by_edge[ii]
                #if properties['parent_group'][index]==1:
                 #   properties['T1_angle_pMN'] =np.append(properties['T1_angle_pMN'],cells.mesh.edge_angle[ii])
                #else:
                 #   properties['T1_angle_pD'] =np.append(properties['T1_angle_pD'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 

        #ADDING THE FLOOR PLATE
        len_modified=np.matrix.copy(cells.mesh.length)
        #this next line can be use for modifying all edges at the boundary
        #len_modified[np.where( ((properties['parent_group'][cells.mesh.face_id_by_edge]==3) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==3)) | ((properties['parent_group'][cells.mesh.face_id_by_edge]==2) & (properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==3)) | ((properties['parent_group'][cells.mesh.face_id_by_edge]==4) & (properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==3)) )]*=0.004
        
        #modifying tension w.r.t. area
        for n in np.where( ((properties['parent_group'][cells.mesh.face_id_by_edge]==3) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==3)))[0]:
            if abs(cells.mesh.area[cells.mesh.face_id_by_edge[n]]-cells.mesh.area[cells.mesh.face_id_by_edge[cells.mesh.edges.reverse[n]]])>0.4:
                len_modified[n]*=0.002
   
        tens = (0.5*cells.by_edge('Lambda', 'Lambda_boundary')/len_modified)*cells.mesh.edge_vect
        tens= tens - tens.take(cells.mesh.edges.prev, 1)
        F+=tens

        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        #properties['force_x'] = F[0]*viscosity
        #properties['force_y'] = F[1]*viscosity
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        
        yield cells         
#simulation with division & INM with new formulation by Aaryan tldr: is that k now has time dependence, this is to allow different forces for g1 and g2 without having cells move to target positions blindly eg a cell that did not make it to the basal end by S phase would drop immediately to the basal end at the start of G2 where downard movement has not been described formulation below:
'''
dz/dt=K(t)(z_0 (t)−z_t )dt+√((2D)normal) (t)

z_(t+1)=z_t+K(t)(z_0 (t)−z_t )dt

-	G1	S	G2	M
K	12	0	24	24
z_0  	0	-	1	1
Expected duration	0.4	0.33	0.16	0.11


Active downard force by cytoskeleton^
'''
def simulation_with_division_model_5(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None):
    lifespan=46800/460
    random.seed(1999)
    properties=cells.properties
    properties['parent']=cells.mesh.face_ids
    #save the ids to control division parents-daughters
    properties['ageingrate']=np.random.normal(1.0/lifespan,0.2/lifespan,len(cells))
    #degradation rate per each cell which varies by cell
    properties['ids_division']=[] #save ids of the cell os the division when its ready per each time step
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pD'] = []
    
    properties['k(t)']=[0]*len(properties['zposn'])
    #contains the current k term (inertia of INM) as a function of time, it is selected from the k array fed to this type of simulation

    
    
    properties['force_z'] = []
    
    expansion = np.array([0.0,0.0])
    iteration_tracker=0
    while True:
        if iteration_tracker%1000==0:
            print(f"at timepoint {iteration_tracker/1000}")        
        iteration_tracker+=1
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        ready = np.where(~cells.empty() & (cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)) & (cells.properties['nucl_pos']>=0.75))[0]  # divides if nucleus pos > 0.75
        if len(ready): #these are the cells ready to undergo division at the current timestep
            properties['ageingrate'] =np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2*len(ready))))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['k(t)']=np.append(properties['k(t)'],np.zeros(2*len(ready)))
            #appending 0 because it is determined outside of current if function
            properties['zposn']=np.append(properties['zposn'],np.zeros(2*len(ready))) #appending 0 because it is determined outside of if function 
            
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  
            
            # daughter cells = same positions as parents - do this NEXT
            properties['nucl_pos'] = np.append(properties['nucl_pos'], np.ones(2*len(ready)))

            properties['ids_division'] = ready
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready] #New edges after division 
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                properties['Division_angle_pD']= np.append(properties['Division_angle_pD'],cells.mesh.edge_angle[division_new_edge][0])
        #properties['age'] = properties['age']+dt*properties['ageingrate'] #add time step depending of the#properties['age_ready'] = properties['age'][ready] # get age of cells which are READY to divide

        # IMPORTANT: add age ingrate only for alive cells
        alivecells = np.where(~cells.empty())[0]
        properties['age'][alivecells] += dt*properties['ageingrate'][alivecells] #add time step depending of the degradation rate 
        
        #calculating z nuclei position now target position is either 0 or 1 and K has a time dependence depending on age 
        N_G1=0
        N_S=0
        N_G2=1
        N_M=1
        #need to scale k between biological time and simulation time ; the scale is lifepan one cell divsion is 100 seconds
        k_G1=properties['k'][0]/lifespan
        k_S=properties['k'][1]/lifespan
        k_G2=properties['k'][2]/lifespan
        k_M=properties['k'][3]/lifespan
        #now updating cells.properties['zposn'] by their age
        #nucl_pos is now zposn and zposn is now either 0 or 1
        n_of_cells=len(cells.properties['zposn'])
        cells.properties['k(t)']=np.zeros(n_of_cells)
        G1_index=np.where(np.logical_and(0<=cells.properties['age'],cells.properties['age']<t_G1))
        S_index=np.where(np.logical_and(t_G1<=cells.properties['age'],cells.properties['age']<t_G1+t_S))
        G2_index=np.where(np.logical_and(t_G1+t_S<=cells.properties['age'],cells.properties['age']<t_G1+t_S+t_G2))
        M_index=np.where(cells.properties['age']>=t_G1+t_S+t_G2)
        cells.properties['zposn'][G1_index]=N_G1
        cells.properties['k(t)'][G1_index]=k_G1
        cells.properties['zposn'][S_index]=N_S
        cells.properties['k(t)'][S_index]=k_S
        cells.properties['zposn'][G2_index]=N_G2
        cells.properties['k(t)'][G2_index]=k_G2
        cells.properties['zposn'][M_index]=N_M
        cells.properties['k(t)'][M_index]=k_M
        
        # if iteration_tracker%1000==0:
        #     plt.hist2d(properties['age'],properties['nucl_pos'],norm=mpl.colors.LogNorm(),bins=500)
        #     plt.show()
        properties['nucl_pos']=properties['nucl_pos']+properties['k(t)']*(properties['zposn']-properties['nucl_pos'])*dt + crowding_force(cells,a=properties['a'],s=properties['s'])*dt + np.sqrt(2*properties['D']*dt)*np.random.randn(len(properties['zposn']))
        
        max_nucl=np.array([min(i,1)for i in properties['nucl_pos']])
        #nucl_pos cannot exceed 1
        properties['nucl_pos']=max_nucl
        #nucl_pos cannot be less than 0 
        max_nucl=np.array([max(i,0)for i in properties['nucl_pos']])
        properties['nucl_pos']=max_nucl
        """Target area function depending age and z nuclei position"""
        properties['A0'] = (properties['age']+1.0)*0.5*(1.0+max_nucl**2) # target area now depends on nucl_pos
        #properties['A0_initial'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2) # initial
        #CHANGE BY AARYAN TO MAKE it so that
        '''
        instead of having A0 determined by (1+nucl_pos^2) which has a shape like 
       ┌────────────────────────────────────────────────┐
A0 104 │⠀⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠀│ 
       │⠀⠘⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠇⠀│
       │⠀⠀⢱⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡜⠀⠀│
       │⠀⠀⠀⢇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡸⠀⠀⠀│
       │⠀⠀⠀⠘⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠃⠀⠀⠀│
       │⠀⠀⠀⠀⠸⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠇⠀⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠱⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠎⠀⠀⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠣⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠎⠀⠀⠀⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠱⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠎⠀⠀⠀⠀⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠀⠱⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠎⠀⠀⠀⠀⠀⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⢀⠎⠀⠀⠀⠀⠀⠀⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⢄⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⡠⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠣⡀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⢀⠜⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠒⢄⠀⠀⠀⠀⡇⠀⠀⠀⡠⠒⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀ │
    -2 │⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠭⠵⠶⠤⡧⠴⠮⠭⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤│
       └────────────────────────────────────────┘
       ⠀-10.6⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀nucl_pos               10.6⠀

       we want the following shape from (1+max(0,nucl_pos)^2))
       ┌────────────────────────────────────────────────┐
A0 104 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠀│ 
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠇⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡜⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡸⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠃⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠇⠀⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠎⠀⠀⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠎⠀⠀⠀⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠎⠀⠀⠀⠀⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠎⠀⠀⠀⠀⠀⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⢀⠎⠀⠀⠀⠀⠀⠀⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⡠⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⢀⠜⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
       │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⡠⠒⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
    -2 │⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⡧⠴⠮⠭⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤│
       └────────────────────────────────────────┘
       ⠀-10.6⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀Nucl_Pos⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀10.6⠀
        '''

         #############BAETTI VID
        
        #ath BAETTI D VID
        cells.mesh , number_T1, d= cells.mesh.transition(T1_eps)  #check edges verifing T1 transition
        #UNKNOWN BUG HERE NEEDS FIXING BY AARYAN
        #if len(number_T1)>0:
         #   for ii in number_T1:
          #      index = cells.mesh.face_id_by_edge[ii]
           #     properties['T1_angle_pD'] #=np.append(properties['T1_angle_pD'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 

        #EG BAETTI VID
        len_modified=np.matrix.copy(cells.mesh.length)
        #len_modified[np.where((properties['parent_group'][cells.mesh.face_id_by_edge]==1) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==0))]*=0.12
        
        tens = (0.5*cells.by_edge('Lambda', 'Lambda_boundary')/len_modified)*cells.mesh.edge_vect
        tens= tens - tens.take(cells.mesh.edges.prev, 1)
        F+=tens

        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        properties['force_x'] = F[0]*viscosity
        properties['force_y'] = F[1]*viscosity   
        
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        
        yield cells       



    #simulation with division & INM with new formulation & new age equation to limit increase of target area by Aaryan tldr: same as simulation_with_division_model_5 except we are now changing the equation: for target area:

'''
BEFORE:
A_0=((age+1)(1+z^2))/2
         ┌────────────────────────────────────────┐
   11.33 │⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⠄⠀│ y1
         │⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡤⠊⠀⠀⠀│
         │⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⠔⠉⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠤⠊⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠒⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠔⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠔⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡠⠊⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠒⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⣀⠔⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⢀⠤⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⢸⠀⠀⠀⡠⠒⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⢸⢀⠤⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⡠⢺⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
   -0.33 │⠤⠴⠭⠤⢼⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤│
         └────────────────────────────────────────┘
         ⠀-1.33⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀10.33⠀

AFTER:
A_0=max*(1+z^2)/(1+exp(-r*age))
where max is the maximum value of the logistic function
and r is the rate of growth of said function
r is determined to be = -np.log((1+0.75**2)/(A_c)-1)/(1-t_M)
this is chosen such that eqn 2 has coordinates ((1-t_M),A_c) when z=0.75
the rationale is that cells should have meet their division requirements at roughly the same point.
            ┌────────────────────────────────────────┐
    1.02188 │⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⡠⠤⠤⠤⠤⠤⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠂⠀│ y1
            │⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⡠⠔⠊⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
            │⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⡤⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
            │⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⡜⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
            │⠀⠀⠀⠀⢸⠀⠀⠀⢀⠎⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
            │⠀⠀⠀⠀⢸⠀⠀⢠⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
            │⠀⠀⠀⠀⢸⠀⠀⡎⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
            │⠀⠀⠀⠀⢸⠀⡜⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
            │⠀⠀⠀⠀⢸⢰⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
            │⠀⠀⠀⠀⢸⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
            │⠀⠀⠀⠀⣾⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
            │⠀⠀⠀⡜⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
            │⠀⠀⡰⠁⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
            │⠀⢠⠃⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
   0.247011 │⠀⠎⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
            └────────────────────────────────────────┘
            ⠀-1.33⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀10.33⠀

dz/dt=K(t)(z_0 (t)−z_t )dt+√((2D)normal) (t)

z_(t+1)=z_t+K(t)(z_0 (t)−z_t )dt

-	G1	S	G2	M
K	12	0	24	24
z_0  	0	-	1	1
Expected duration	0.4	0.33	0.16	0.11


Active downard force by cytoskeleton^
'''
def simulation_with_division_model_6(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None):
    lifespan=46800/460
    random.seed(1999)
    properties=cells.properties
    properties['parent']=cells.mesh.face_ids
    #save the ids to control division parents-daughters
    properties['ageingrate']=np.random.normal(1.0/lifespan,0.2/lifespan,len(cells))
    #degradation rate per each cell which varies by cell
    properties['ids_division']=[] #save ids of the cell os the division when its ready per each time step
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pD'] = []
    
    properties['k(t)']=[0]*len(properties['zposn'])
    #contains the current k term (inertia of INM) as a function of time, it is selected from the k array fed to this type of simulation

    
    
    properties['force_z'] = []
    
    expansion = np.array([0.0,0.0])
    iteration_tracker=0
    while True:
        if iteration_tracker%1000==0:
            print(f"at timepoint {iteration_tracker/1000}")        
        iteration_tracker+=1
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        ready = np.where(~cells.empty() & (cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)) & (cells.properties['nucl_pos']>=0.75))[0]  # divides if nucleus pos > 0.75
        if len(ready): #these are the cells ready to undergo division at the current timestep
            properties['ageingrate'] =np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2*len(ready))))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['k(t)']=np.append(properties['k(t)'],np.zeros(2*len(ready)))
            #appending 0 because it is determined outside of current if function
            properties['zposn']=np.append(properties['zposn'],np.zeros(2*len(ready))) #appending 0 because it is determined outside of if function 
            
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  
            
            # daughter cells = same positions as parents - do this NEXT
            properties['nucl_pos'] = np.append(properties['nucl_pos'], np.ones(2*len(ready)))

            properties['ids_division'] = ready
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready] #New edges after division 
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                properties['Division_angle_pD']= np.append(properties['Division_angle_pD'],cells.mesh.edge_angle[division_new_edge][0])
        #properties['age'] = properties['age']+dt*properties['ageingrate'] #add time step depending of the#properties['age_ready'] = properties['age'][ready] # get age of cells which are READY to divide

        # IMPORTANT: add age ingrate only for alive cells
        alivecells = np.where(~cells.empty())[0]
        properties['age'][alivecells] += dt*properties['ageingrate'][alivecells] #add time step depending of the degradation rate 
        
        #calculating z nuclei position now target position is either 0 or 1 and K has a time dependence depending on age 
        N_G1=0
        N_S=0
        N_G2=1
        N_M=1
        #need to scale k between biological time and simulation time ; the scale is lifepan one cell divsion is 100 seconds
        k_G1=properties['k'][0]/lifespan
        k_S=properties['k'][1]/lifespan
        k_G2=properties['k'][2]/lifespan
        k_M=properties['k'][3]/lifespan
        #now updating cells.properties['zposn'] by their age
        #nucl_pos is now zposn and zposn is now either 0 or 1
        n_of_cells=len(cells.properties['zposn'])
        cells.properties['k(t)']=np.zeros(n_of_cells)
        G1_index=np.where(np.logical_and(0<=cells.properties['age'],cells.properties['age']<t_G1))
        S_index=np.where(np.logical_and(t_G1<=cells.properties['age'],cells.properties['age']<t_G1+t_S))
        G2_index=np.where(np.logical_and(t_G1+t_S<=cells.properties['age'],cells.properties['age']<t_G1+t_S+t_G2))
        M_index=np.where(cells.properties['age']>=t_G1+t_S+t_G2)
        cells.properties['zposn'][G1_index]=N_G1
        cells.properties['k(t)'][G1_index]=k_G1
        cells.properties['zposn'][S_index]=N_S
        cells.properties['k(t)'][S_index]=k_S
        cells.properties['zposn'][G2_index]=N_G2
        cells.properties['k(t)'][G2_index]=k_G2
        cells.properties['zposn'][M_index]=N_M
        cells.properties['k(t)'][M_index]=k_M
        
        # if iteration_tracker%1000==0:
        #     plt.hist2d(properties['age'],properties['nucl_pos'],norm=mpl.colors.LogNorm(),bins=500)
        #     plt.show()
        properties['nucl_pos']=properties['nucl_pos']+properties['k(t)']*(properties['zposn']-properties['nucl_pos'])*dt + crowding_force(cells,a=properties['a'],s=properties['s'])*dt + np.sqrt(2*properties['D']*dt)*np.random.randn(len(properties['zposn']))
        
        max_nucl=np.array([min(i,1)for i in properties['nucl_pos']])
        #nucl_pos cannot exceed 1
        properties['nucl_pos']=max_nucl
        #nucl_pos cannot be less than 0 
        max_nucl=np.array([max(i,0)for i in properties['nucl_pos']])
        properties['nucl_pos']=max_nucl
        """Target area function depending age and z nuclei position"""


        r=1.6
        max_value=2.5
        c=0.25
        properties['A0'] = calc_target_area(age=properties['age'],nucl_pos=max_nucl,r=r,max_value=max_value,c=c) # target area now depends on nucl_pos
        #properties['A0_initial'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2) 

         #############BAETTI VID
        
        #ath BAETTI D VID
        cells.mesh , number_T1, d= cells.mesh.transition(T1_eps)  #check edges verifing T1 transition
        #UNKNOWN BUG HERE NEEDS FIXING BY AARYAN
        #if len(number_T1)>0:
         #   for ii in number_T1:
          #      index = cells.mesh.face_id_by_edge[ii]
           #     properties['T1_angle_pD'] #=np.append(properties['T1_angle_pD'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 

        #EG BAETTI VID
        len_modified=np.matrix.copy(cells.mesh.length)
        #len_modified[np.where((properties['parent_group'][cells.mesh.face_id_by_edge]==1) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==0))]*=0.12
        
        tens = (0.5*cells.by_edge('Lambda', 'Lambda_boundary')/len_modified)*cells.mesh.edge_vect
        tens= tens - tens.take(cells.mesh.edges.prev, 1)
        F+=tens

        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        properties['force_x'] = F[0]*viscosity
        properties['force_y'] = F[1]*viscosity   
        
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        
        yield cells       

   #simulation with division & INM with new formulation & new age equation to limit increase of target area by Aaryan tldr: same as simulation_with_division_model_5 except we are now changing the equation: for target area:

'''
Same as model_6 except we are using new crowding force calculated via crowding_force_2()
'''
def simulation_with_division_model_7(cells,force,dt=dt,T1_eps=T1_eps,lifespan=100.0,rand=None):
    lifespan=46800/460
    random.seed(1999)
    properties=cells.properties
    properties['parent']=cells.mesh.face_ids
    #save the ids to control division parents-daughters
    properties['ageingrate']=np.random.normal(1.0/lifespan,0.2/lifespan,len(cells))
    #degradation rate per each cell which varies by cell
    properties['ids_division']=[] #save ids of the cell os the division when its ready per each time step
    properties['force_x'] = []
    properties['force_y'] = []
    properties['T1_angle_pD'] = []
    properties['Division_angle_pD'] = []
    
    properties['k(t)']=[0]*len(properties['zposn'])
    #contains the current k term (inertia of INM) as a function of time, it is selected from the k array fed to this type of simulation

    
    
    properties['force_z'] = []
    
    expansion = np.array([0.0,0.0])
    iteration_tracker=0
    while True:
        if iteration_tracker%10==0:
            print(f"at timepoint {iteration_tracker/1000}")        
        iteration_tracker+=1
        #cells id where is true the division conditions: living cells & area greater than 2 & age cell in mitosis 
        ready = np.where(~cells.empty() & (cells.mesh.area>=A_c) & (cells.properties['age']>=(t_G1+t_S+t_G2)) & (cells.properties['nucl_pos']>=0.75))[0]  # divides if nucleus pos > 0.75
        if len(ready): #these are the cells ready to undergo division at the current timestep
            properties['ageingrate'] =np.append(properties['ageingrate'], np.abs(np.random.normal(1.0/lifespan,0.2/lifespan,2*len(ready))))
            properties['age'] = np.append(properties['age'],np.zeros(2*len(ready)))
            properties['k(t)']=np.append(properties['k(t)'],np.zeros(2*len(ready)))
            #appending 0 because it is determined outside of current if function
            properties['zposn']=np.append(properties['zposn'],np.zeros(2*len(ready))) #appending 0 because it is determined outside of if function 
            
            properties['parent'] = np.append(properties['parent'],np.repeat(properties['parent'][ready],2))  
            
            # daughter cells = same positions as parents - do this NEXT
            properties['nucl_pos'] = np.append(properties['nucl_pos'], np.ones(2*len(ready)))

            properties['ids_division'] = ready
            edge_pairs = [division_axis(cells.mesh,cell_id,rand) for cell_id in ready] #New edges after division 
            cells.mesh = cells.mesh.add_edges(edge_pairs) #Add new edges in the mesh
            for i in range(len(ready)):
                commun_edges = np.intersect1d(cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-2*(i+1)))[0]],cells.mesh.length[np.where(cells.mesh.face_id_by_edge==(cells.mesh.n_face-1-2*i))[0]])
                division_new_edge=np.where(cells.mesh.length==np.max(commun_edges))[0]
                properties['Division_angle_pD']= np.append(properties['Division_angle_pD'],cells.mesh.edge_angle[division_new_edge][0])
        #properties['age'] = properties['age']+dt*properties['ageingrate'] #add time step depending of the#properties['age_ready'] = properties['age'][ready] # get age of cells which are READY to divide

        # IMPORTANT: add age ingrate only for alive cells
        alivecells = np.where(~cells.empty())[0]
        properties['age'][alivecells] += dt*properties['ageingrate'][alivecells] #add time step depending of the degradation rate 
        
        #calculating z nuclei position now target position is either 0 or 1 and K has a time dependence depending on age 
        N_G1=0
        N_S=0
        N_G2=1
        N_M=1
        #need to scale k between biological time and simulation time ; the scale is lifepan one cell divsion is 100 seconds
        k_G1=properties['k'][0]/lifespan
        k_S=properties['k'][1]/lifespan
        k_G2=properties['k'][2]/lifespan
        k_M=properties['k'][3]/lifespan
        #now updating cells.properties['zposn'] by their age
        #nucl_pos is now zposn and zposn is now either 0 or 1
        n_of_cells=len(cells.properties['zposn'])
        cells.properties['k(t)']=np.zeros(n_of_cells)
        G1_index=np.where(np.logical_and(0<=cells.properties['age'],cells.properties['age']<t_G1))
        S_index=np.where(np.logical_and(t_G1<=cells.properties['age'],cells.properties['age']<t_G1+t_S))
        G2_index=np.where(np.logical_and(t_G1+t_S<=cells.properties['age'],cells.properties['age']<t_G1+t_S+t_G2))
        M_index=np.where(cells.properties['age']>=t_G1+t_S+t_G2)
        cells.properties['zposn'][G1_index]=N_G1
        cells.properties['k(t)'][G1_index]=k_G1
        cells.properties['zposn'][S_index]=N_S
        cells.properties['k(t)'][S_index]=k_S
        cells.properties['zposn'][G2_index]=N_G2
        cells.properties['k(t)'][G2_index]=k_G2
        cells.properties['zposn'][M_index]=N_M
        cells.properties['k(t)'][M_index]=k_M
        
        # if iteration_tracker%1000==0:
        #     plt.hist2d(properties['age'],properties['nucl_pos'],norm=mpl.colors.LogNorm(),bins=500)
        #     plt.show()
        properties['nucl_pos']=properties['nucl_pos']+properties['k(t)']*(properties['zposn']-properties['nucl_pos'])*dt + crowding_force_2(cells,a=properties['a'],s=properties['s'])*dt + np.sqrt(2*properties['D']*dt)*np.random.randn(len(properties['zposn']))
        
        max_nucl=np.array([min(i,1)for i in properties['nucl_pos']])
        #nucl_pos cannot exceed 1
        properties['nucl_pos']=max_nucl
        #nucl_pos cannot be less than 0 
        max_nucl=np.array([max(i,0)for i in properties['nucl_pos']])
        properties['nucl_pos']=max_nucl
        """Target area function depending age and z nuclei position"""


        r=1.6
        max_value=2.5
        c=0.25
        properties['A0'] = calc_target_area(age=properties['age'],nucl_pos=max_nucl,r=r,max_value=max_value,c=c) # target area now depends on nucl_pos
        #properties['A0_initial'] = (properties['age']+1.0)*0.5*(1.0+properties['zposn']**2) 

         #############BAETTI VID
        
        #ath BAETTI D VID
        cells.mesh , number_T1, d= cells.mesh.transition(T1_eps)  #check edges verifing T1 transition
        #UNKNOWN BUG HERE NEEDS FIXING BY AARYAN
        #if len(number_T1)>0:
         #   for ii in number_T1:
          #      index = cells.mesh.face_id_by_edge[ii]
           #     properties['T1_angle_pD'] #=np.append(properties['T1_angle_pD'],cells.mesh.edge_angle[ii])
        F = force(cells)/viscosity  #force per each cell force= targetarea+Tension+perimeter+pressure_boundary 

        #EG BAETTI VID
        len_modified=np.matrix.copy(cells.mesh.length)
        #len_modified[np.where((properties['parent_group'][cells.mesh.face_id_by_edge]==1) & np.logical_not(properties['parent_group'][cells.mesh.face_id_by_edge[cells.mesh.edges.reverse]]==0))]*=0.12
        
        tens = (0.5*cells.by_edge('Lambda', 'Lambda_boundary')/len_modified)*cells.mesh.edge_vect
        tens= tens - tens.take(cells.mesh.edges.prev, 1)
        F+=tens

        dv = dt*model.sum_vertices(cells.mesh.edges,F) #movement of the vertices using eq: viscosity*dv/dt = F
        properties['force_x'] = F[0]*viscosity
        properties['force_y'] = F[1]*viscosity   
        
        if hasattr(cells.mesh.geometry,'width'):
            expansion[0] = expansion_constant*np.average(F[0]*cells.mesh.vertices[0])*dt/(cells.mesh.geometry.width**2)
        if hasattr(cells.mesh.geometry,'height'): #Cylinder mesh doesn't have 'height' argument
            expansion[1] = np.average(F[1]*cells.mesh.vertices[1])*dt/(cells.mesh.geometry.height**2)
        cells.mesh = cells.mesh.moved(dv).scaled(1.0+expansion)
        
        yield cells       



#_clones
def definecolors(cells):
    peach = '#eed5b7'
    light_blue ='#87cefa'
    pink = '#ffc0cb'
    light_green = '#98fb98'
    import matplotlib.colors as colors
    vv=sns.color_palette("hls", 10)
    v=[colors.rgb2hex(colorrgb) for colorrgb in vv]
    palette = np.array([light_green, pink,light_green,'g','r','g','m','c','',peach])
    palette = np.array([v[1],v[0],v[1], v[1],v[4],v[5],v[6],v[7],v[8],v[9],peach])
    # palette = np.array([peach, 'b', 'r','g','r','b', peach,'w',light_blue,pink,light_green,pink,light_blue,'w'])
    #colors = np.array([1 if x==71 else 2 if x==80 else 3 if x==49 else 0 for x in cells.properties['parent']])
    #colors = cells.properties['parent_group']*np.array([1 if x in sampleset else 0 for x in cells.properties['parent']])
    colors = cells.properties['parent_group']
    return palette[colors]

"""Run simulation and save data functions"""
"""
def run_simulation(x):
    K=x[0]
    G=x[1]
    L=x[2]
    rand = np.random #np.random.RandomState(123456) #I have modified the random function because RamdomState takes always the same numbers
    mesh = init.cylindrical_hex_mesh(10,10,noise=0.2,rand=rand)
    cells = model.Cells(mesh,properties={'K':K,'Gamma':G,'P':0.0,'boundary_P':P,'Lambda':L, 'Lambda_boundary':0.5})
    cells.properties['age'] = np.random.rand(len(cells))
    force = TargetArea() + Tension() + Perimeter() + Pressure()
    history = run(simulation_with_division(cells,force,rand=rand),500.0/dt,1.0/dt)
    # model.animate_video_mpg(history)
    return history
"""
def run_simulation_INM(x, timend,rand, sim_type):
    global dt
    #sim_type 0 simulation_with_division_clone (no differentiation rate)
    #sim_type 1 simulation_with_division_clone_differentiation (all differentiation rate)
    #sim_type 2 simulation_with_division_clone_differenciation_3stripes (2 population with and without diffentiation rate)
    #sim_type 3 simulation_with_division_clone_whole_tissue_differenciation (differentiation rate everywhere)
    # print(dt)

    # parameters of the vertex model
    K=x[0]
    G=x[1]
    L=x[2]

    # parameters of the nucleus A-B stochastic dynamics
    k=x[3] #NB in aaryan_simulations this is an array
    D=x[4]

    # parameters of the crowding force
    s=x[5]
    a=x[6]

    rand1 = np.random.RandomState(123456) #I have modified the random function because RamdomState takes always the same numbers
    rand.seed(1999)
    #mesh = init.cylindrical_hex_mesh(2,2,noise=0.2,rand=rand1)
    mesh = init.toroidal_hex_mesh(20,20,noise=0.2,rand=rand1)
    cells = model.Cells(mesh,properties={'K':K,'Gamma':G,'P':0.0,'boundary_P':P,'Lambda':L, 'Lambda_boundary':0.5})
  
    cells.properties['age'] = np.random.rand(len(cells))

    force = TargetArea()  + Perimeter() + Pressure()
    
    

####eg er ad baeta thessu vid
    cells.properties['parent_group'] = np.zeros(len(cells),dtype=int) #use to draw clone
    cells.properties['parent_group'] = bin_by_xpos(cells,np.cumsum([0.475,0.5,0.475]))
 
    print("Start thermalization...")
    history1 = run(simulation_with_division(cells,force,rand=rand1),200/dt,50.0/dt)
    cells = history1[-1].copy() #last timestep in the 'thermalization' phase -> use this as the initial condition for the actual run below
    # sampleset = np.random.choice(cells.mesh.face_ids,20, replace=False) #take 7 different ramdon cells to follow clones

 
    cells.properties['parent_group'] = np.zeros(len(cells),dtype=int) #use to draw clone
    cells.properties['parent'] = cells.mesh.face_ids #save the ids to control division parents-daugthers 
    cells.properties['parent_group'] = bin_by_xpos(cells,np.cumsum([0.225, 0.15, 0.075, 0.1, 0.075, 0.15, 0.225]))

    print(f"Start simulation of type {sim_type}...")
    # sampleset = np.random.choice(cells.mesh.face_ids[cells.properties['parent_group']==1],10, replace=False) #take 10 different ramdon cells to follow clones
    if sim_type == 100:
        history = run(simulation_with_division(cells,force,rand=rand),(timend)/dt,1.0/dt)
    if sim_type == 0:
        history = run(simulation_with_division_clone(cells,force,rand=rand),(timend)/dt,1.0/dt)
        history[-1].properties['parent_group'] = np.zeros(len(history[-1].properties['parent_group']),dtype=int)
    if sim_type == 1:
        history = run(simulation_with_division_clone_differentiation(cells,force,rand=rand),(timend)/dt,1.0/dt)
        history[-1].properties['parent_group'] = np.zeros(len(history[-1].properties['parent_group']),dtype=int)+1
    if sim_type == 2:
        #we take ventral and dorsal time per phase cell cycle if we are in the 2 pop part, because pNM are ventral and pD are dorsal
        history = run(simulation_with_division_clone_differenciation_3stripes(cells,force,rand=rand),(timend)/dt,1.0/dt)
        cells.properties['parent_group'] = cells.properties['parent_group'] #+np.array([3 if x in sampleset else 0 for x in cells.properties['parent']])
    if sim_type == 3:
        history = run(simulation_with_division_clone_whole_tissue_differenciation(cells,force,rand=rand),(timend)/dt,1.0/dt)
        history[-1].properties['parent_group'] = np.zeros(len(history[-1].properties['parent_group']),dtype=int)
    # added by Rebeca:
    if sim_type == 4:
        cells.properties['k']=k
        cells.properties['D']=D
        history = run(simulation_with_division_model_1(cells,force,rand=rand),(timend)/dt,1.0/dt)
    if sim_type == 5:
        cells.properties['k']=k
        cells.properties['D']=D
        history = run(simulation_with_division_model_2(cells,force,rand=rand),(timend)/dt,1.0/dt)
    if sim_type == 6:
        cells.properties['k']=k
        cells.properties['D']=D
        cells.properties['s']=s
        cells.properties['a']=a
        history = run(simulation_with_division_model_3(cells,force,rand=rand),(timend)/dt,1.0/dt)
    if sim_type == 7:
        cells.properties['k']=k
        cells.properties['D']=D
        cells.properties['s']=s
        cells.properties['a']=a
        history = run(simulation_with_division_model_4(cells,force,rand=rand),(timend)/dt,1.0/dt)
#added by aaryan
    if sim_type==8:
        cells.properties['k']=k
        cells.properties['D']=D
        cells.properties['s']=s
        cells.properties['a']=a
        zposn_list=list(cells.properties['zposn'])
        cells.properties['nucl_pos']=np.array(zposn_list)
        history=run(simulation_with_division_model_5(cells,force,rand=rand),(timend)/dt,1.0/dt)     
    if sim_type==9: #same as type 8 but with new age eqn
        cells.properties['k']=k
        cells.properties['D']=D
        cells.properties['s']=s
        cells.properties['a']=a
        zposn_list=list(cells.properties['zposn'])
        cells.properties['nucl_pos']=np.array(zposn_list)
        history=run(simulation_with_division_model_6(cells,force,rand=rand),(timend)/dt,1.0/dt)     
    if sim_type==10: #same as type 9 but with new crowding eqn
        cells.properties['k']=k
        cells.properties['D']=D
        cells.properties['s']=s
        cells.properties['a']=a
        zposn_list=list(cells.properties['zposn'])
        cells.properties['nucl_pos']=np.array(zposn_list)
        history=run(simulation_with_division_model_7(cells,force,rand=rand),(timend)/dt,1.0/dt)
    return history
   

def calc_target_area(age,nucl_pos,r=1.6,max_value=2.5,c=0.25):
    #should work with arrays?
    return((max_value/2)/(1+np.exp(-r*(age-c)))*(1+nucl_pos**2))
    
def save_data(I,history,outputdir):

    # """Save create the folder to save the data"""
    # outputdir="output_Dorsal%s"%J+"_A_c%s_viscosity003"%int(10*A_c)
    # if not os.path.exists(outputdir): # if the folder doesn't exist create it
    #     os.makedirs(outputdir)
    """Information of Area, Perimeter, Neigbours: time_mean and end. And force in time and final age distribution"""
    force = TargetArea()  + Perimeter() + Pressure()
    generation_a=[cells.mesh.area for cells in history] # collects a list of arrays each of which contains area of each cell at any given time 
    generation_p=[cells.mesh.perimeter for cells in history]
    generation_n=[np.bincount(cells.mesh.face_id_by_edge) for cells in history] # takes a vector of faces ids -> how many edges of a particular id are = count the neighbors
    generation_f=[force(cells) for cells in history] 
    # added by Rebeca:
    generation_nucl_pos = [cells.properties['nucl_pos'] for cells in history]
    generation_nucl_pos_2 = [cells.properties['nucl_pos_2'] for cells in history]
    generation_final_nucl_pos = [cells.properties['final_nucl_pos'] for cells in history]
    generation_final_nucl_pos_2 = [cells.properties['final_nucl_pos_2'] for cells in history]
    death=[cells.empty() for cells in history]
    #mean variables in time
    area_mean=[]
    perimeter_mean=[]
    neigh_mean=[]
    force_mean=[]
    force_units=[]
    area_total=[]
    number_cells = []
    ids_face_by_edge = []
    for i in range(0,len(history)-1):
        valid=np.where(~death[i] & (generation_a[i]>0))[0] # -> valid checks which cells haven't died
        number_cells = np.append(number_cells, len(valid))
        area_mean=np.append(area_mean,np.mean(generation_a[i][valid]))
        area_total=np.append(area_total,np.sum(generation_a[i][valid]))
        perimeter_mean=np.append(perimeter_mean,np.mean(generation_p[i][valid]))
        neigh_mean=np.append(neigh_mean,np.mean(generation_n[i][valid]))
        # force_mean=np.append(force_mean,np.mean(np.sqrt(generation_f[i][0][valid]**2+generation_f[i][1][valid]**2)))
        # force_units=np.append(force_units, np.sum(np.sqrt(generation_f[i][0][valid]**2+generation_f[i][1][valid]**2))/9.0)


    #end value of Area, Perimeter and Neighbour of the tissue
    """Define files per different values of area division condition"""
    outputdirname=outputdir+'/'
    outfile_a=outputdirname+"area_mean_%0.3f"%I # Is I meant to be area_mean[i] ???
    outfile_a_total=outputdirname+"area_total_%0.3f"%I
    outfile_p=outputdirname+"perimeter_mean_%0.3f"%I
    outfile_n=outputdirname+"neigh_mean_%0.3f"%I
    # outfile_f=outputdirname+"force_mean_%0.3f"%I
    # outfile_f_units=outputdirname+"force_units_%0.3f"%I

    with open(outfile_a,"w") as tfile:
        np.savetxt(tfile,area_mean)
    with open(outfile_a_total,"w") as tfile:
        np.savetxt(tfile,area_total)
    with open(outfile_p,"w") as tfile:
        np.savetxt(tfile,perimeter_mean)
    with open(outfile_n,"w") as tfile:
        np.savetxt(tfile,neigh_mean)
    # with open(outfile_f,"w") as tfile:
    #     np.savetxt(tfile,force_mean)
    # with open(outfile_f_units,"w") as tfile:
    #     np.savetxt(tfile,force_units)
        
    vert_x = cells.mesh.vertices[0]
    vert_y = cells.mesh.vertices[1]
    ids_face_by_edge = cells.mesh.face_id_by_edge  
    outfile_a_end=outputdirname+"area_end_%0.3f"%I
    outfile_p_end=outputdirname+"perimeter_end_%0.3f"%I
    outfile_n_end=outputdirname+"neigh_end_%0.3f"%I
    outfile_age_end=outputdirname+"age_end_%0.3f"%I
    outfile_nc=outputdirname+"number_cells_%0.3f"%I
    outfile_x_vert = outputdirname + "vertices_x_end_%0.3f"%I
    outfile_y_vert = outputdirname + "vertices_y_end_%0.3f"%I
    outfile_ids_face_edge = outputdirname + "ids_face_edge_%0.3f"%I
    
    valid=np.where(~death[-1] & (generation_a[-1]>0))[0]    
    with open(outfile_a_end,"w") as tfile:
        np.savetxt(tfile,np.array(generation_a[-1][valid]))
    with open(outfile_p_end,"w") as tfile:
        np.savetxt(tfile,np.array(generation_p[-1][valid]))
    with open(outfile_n_end,"w") as tfile:
        np.savetxt(tfile,np.array(generation_n[-1][valid]))
    with open(outfile_age_end,"w") as tfile:
        np.savetxt(tfile,np.array(cells.properties['age'][valid]))
    with open(outfile_nc,"w") as tfile:
        np.savetxt(tfile,number_cells)
    with open(outfile_x_vert,"w") as tfile:
        np.savetxt(tfile, vert_x)
    with open(outfile_y_vert,"w") as tfile:
        np.savetxt(tfile, vert_y)
    with open(outfile_ids_face_edge,"w") as tfile:
        np.savetxt(tfile, ids_face_by_edge)
