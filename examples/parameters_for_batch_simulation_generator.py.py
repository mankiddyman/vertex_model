from cProfile import run
import numpy as np
#types of experiments
#1) microcrowding
#2) passive
#3) active
#4) mixed constant k_g1
#5) mixed constant a

#to generate

#specify k_g1,k_s,k_g2,k_m, D, a, run, duration, sim_type

#specify sim type(s) to be run

sim_types=[12,14,13]

sim_types=list(np.repeat(sim_types,3))
runs=1
microcrowding={"k_g1":[0,0,0]*int(len(sim_types)/3),"k_s":[0,0,0]*int(len(sim_types)/3),"k_g2":[29,29,29]*int(len(sim_types)/3),"k_m":[29,29,29]*int(len(sim_types)/3),"D":[1e-5,1e-5,1e-5]*int(len(sim_types)/3),"a":[0.0001,0.001,0.001]*int(len(sim_types)/3),"run":[1,1,1]*int(len(sim_types)/3),"duration":[306,306,306]*int(len(sim_types)/3),"sim_type":sim_types}


passive_model={"k_g1":[0,0,0]*int(len(sim_types)/3),"k_s":[0,0,0]*int(len(sim_types)/3),"k_g2":[29,29,29]*int(len(sim_types)/3),"k_m":[29,29,29]*int(len(sim_types)/3),"D":[1e-5,1e-5,1e-5]*int(len(sim_types)/3),"a":[0.1,1,10]*int(len(sim_types)/3),"run":[1,1,1]*int(len(sim_types)/3),"duration":[306,306,306]*int(len(sim_types)/3),"sim_type":sim_types}


active_model={"k_g1":[1,6,12]*int(len(sim_types)/3),"k_s":[0,0,0]*int(len(sim_types)/3),"k_g2":[29,29,29]*int(len(sim_types)/3),"k_m":[29,29,29]*int(len(sim_types)/3),"D":[1e-5,1e-5,1e-5]*int(len(sim_types)/3),"a":[0,0,0]*int(len(sim_types)/3),"run":[1,1,1]*int(len(sim_types)/3),"duration":[306,306,306]*int(len(sim_types)/3),"sim_type":sim_types}


mixed_constant_k_g1={"k_g1":[6,6,6]*int(len(sim_types)/3),"k_s":[0,0,0]*int(len(sim_types)/3),"k_g2":[29,29,29]*int(len(sim_types)/3),"k_m":[29,29,29]*int(len(sim_types)/3),"D":[1e-5,1e-5,1e-5]*int(len(sim_types)/3),"a":[0.1,0.1,0.1]*int(len(sim_types)/3),"run":[1,1,1]*int(len(sim_types)/3),"duration":[306,306,306]*int(len(sim_types)/3),"sim_type":sim_types}


mixed_constant_a={"k_g1":[1,6,12]*int(len(sim_types)/3),"k_s":[0,0,0]*int(len(sim_types)/3),"k_g2":[29,29,29]*int(len(sim_types)/3),"k_m":[29,29,29]*int(len(sim_types)/3),"D":[1e-5,1e-5,1e-5]*int(len(sim_types)/3),"a":[0.1,1,10]*int(len(sim_types)/3),"run":[1,1,1]*int(len(sim_types)/3),"duration":[306,306,306]*int(len(sim_types)/3),"sim_type":sim_types}



params={"microcrowding":microcrowding,"passive_model":passive_model,"active_model":active_model,"mixed_constant_k_g1":mixed_constant_k_g1,"mixed_constant_a":mixed_constant_a}


k_g1_array=[]
k_s_array=[]
k_g2_array=[]
k_m_array=[]
D_array=[]
a_array=[]
run_array=[]
duration_array=[]
sim_type_array=[]
for dict_index in range(0,len(params)):
    dict=list(params.values())[dict_index]
    k_g1_array.extend(dict['k_g1'])
    k_s_array.extend(dict['k_s'])
    k_g2_array.extend(dict['k_g2'])
    k_m_array.extend(dict['k_m'])
    D_array.extend(dict['D'])
    a_array.extend(dict['a'])
    run_array.extend(dict['run'])
    duration_array.extend(dict['duration'])
    sim_type_array.extend(dict['sim_type'])


with open('parameters.csv', 'w') as f:
    for i in range(0,len(k_g1_array)):
        f.write(str(k_g1_array[i])+" " +str(k_s_array[i])+" "+str(k_g2_array[i])+" "+str(k_m_array[i])+" "+str(D_array[i])+" "+str(a_array[i])+" "+str(run_array[i])+" "+str(duration_array[i])+" "+str(sim_type_array[i])+"\n")
        




