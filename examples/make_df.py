import sys
import vertex_model
from vertex_model import simulation_parser
from vertex_model.simulation_parser import *
import time
import vertex_model.timer as timer
timestart = time.time()
i=int(sys.argv[1])
make_df(i)
print("aaryan",i)
timeend = time.time()
timer.timer(timestart,timeend)
