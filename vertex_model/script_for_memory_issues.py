
from __future__ import print_function  # for Python2
import ctypes
import sys
import psutil
#%%
def free_memory():

    #to free memory after del
    
    libc = ctypes.CDLL("libc.so.6")
    libc.malloc_trim(0)
#to see used
#%%
# Importing the library

def check_memory():

# Getting % usage of virtual_memory ( 3rd field)
    print('RAM memory % used:', psutil.virtual_memory()[2])
    local_vars = list(locals().items())
    record=0
    for var, obj in local_vars:
        if int(sys.getsizeof(obj))>int(record):
            record=sys.getsizeof(obj)
            print(var, sys.getsizeof(obj))
# %%
free_memory()
check_memory()