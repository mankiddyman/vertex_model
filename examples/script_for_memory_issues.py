del #to del
#to free memory after del
import ctypes
libc = ctypes.CDLL("libc.so.6")
libc.malloc_trim(0)
#to see used
# Importing the library
import psutil

# Getting % usage of virtual_memory ( 3rd field)
print('RAM memory % used:', psutil.virtual_memory()[2])

#to see what is using
# %%
from __future__ import print_function  # for Python2
import sys

local_vars = list(locals().items())
record=0
for var, obj in local_vars:
    if int(sys.getsizeof(obj))>int(record):
        record=sys.getsizeof(obj)
        print(var, sys.getsizeof(obj))