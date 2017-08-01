"""This task demonstrates how to compute slant total electron count"""
from NequickG import NEQTime, Position, GalileoBroadcast
from NequickG_global import NequickG_global, Ray

# Input objects
TX = NEQTime(10,12)
BX = GalileoBroadcast(80,0,0)

# Nequick objects
NEQ_global = NequickG_global(TX, BX)
ray = Ray(100,40,0,1000,40,0)

print NEQ_global.sTEC(ray)
print NEQ_global.sTEC2(ray)

NEQ_global.slant_table(91,ray, 'slant_ray.dat', ex_attr=['foF1', 'foF2', 'foE'])

