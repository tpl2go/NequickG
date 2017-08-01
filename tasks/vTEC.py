"""This task demonstrates how to obtain vertical total electron count
and the ratio of total electrons in bottomside to that in topside"""

from NequickG import NEQTime, Position, GalileoBroadcast
from NequickG_global import NequickG_global

# Input objects
TX = NEQTime(4,12)
RX = Position(40,0)
BX = GalileoBroadcast(236.831,0,0)

# Nequick objects
NEQ_global = NequickG_global(TX, BX)
NEQ, Para = NEQ_global.get_Nequick_local(RX)

# vTEC quantities of interest
print NEQ.vTEC(100, 2000)
print NEQ.vTEC_ratio()
