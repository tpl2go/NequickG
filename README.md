# Nequick-G
This is a python implementation of the Nequick-G ionospheric correction model described in European Space Agency's Galileo Ionospheric Model. Ionospheric disturbance flags are not implemented

## Installation
This module depends on numpy for computation, and matplotlib and Basemap for plotting.
They can be installed using pip via:


## Main classes
**NequickG_global:**

The NequickG model is fully defined with two types of inputs: time `(month and universal time)` and broadcast coefficients `(a01, a02, a03)`. Once the `NequickG_global` object is instantiated with these two types of input parameters, one can:

- calculate slant TEC from point A to point B

```python
TX = NEQTime(10,12)
BX = GalileoBroadcast(80,0,0)
NEQ_global = NequickG_global(TX, BX)
print NEQ_global.sTEC(100,40,0, 2000, 55, 0)
print NEQ_global.sTEC2(100,40,0, 2000, 55, 0)
 ```
 
 - instantiate local NequickG objects
```python
TX = NEQTime(10,12)
BX = GalileoBroadcast(80,0,0)
pos = Position(40,0)
NEQ_global = NequickG_global(TX, BX)
neq = NEQ_global.get_Nequick_local(pos)
```

**NequickG:**

NequickG_global is a factory for NequickG objects. It requires a position input and produces an object containing the topside and bottomside models that constitutes the NequickG model at a ground point. With this object, one can:

- view the vertical electron density profile
```python
TX = NEQTime(10,12)
BX = GalileoBroadcast(80,0,0)
pos = Position(40,0)
NEQ_global = NequickG_global(TX, BX)
neq = NEQ_global.get_Nequick_local(pos)
h = np.arange(100,1000)
print neq.electrondensity(h)
```


2.- calculate the vertical total electron content

```python
TX = NEQTime(10,12)
BX = GalileoBroadcast(80,0,0)
pos = Position(40,0)
NEQ_global = NequickG_global(TX, BX)
neq = NEQ_global.get_Nequick_local(pos)
print neq.vTEC(100,1000)
```

NequickG_global does its work by spawning many NequickG objects and collating their results. There might be room for optimisations in the future. 


## Validation
ESA's Galileo Ionospheric document comes attached with a validation data table. The data has been copied into three .dat files for convenience
