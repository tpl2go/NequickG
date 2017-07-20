# Nequick-G
This is a python implementation of the Nequick-G ionospheric correction model described in European Space Agency's Galileo Ionospheric Model. Ionospheric disturbance flags are not implemented

## Main classes
**NequickG_global:**

The NequickG model is fully defined with two types of inputs: time (month and universal time) and broadcast coefficients (a01, a02, a03). Once the NequickG_global object is instantiated with these two input parameters, one can:

1. view the slant TEC view of the sky from a ground viewer 
2. view the vertical TEC map across the globe 
3. calculate slant TEC from point A to point B

**NequickG:**

NequickG_global is a factory for NequickG objects. It requires a position input and produces an object containing the topside and bottomside model that constitutes the NequickG model at a ground point. With this object, one can:

1. view the vertical electron density profile
2. calculate the vertical TEC at that point

NequickG_global does its work by spawning many NequickG objects and collating its results



