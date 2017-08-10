ESA's Galileo Ionospheric Correction document provides the 13 data files
required for the NequickG algorithm. However, it will be slow if python
had to perform file I/O each time it needed a value from the data set.

A script `compile.py` has been written to convert the data text file into
python files to be imported before running the Nequick-G algorithm.

3 excel spreadsheets are provided to help associate longitudes and lattiudes
to each value in the standard modip table
