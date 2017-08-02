# Validation
ESA's Galileo Ionospheric Correction document comes with three validation data tables (High, Medium, and Low solar activities)

They have been copied into three text files for convenience: High_reference.dat, Medium_reference.dat Low_reference.dat

Each table describes inputs to the slant TEC calculation and the expected output.

This validation scripts calls the NequickG routine and compares the computed values against the expected value.
