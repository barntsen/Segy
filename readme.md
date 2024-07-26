
# SEGY reading and writing  data in the SEGY format.

SEGY is a collection of classes and functions for reading
and writing seismic trace data in the Seismic Unix and SEGY formats.
SEGY gives full control of header words and data formats through
a small number of primitives for reading/writing traces.

The package is written completely in python making it very portable,
but, unfortunately, lacking in speed. For large amounts of
data it is too slow, but for smaller datasets it is running ok.

A couple of simple examples demonstrating SEGY is found in the
Exampls directory.

The source code is in (sic) Src, and the Doc directory contains
the latest SEGY standard.



