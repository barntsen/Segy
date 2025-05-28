''' Example of reading an su data set from file '''

## Imports
import matplotlib.pyplot as pl
import numpy as np
import sys
import segy

def main() :
  ''' Read a segy file '''

  #Open input file
  #fp = open("sect.su","rb")
  fp = open("data.su","rb")

  #Read first 240 traces 
  trhds,data=segy.readtrs(fp,240,"su") 
  print("no of samples : ", trhds[0].ns)
  print("sampling interv (microSec)  : ", trhds[0].dt)
  print("shot no: ", trhds[0].fldr)
  print("offset : ", trhds[0].offset)
  fp.close()
  
  # Plot the data
  pl.imshow(data.transpose(),cmap='gray',extent=[0,252,4.5,0])
  pl.xlabel("Trace no")
  pl.ylabel("Time (sec)")
  #Set aspect ratio
  ax=pl.gca()
  ar=3.0
  asr = 1.0/(ax.get_data_ratio()*ar)
  pl.Axes.set_aspect(ax,asr)
  pl.show()

#--Usual python magic
if __name__ == "__main__":
    main()


