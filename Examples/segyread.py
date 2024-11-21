''' Example of reading a segy data set from file '''

## Imports
import matplotlib.pyplot as pl
import numpy as np
import sys
import segy

def main() :
  ''' Read a segy file '''

  #Open input file
  fp = open("sect.sgy","rb")

  #Read the data
  texthd,bhd,trhds,data=segy.readtrs(fp,-1,"segy") 
  fp.close()

  print("Binary header no of samples: ", bhd.ns )
  print("Binary header sampling: ",      bhd.dt )
  print("Data format: ",                 bhd.fmt)
  print("No of traces: ",                len(trhds))
  
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



