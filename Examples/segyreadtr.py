''' Example of writing a numpy array to a segy file '''

## Imports
import numpy as np
import sys
import segy

def main() :
  ''' Example of a segy file and converting to numpy'''


  # Create a segy binary header
  bhd = segy.segybhd("segy1")
   
  # Create the ebcdic text header
  texthd = bytearray(3200)

  # Create a segy trace header
  trhd = segy.segyhd("segy1")
 
  #Open input file
  fp = open("data.sgy","rb")

  #Read the text header
  texthd=segy.readtexthd(fp)

  #Read the binary header
  segy.readbhd(fp,bhd)

  print("Binary header job id:  ", bhd.jid )
  print("Binary header job lno: ", bhd.lno )
  print("Binary header no of samples: ", bhd.ns )
  print("Binary header sampling: ",      bhd.dt )
  print("Data format: ",                 bhd.fmt)

  #Read the data
  trhd,data=segy.readtr(fp,trhd) 
  fp.close()

  print("*** Dump of data: ")
  print(data)

#--Usual python magic
if __name__ == "__main__":
    main()



