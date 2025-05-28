''' Example of writing a numpy array to a segy file '''

## Imports
import numpy as np
import sys
import segy

def main() :
  ''' Example of a segy file and converting to numpy'''


  # Create a segy binary header
  #bhd = segy.segybhd()
   
  # Create the ebcdic text header
  #texthd = bytearray(3200)

  # Create a segy/su trace header
  trhd = segy.segyhd("su")
 
  #Open input file
  fp = open("data.su","rb")

  #Read the text header
  #texthd=segy.readtexthd(fp)

  #Read the binary header
  #segy.readbhd(fp,bhd)
  #Read the data
  trhd,data=segy.readtrs(fp,-1,"su") 

  fp.close()

  print(data)

#--Usual python magic
if __name__ == "__main__":
    main()



