''' Example of writing a numpy array to a segy file 

'''

## Imports
import numpy as np
import sys
import segy

def main() :
  ''' Example of writing a numpy array to a segy file '''

  #Create a numpy array
  nt=1;
  ns=10;
  data=np.ones((nt,ns))

  print("*** dump of data: ")
  print(data)

  # Create a segy binary header
  bhd = segy.segybhd("segy1")
  bhd.jid = 1
  bhd.lno = 8
  bhd.dt = 1
  bhd.ns = ns
  bhd.fmt = 1 #IBM float data format
   
  # Create the ebcdic text header
  texthd = bytearray(3200)

  # Create a segy trace header
  trhd = segy.segyhd("segy1")
  trhd.ns = ns
  trhd.dt = 1
  trhd.float = 1
 
  #Open output file
  fp = open("data.sgy","wb")

  #Write the text header
  segy.writetexthd(fp,texthd)

  #Write the binary header
  segy.writebhd(fp,bhd)

  #Write the data
  for i in range(0,nt) :  
    segy.writetr(fp,trhd,data[i,:]) 

  fp.close()

#--Usual python magic
if __name__ == "__main__":
    main()



