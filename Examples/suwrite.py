''' Example of writing a numpy array to a segy file '''

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

  print("*** Dump of data: ")
  print(data)

  # Create a segy/su trace header
  trhd = segy.segyhd("su")
  trhd.ns = ns     #No of samples in trace
  trhd.dt = 2000   #Sampling interval (in micro sec)
 
  #Open output file
  fp = open("data.su","wb")


  #Write the data
  for i in range(0,nt) :  
    segy.writetr(fp,trhd,data[i,:]) 

  fp.close()

#--Usual python magic
if __name__ == "__main__":
    main()
