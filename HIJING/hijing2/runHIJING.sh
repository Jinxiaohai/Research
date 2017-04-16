#!/bin/sh

for prod in `seq 9001 9999`
  do 
      ./test0.exe < InputFiles/input_p${prod}.txt > OutputFiles/output_${prod}.txt
      echo ${prod}
  done
   