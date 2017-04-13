#!/bin/bash

make clean
make
./ampt.exe ./list/data.list ./out/file.root > log.dat 
