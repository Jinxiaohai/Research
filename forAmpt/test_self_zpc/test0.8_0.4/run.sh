#!/bin/bash
make clean
make 
./zpc.exe ./list/data.list ./out/file.root ./out/result > log.dat 
