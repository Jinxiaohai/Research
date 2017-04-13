#!/bin/bash

make clean
make

./eccentricity.exe ./list/data.list ./out/epsilon.txt
