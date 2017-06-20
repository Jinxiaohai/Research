#!/bin/bash

make clean
make

./bin/analysis ./list/data.list ./out/file
