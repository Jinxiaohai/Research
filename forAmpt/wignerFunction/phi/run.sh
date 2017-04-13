#!/bin/bash

make clean
make
./meson.exe ./list/data.list ./out/file
