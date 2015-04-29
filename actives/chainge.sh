#!/bin/bash

for file in ./*
do
 /usr/local/gromacs/bin/editconf -f $file -label A -o $file 
done
