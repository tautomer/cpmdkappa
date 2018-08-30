#!/bin/sh

if [[ $1 -eq 1 ]]
then
    FC='ifort'
    FLAGS='-traceback'
	flag1='-qopenmp'
	flag2='-g -check all -fpe0 -warn -debug extended'
else
    FC='gfortran'
    FLAGS='-fbacktrace'
	flag1='-fopenmp'
	flag2='-g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none \
	       -ffree-line-length-0 -fcheck=all \
           -ffpe-trap=zero,overflow,underflow -finit-real=nan'
fi
if [[ $2 -eq 1 ]]
then
    FLAGS="$FLAGS $flag1"
elif [[ $2 -eq 4 ]]
then
    FLAGS="$FLAGS $flag2"
elif [[ $2 -eq 3 ]]
then
    FLAGS="$FLAGS $flag1 $flag2"
fi
make FC="$FC" FLAGS="$FLAGS"
rm -f *.o *genmod*
