#!/bin/bash

gfortran -fbounds-check -Wall -fdefault-real-8 -fdefault-double-8 -cpp -o test test.f90 -llapack -lm