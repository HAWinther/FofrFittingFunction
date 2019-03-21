#!/bin/bash
gfortran-mp-6 fit.f90 main.f90 -o test
gfortran-mp-6 fit_lin.f90 main_lin.f90 -o test_lin

