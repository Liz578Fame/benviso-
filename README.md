BENVISO Program
===============

# benviso-
Program for computing binding energy and isovalue properties of nuclear EoS; printing options w/ focus on modularization


# Folders

1. routs: Contains fortran routines
   * iso_routines.f : Contains the subroutines necessary to compute EoS isovalues 
   * eb_v4.f        : Contains the subroutines necessary to compute the Binding
                      Energy and radii values based on the Semi-Emperical Mass
                      Formula
                      
2. py_libs: Contains python libraries
   * __init__.py   : Python initialization file
   * PMOD_26.py    : Python Parse library 
   * ...pyc        : Python compilied files
   
3. data: Contains data ('.don') files
   
# Files

1. isobenv.py: Main script, contains parsing and Nelder-Mead optimization routines

2. main_server.f: Main program, interfaces with the isobenv pipeline through std io

3. Makefile: the Makefile for compiling the fortran code

4. Server: fortran binary
   
