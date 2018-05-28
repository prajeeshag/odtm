# Makefile for the "mppnccombine" program
#
# Written by Hans Vahlenkamp (Hans.Vahlenkamp@noaa.gov)
# Geophysical Fluid Dynamics Laboratory / NOAA
# Princeton Forrestal Campus
# Last updated: 05/15/08

NETCDF = /gpfs1/home/Libs/INTEL/NETCDF4/netcdf-4.2.1
SOURCES = mppnccombine.c
OBJECTS = mppnccombine.o
CC = icc
##CFLAGS = -O2 -I$(NETCDF)/include
CFLAGS = -O2 -I$(NETCDF)/include -Duse_netCDF4
LFLAGS = -L$(NETCDF)/lib
LIBS = -lnetcdff -lnetcdf

mppnccombine: $(OBJECTS)
	$(CC) -o mppnccombine -s $(OBJECTS) $(CFLAGS) $(LFLAGS) $(LIBS)

clean:
	rm -f *.o mppnccombine
