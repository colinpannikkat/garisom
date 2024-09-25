# Make file for GARISOM
#
# Author: German Vargas
# Date: June 2023.
#
# This files loads and set the instructions for compiling GARISOM.
# Each line is an compilation instruction or the location of a module:

all: garisom205
 
garisom205: 01gainriskSOM.2.0.5.cpp
	g++ -std=c++11 -O3 -ffast-math 01gainriskSOM.2.0.5.cpp \
	02Soils.cpp \
	-o garisom205

clean:
	rm -rf garisom205