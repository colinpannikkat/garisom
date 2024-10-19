CXX = g++
CXXFLAGS = -Wall -pedantic -std=c++17 -g -fsanitize=address
LDFLAGS =  -fsanitize=address

SRC = main.cpp 01Utilities.cpp 09Plant.cpp 11Parameters.cpp 08Xylem.cpp 04Component.cpp 02Soils.cpp
OBJ = $(SRC:.cpp=.o)
EXEC = main

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $(OBJ) $(LBLIBS)

garisom205: 01gainriskSOM.2.0.5.cpp
	g++ ${DEBUGFLAGS} 01gainriskSOM.2.0.5.cpp \
	02Soils.cpp \
	04Component.cpp \
	-o garisom205

clean:
	rm -rf $(OBJ) $(EXEC)