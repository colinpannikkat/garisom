CXX = g++
CXXFLAGS = -Wall -pedantic -std=c++17 -g -fsanitize=address -O3

SRC = main.cpp 01Utilities.cpp 09Plant.cpp 11Parameters.cpp 08Xylem.cpp 04Component.cpp 02Soils.cpp 07Leaf.cpp 10Carbon.cpp
OBJ = $(SRC:.cpp=.o)
EXEC = run

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ)

garisom205: 01gainriskSOM.2.0.5.cpp
	g++ ${CXXFLAGS} 01gainriskSOM.2.0.5.cpp \
	02Soils.cpp \
	04Component.cpp \
	-o garisom205

clean:
	rm -rf $(OBJ) $(EXEC)
	rm -rf garisom205