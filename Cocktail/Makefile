OBJS = mesondecay.o VecMeson.o
EXE = mesondecay

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

INCFLAGS = -I$(ROOTSYS)/include
LDFLAGS = -L$(ROOTSYS)/lib

CXX = g++
FLAGS = -std=c++11 -m32 -fno-inline -Wall -g $(INCFLAGS) $(LDFLAGS)

COMPILE = $(CXX) $(FLAGS) -c 

all: $(EXE)

$(EXE): $(OBJS)
	$(CXX) -m32 -fno-inline -o $(EXE) $(OBJS) $(ROOTFLAGS) $(ROOTLIBS)

%.o: %.cxx
	$(COMPILE)  $< 

