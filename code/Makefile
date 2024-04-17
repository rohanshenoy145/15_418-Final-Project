APP_NAME=boruvka_ds

OBJS=boruvka_ds.cpp

HDRS = DisjointSets.h


CXX = g++
CXXFLAGS = -Wall -O3 -std=c++17 -m64 -I. -fopenmp -Wno-unknown-pragmas

all: $(APP_NAME)

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $<

clean:
	/bin/rm -rf *~ *.o $(APP_NAME) *.class
