include ../Makefile.in

BOOST=/opt/homebrew/include

ifeq ($(BOOST),)
    CXXFLAGS = -c $(CFLAGS) -std=c++0x -fPIC
else
    CXXFLAGS = -c $(CFLAGS) -std=c++0x -fPIC -I$(BOOST)
endif

ADOLC_INC = /Users/john/src/idaes-ext-build/coinbrew/dist/include
ADOLC_LIB = /Users/john/src/idaes-ext-build/coinbrew/dist/lib64

OBJECTS = phi.o props.o h2o.o sat.o testing.o solver.o delta.o

.PHONY: clean

phi.o: phi.cpp
	  $(CXX) $(CXXFLAGS) phi.cpp -o phi.o -I$(ADOLC_INC)

solver.o: solver.cpp
	  $(CXX) $(CXXFLAGS) solver.cpp -o solver.o

delta.o: delta.cpp
	  $(CXX) $(CXXFLAGS) delta.cpp -o delta.o

props.o: props.cpp
	  $(CXX) $(CXXFLAGS) props.cpp -o props.o

sat.o: sat.cpp
	  $(CXX) $(CXXFLAGS) sat.cpp -o sat.o

h2o.o: h2o.cpp
	  $(CXX) $(CXXFLAGS) h2o.cpp -o h2o.o -I$(ADOLC_INC)

testing.o: testing.cpp
	  $(CXX) $(CXXFLAGS) testing.cpp -o testing.o

testing: $(OBJECTS)
	  $(LINK) $(OBJECTS) -o testing -ladolc -L$(ADOLC_LIB)

clean:
	rm -f *.o
	rm -f *.so
	rm -f testing
