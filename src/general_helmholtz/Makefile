include ../Makefile.in

#BOOST=/opt/homebrew/include
#ASL=/Users/john/src/idaes-ext-build/coinbrew/dist/include/coin-or/asl
BOOST=/home/johne/local/include
ASL=/home/johne/local/include/asl


ifeq ($(BOOST),)
    CXXFLAGS = -c $(CFLAGS) -std=c++0x -fPIC
else
    CXXFLAGS = -c $(CFLAGS) -std=c++0x -fPIC -I$(BOOST)
endif

#ADOLC_INC = /Users/john/src/idaes-ext-build/coinbrew/dist/include
#ADOLC_LIB = /Users/john/src/idaes-ext-build/coinbrew/dist/lib64
ADOLC_INC = /home/johne/local/include
ADOLC_LIB = /home/johne/local/lib64

OBJECTS = phi.o props.o sat.o solver.o delta.o state.o ampl_wrap.o \
	h2o.o r1234ze.o co2.o
.PHONY: clean

ALL: general_helmholtz_external.so

run_test:
	./testing

test: read_data.o testing.o testing run_test

h2o.o: components/h2o.cpp
	  $(CXX) $(CXXFLAGS) components/h2o.cpp -o h2o.o -I$(ADOLC_INC)

r1234ze.o: components/r1234ze.cpp
	  $(CXX) $(CXXFLAGS) components/r1234ze.cpp -o r1234ze.o -I$(ADOLC_INC)

co2.o: components/co2.cpp
	  $(CXX) $(CXXFLAGS) components/co2.cpp -o co2.o -I$(ADOLC_INC)

phi.o: phi.cpp
	  $(CXX) $(CXXFLAGS) phi.cpp -o phi.o -I$(ADOLC_INC)

solver.o: solver.cpp
	  $(CXX) $(CXXFLAGS) solver.cpp -o solver.o

state.o: state.cpp
	  $(CXX) $(CXXFLAGS) state.cpp -o state.o

delta.o: delta.cpp
	  $(CXX) $(CXXFLAGS) delta.cpp -o delta.o

props.o: props.cpp
	  $(CXX) $(CXXFLAGS) props.cpp -o props.o

sat.o: sat.cpp
	  $(CXX) $(CXXFLAGS) sat.cpp -o sat.o

ampl_wrap.o: ampl_wrap.cpp
	  $(CXX) $(CXXFLAGS) ampl_wrap.cpp -o ampl_wrap.o -I$(ASL)

read_data.o: test_data/read_data.cpp
	$(CXX) $(CXXFLAGS) test_data/read_data.cpp -o read_data.o

testing.o: testing.cpp
	  $(CXX) $(CXXFLAGS) testing.cpp -o testing.o

testing: $(OBJECTS)
	  $(LINK) $(OBJECTS) read_data.o testing.o -o testing -ladolc -L$(ADOLC_LIB)

general_helmholtz_external.so: $(OBJECTS)
	  $(LINK) $(LDFLAGS) $(OBJECTS) -o general_helmholtz_external.so -ladolc -L$(ADOLC_LIB)

clean:
	rm -f *.o
	rm -f *.so
	rm -f testing
