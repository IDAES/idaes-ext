include ../Makefile.in

BOOST=/opt/homebrew/include
ASL=/Users/john/src/idaes-ext-build/coinbrew/dist/include/coin-or/asl
ASL_LIB=/Users/john/src/idaes-ext-build/coinbrew/dist/lib
#IDAES_HELMHOLTZ_DATA_PATH=/Users/john/git/idaes-ext/src/general_helmholtz/param_data/
#IDAES_HELMHOLTZ_TEST_DATA_PATH=/Users/john/git/idaes-ext/src/general_helmholtz/test_data/

#BOOST=/home/johne/local/include
#ASL=/home/johne/local/include/asl


ifeq ($(BOOST),)
    CXXFLAGS = -c $(CFLAGS) -std=c++0x -fPIC
else
    CXXFLAGS = -c $(CFLAGS) -std=c++0x -fPIC -I$(BOOST)
endif

OBJECTS = read_params.o phi.o props.o sat.o delta.o state.o props_hp.o \
	props_sp.o props_up.o props_tp.o ampl_wrap.o
.PHONY: clean

ALL: general_helmholtz_external.so test_exec

test: test_exec
	./test_exec

read_params.o: read_params.cpp
	  $(CXX) $(CXXFLAGS) read_params.cpp -o read_params.o -I$(ASL)

phi.o: phi.cpp
	  $(CXX) $(CXXFLAGS) phi.cpp -o phi.o -I$(ASL)

props.o: props.cpp
	  $(CXX) $(CXXFLAGS) props.cpp -o props.o

sat.o: sat.cpp
	  $(CXX) $(CXXFLAGS) sat.cpp -o sat.o

delta.o: delta.cpp
	  $(CXX) $(CXXFLAGS) delta.cpp -o delta.o

state.o: state.cpp
	  $(CXX) $(CXXFLAGS) state.cpp -o state.o

read_data.o: read_data.cpp
	$(CXX) $(CXXFLAGS) read_data.cpp -o read_data.o

props_hp.o: props_hp.cpp
	  $(CXX) $(CXXFLAGS) props_hp.cpp -o props_hp.o

props_sp.o: props_sp.cpp
	  $(CXX) $(CXXFLAGS) props_sp.cpp -o props_sp.o

props_up.o: props_up.cpp
	  $(CXX) $(CXXFLAGS) props_up.cpp -o props_up.o

props_tp.o: props_tp.cpp
	  $(CXX) $(CXXFLAGS) props_tp.cpp -o props_tp.o

ampl_wrap.o: ampl_wrap.cpp
	  $(CXX) $(CXXFLAGS) ampl_wrap.cpp -o ampl_wrap.o -I$(ASL)

test_util.o: test_util.cpp
	  $(CXX) $(CXXFLAGS) test_util.cpp -o test_util.o

test_exec.o: test_exec.cpp
	  $(CXX) $(CXXFLAGS) test_exec.cpp -o test_exec.o 

test_exec: $(OBJECTS) read_data.o test_exec.o test_util.o
	  $(LINK) $(OBJECTS) read_data.o test_exec.o test_util.o -o test_exec -lcoinasl -L$(ASL_LIB)

general_helmholtz_external.so: $(OBJECTS)
	  $(LINK) $(LDFLAGS) $(OBJECTS) -o general_helmholtz_external.so -lcoinasl -L$(ASL_LIB)

clean:
	rm -f *.o
	rm -f *.so
	rm -f test_exec
	rm -f test_exec.exe
