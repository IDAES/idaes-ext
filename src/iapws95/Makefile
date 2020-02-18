SOLIB = iapws95_external.so
include ../Makefile.in

ifeq ($(BOOST),)
    CXXFLAGS = -c $(CFLAGS) -std=c++11 -fPIC
else
    CXXFLAGS = -c $(CFLAGS) -std=c++11 -fPIC -I$(BOOST)
endif

OBJECTS = helmholtz_external.o helmholtz_solve.o helmholtz_phi.o helmholtz_asl_funcs.o helmholtz_memo.o

ALL: helmholtz_external.so

helmholtz_external.o: helmholtz_external.cpp
		$(CXX) $(CXXFLAGS) helmholtz_external.cpp -o helmholtz_external.o

helmholtz_solve.o: helmholtz_solve.cpp
		$(CXX) $(CXXFLAGS) helmholtz_solve.cpp -o helmholtz_solve.o

helmholtz_phi.o: helmholtz_phi.cpp
	  $(CXX) $(CXXFLAGS) helmholtz_phi.cpp -o helmholtz_phi.o

helmholtz_memo.o: helmholtz_memo.cpp
	  $(CXX) $(CXXFLAGS) helmholtz_memo.cpp -o helmholtz_memo.o

helmholtz_asl_funcs.o: helmholtz_asl_funcs.cpp
	  $(CXX) $(CXXFLAGS) helmholtz_asl_funcs.cpp -o helmholtz_asl_funcs.o

helmholtz_external.so: $(OBJECTS)
	  $(LINK) $(LDFLAGS) $(OBJECTS) -o $(SOLIB)

clean:
	rm -f *.o
	rm -f *.so
