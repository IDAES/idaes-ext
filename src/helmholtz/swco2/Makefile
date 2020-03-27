SOLIB = swco2_external.so
GUESS = swco2_guess.o
include ../../Makefile.in

ifeq ($(BOOST),)
    CXXFLAGS = -c $(CFLAGS) -std=c++0x -fPIC
else
    CXXFLAGS = -c $(CFLAGS) -std=c++0x -fPIC -I$(BOOST)
endif

OBJECTS = helmholtz_external.o helmholtz_solve.o helmholtz_phi.o \
          helmholtz_asl_funcs.o helmholtz_memo.o helmholtz_param_common.o \
					$(GUESS)

ALL: get_common helmholtz_external.so

.PHONY: get_common clean

get_common:
	cp ../common/helmholtz* ./
	echo '#include"swco2_param.h"' > helmholtz_param.h

swco2_guess.o: swco2_guess.cpp
	  $(CXX) $(CXXFLAGS) swco2_guess.cpp -o swco2_guess.o

iapws95_guess.o: iapws95_guess.cpp
	  $(CXX) $(CXXFLAGS) iapws95_guess.cpp -o iapws95_guess.o

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

helmholtz_param_common.o: helmholtz_param_common.cpp
	  $(CXX) $(CXXFLAGS) helmholtz_param_common.cpp -o helmholtz_param_common.o

helmholtz_external.so: $(OBJECTS)
	  $(LINK) $(LDFLAGS) $(OBJECTS) -o $(SOLIB)

clean:
	rm -f *.o
	rm -f *.so
	rm -f helmholtz*
