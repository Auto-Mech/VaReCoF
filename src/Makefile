
VPATH=${HOME}/src/rotd

# ****************** C++ Compiler *********************
#
CXX = g++

# Optimization Flags:
#
CXXFLAGS += -O3  -fno-inline -std=c++11

# Debuging Flags:
#
CXXFLAGS +=

# Warning Flags:
#
CXXFLAGS += -pedantic -Wno-long-long

# Include Directories Flags:
#
CXXFLAGS += 

# Macros
#
CXXFLAGS +=

# Linker Flags:
#
#LDFLAGS = -Wl,--export-dynamic
#LDLIBS =  -ldl -L${HOME}/lib64 -lslatec -L${HOME}/lib64/pgi -lacml -lpgftnrtl -lpgc -lg2c -Wl,-rpath,${HOME}/lib64:${HOME}/lib64/pgi
# LIBDIR=${ROOT}/extra/lib
#LDLIBS =  -ldl -L${LIBDIR} -lslatec -llapack -lblas -Wl,-rpath,${LIBDIR}:${ROOT}/lib64
LDLIBS =  -L${HOME}/lib -lslatec -lmkl_intel_lp64  -lmkl_sequential -lmkl_core -ldl -Wl,-rpath,${HOME}/lib

.PHONY: clean tags all

.o:
	${CXX} -o $@ $^ ${LDLIBS}

%.d: %.cc
	set -e; $(CXX) -MM $(CPPFLAGS) -I../common -I$(HOME)/include $< \
            | sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
            [ -s $@ ] || rm -f $@

MPI_SOURCES = force.cc random.cc comm.cc flux.cc run.cc raninit.cc
	

CXX_SOURCES = rotd.cc system.cc divsur.cc surface.cc math.cc \
	integral.cc error.cc gauss.cc expression.cc multipole.cc \
	molpro.cc sjk.cc pes.cc log.cc units.cc input.cc interflux.cc

EXE_SOURCES = multi.cc rfactor.cc sampling.cc mc_flux.cc \
	tst_test.cc convert_multi.cc ej_flux.cc cut_multi.cc convert_corr.cc

CXX_OBJECTS = $(CXX_SOURCES:.cc=.o)

MPI_OBJECTS = $(MPI_SOURCES:.cc=.o)

all: $(EXE_SOURCES:.cc=) convolute

sjk_pot:
	cd sjk_pot; make

pot_corr:
	cd pot_corr; make


force.o: force.cc
	mpicxx  -c $(CXXFLAGS) -o $@ $^

random.o: random.cc
	mpicxx  -c $(CXXFLAGS) -o $@ $^

comm.o: comm.cc
	mpicxx  -c $(CXXFLAGS) -o $@ $^

flux.o: flux.cc
	mpicxx  -c $(CXXFLAGS) -o $@ $^

run.o: run.cc
	mpicxx  -c $(CXXFLAGS) -o $@ $^

multi.o: multi.cc
	mpicxx  -c $(CXXFLAGS) -o $@ $<

sampling.o: sampling.cc
	mpicxx  -c $(CXXFLAGS) -o $@ $<

rfactor.o: rfactor.cc
	mpicxx  -c $(CXXFLAGS) -o $@ $<

multi: $(MPI_OBJECTS) $(CXX_OBJECTS) multi.o
	mpicxx  -o $@ $^ ${LDLIBS}

sampling: $(MPI_OBJECTS) $(CXX_OBJECTS) sampling.o
	mpicxx  -o $@ $^ ${LDLIBS}

rfactor: $(MPI_OBJECTS) $(CXX_OBJECTS) rfactor.o
	mpicxx  -o $@ $^ ${LDLIBS}


tst_test: $(CXX_OBJECTS) tst_test.o

convert_multi: $(CXX_OBJECTS) convert_multi.o

convert_corr: $(CXX_OBJECTS) convert_corr.o

ej_flux: $(CXX_OBJECTS) ej_flux.o

mc_flux: $(CXX_OBJECTS) mc_flux.o

cut_multi: $(CXX_OBJECTS) cut_multi.o


convolute: convolute.cc
	g++ ${CXXFLAGS} -o $@ $^

include $(CXX_SOURCES:.cc=.d) $(EXE_SOURCES:.cc=.d)

clean:
	rm -f *.o
	rm -f *~
	rm -f *.d
	rm -f a.out

tags:
	etags *.cc *.hh *.f

