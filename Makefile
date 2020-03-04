DEBUG ?= 0
STATIC ?= 0

# Submodules
PWD = $(shell pwd)
SDSL_ROOT ?= ${PWD}/src/sdsl-lite

# Flags
CXX=g++
CXXFLAGS += -std=c++11 -O3 -DNDEBUG -I ${SDSL_ROOT}/include -pedantic -W -Wall
LDFLAGS += -L${SDSL_ROOT}/lib -lsdsl -ldivsufsort -ldivsufsort64 -lboost_iostreams -lboost_filesystem -lboost_system -lboost_program_options -lboost_date_time

ifeq (${STATIC}, 1)
	LDFLAGS += -static -static-libgcc -pthread -lz -llzma -lbz2
else
	LDFLAGS += -lz -llzma -lbz2
endif

ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	CXXFLAGS += -O3 -fno-tree-vectorize -DNDEBUG
endif


# External sources
SDSLSOURCES = $(wildcard src/sdsl-lite/lib/*.cpp)
SOURCES = $(wildcard src/*.h) $(wildcard src/*.cpp) 
PBASE=$(shell pwd)

# Target
TARGETS =.sdsl src/paths

all:  $(TARGETS)

.sdsl: $(SDSLSOURCES)

	cd src/sdsl-lite/ && ./install.sh ${PBASE}/src/sdsl-lite && cd ../../ && touch .sdsl

src/paths: .sdsl $(SOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

clean:
	cd src/sdsl-lite/ && ./uninstall.sh && cd ../../ && rm .sdsl
