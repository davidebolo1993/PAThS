DEBUG ?= 0
STATIC ?= 0

# Submodules
PWD = $(shell pwd)
SDSL_ROOT ?= ${PWD}/src/sdsl-lite

# Flags
CXX=g++
CXXFLAGS += -std=c++11 -O3 -DNDEBUG -I ${SDSL_ROOT}/include -pedantic -W -Wall
LDFLAGS += -L${SDSL_ROOT}/lib -lsdsl -ldivsufsort -ldivsufsort64

ifeq (${STATIC}, 1)
	LDFLAGS += -static -static-libgcc -pthread -lhts -lz
else
	LDFLAGS += -lhts -lz -Wl
endif

ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	CXXFLAGS += -O3 -fno-tree-vectorize -DNDEBUG
endif
ifeq (${SDSL_ROOT}, ${PWD}/src/sdsl-lite/)
	SUBMODULES += .sdsl
endif


# External sources
SDSLSOURCES = $(wildcard src/sdsl-lite/lib/*.cpp)
SOURCES = $(wildcard src/*.cpp) $(wildcard src/*.h)
PBASE=$(shell pwd)

# Targets
BUILT_PROGRAMS = src/fqreader
TARGETS = ${SUBMODULES} ${BUILT_PROGRAMS}

all:  $(TARGETS)

.sdsl: $(SDSLSOURCES)

	cd src/sdsl-lite/ && ./install.sh ${PBASE}/src/sdsl-lite && cd ../../ && touch .sdsl

src/fqreader: ${SUBMODULES} $(SOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

clean:
	cd src/sdsl-lite/ && ./uninstall.sh && cd ../../