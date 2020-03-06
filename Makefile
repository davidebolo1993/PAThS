DEBUG ?= 0
STATIC ?= 0

# Submodules
PWD = $(shell pwd)

# Flags
CXX=g++
CXXFLAGS += -std=c++11 -O3 -DNDEBUG -I ${SDSL_ROOT}/include -pedantic -W -Wall
LDFLAGS += -lboost_iostreams -lboost_filesystem -lboost_system -lboost_program_options -lboost_date_time -lboost_serialization

# Install dir
prefix = ${PWD}
exec_prefix = $(prefix)
bindir ?= $(exec_prefix)/bin


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
SOURCES = $(wildcard src/*.h) $(wildcard src/*.cpp) 
PBASE=$(shell pwd)

# Target
BUILT_PROGRAMS = src/paths
TARGETS =src/paths

all:  $(TARGETS)

src/paths: $(SOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

install: ${BUILT_PROGRAMS}
	mkdir -p ${bindir}
	install -p ${BUILT_PROGRAMS} ${bindir}
