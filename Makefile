CPP=g++
CPPFLAGS=-Wall -std=c++11 -O3
MKDIR=mkdir -p

SRCPATH=src
LIBPATH=lib
BINPATH=bin
INPUTPATH=input
OUTPUTPATH=output
BUILD_DIRS=$(LIBPATH) $(BINPATH) $(INPUTPATH) $(OUTPUTPATH)

EX_LIBS=fftw3f
IN_LIBS=fieldio

CPP_INC_LIBS=$(foreach lib,$(EX_LIBS),-l$(lib)) $(foreach lib,$(IN_LIBS),-l$(lib))
LIB_SO=$(foreach lib,$(IN_LIBS),lib$(lib).so)

EXE=main invert_pres find_min
EXE_OUT=$(foreach exe,$(EXE),$(exe).out)

MKFIELDS=$(shell find $(SOURCEDIR) -name 'makefield-*.cpp');

.DEFAULT_GOAL := all

# VPATH specifies the search path
VPATH=$(SRCPATH)

lib%.so: %.cpp
	$(CPP) $(CPPFLAGS) -shared -fPIC -o $(LIBPATH)/$@ $<

%.out: %.cpp
	$(CPP) $(CPPFLAGS) $< -L$(LIBPATH) $(CPP_INC_LIBS) -o $(BINPATH)/$@

main.out: main-shallow-water-pseudospectral.cpp
	$(CPP) $(CPPFLAGS) $< -L$(LIBPATH) $(CPP_INC_LIBS) -o $(BINPATH)/$@

$(BUILD_DIRS):
	$(MKDIR) $@

.PHONY: dirs
dirs: $(BUILD_DIRS)

.PHONY: clean
clean:
	rm -f log
	for dir in $(BUILD_DIRS); do \
		rm -rf $$dir; \
	done

.PHONY: exe
exe: $(EXE_OUT)

.PHONY: libs
libs: $(LIB_SO)

.PHONY: makefield
makefield: makefield-elliptic-vortex.out makefield-const-vortex.out makefield-gaussian.out makefield-Kuo2004.out

.PHONY: all
all: | dirs libs exe makefield vort_src_input.out
