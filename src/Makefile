EXECUTABLES := CMC_Data
LIBS:= mpfr quadmath fftw3
#PREFIXES:=$(HOME)/local
BUILD_DIR=build/
#SRC_DIR=src/
DEP_DIR=dep/
COMMON_OBJS:= allocate_and_free_derivs FieldAndBoundEqns FuncAndJacobian Get_r_B spectral_utilities
COMMON_OBJS:= $(COMMON_OBJS) Get_Coefficients GetNewResolution IndexRoutines math_functions proper_sizes
COMMON_OBJS:= $(COMMON_OBJS) newton_LU newton_Eigen ScanAndPrint utilities cast functions solution_test 
COMMON_OBJS:= $(COMMON_OBJS) fast_spectral_utilities Laplace horizonfinder BondiMass AngularMomentum

CXX:=g++
INC:=$(patsubst %,-I%/include,$(PREFIXES))
LIB:=$(patsubst %,-I%/lib,$(PREFIXES))

CXXFLAGS:= -MMD $(INC) -O3 -Wextra -Wall -fopenmp
#CXXFLAGS:= -g -MMD $(INC) -Wextra -Wall -fopenmp #Debug
CPPFLAGS:=-std=c++0x
CFLAGS:=
LDFLAGS:= -fopenmp $(LIB)

LDLIBS:=$(patsubst %,-l%,$(LIBS))

SRCS:= $(wildcard *.c)

DEPS := $(patsubst %.c,$(DEP_DIR)%.d,$(SRCS))

.SUFFIXES:

all: $(patsubst %,%,$(EXECUTABLES))

CMC_Data: $(BUILD_DIR)main.o $(patsubst %,$(BUILD_DIR)%.o,$(COMMON_OBJS))
	$(CXX) $(LDFLAGS) $< $(patsubst %,$(BUILD_DIR)%.o,$(COMMON_OBJS)) $(LDLIBS) -o $@
		notify-send -t 10000 "Information" "Compiler fertig."

$(BUILD_DIR)%.o: %.c
	mkdir -p dep build ../run/physical_data ../run/plots ../run/solution
	$(CXX) -c $(CXXFLAGS) $(CFLAGS) $< -MF $(patsubst %.c,$(DEP_DIR)%.d,$<) -o $@

-include $(DEPS)

clean:
	rm -f -r $(BUILD_DIR) $(patsubst %,%,$(EXECUTABLES))

veryclean: clean
	rm -f $(DEPS)

	