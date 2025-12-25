# -----------------------
# Main executable
# -----------------------

SRC_BASE = ./gauge_field.cpp ./fermions_wd.cpp ./spinor_algebra.cpp ./linalg.cpp ./utils.cpp ./su3maximization.cpp ./minmax.cpp ./fermions_real_time.cpp
SRC = $(SRC_BASE) ./momentum.cpp
HDR = $(SRC_BASE:.cpp=.hpp)
HDR += ./color_spinor.hpp ./lattice.hpp ./color_algebra.hpp
OBJ = $(SRC:.cpp=.o)

CXX = g++
CXXFLAGS = -std=c++23 -D_GNU_SOURCE -O0 -I./ -fmax-errors=1 -fopenmp 
CXXFLAGS += -I /opt/OpenBLAS/include/ -L /opt/OpenBLAS/lib/ -I./ -L./
CXXFLAGS += -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP

#This trick is needed for Barkla, because the module system on Barkla only sets LD_LIBRARY_PATH
export LIBRARY_PATH := $(LD_LIBRARY_PATH):$(LIBRARY_PATH)

LIB = -lm -lgfortran -lboost_program_options -lopenblas -larpack

#Ben's local ARPACK++
#CXXFLAGS += -I./arpackpp/external/include
#LIB += -L./arpackpp/external/lib64 -larpack

#This trick was also needed for Barkla - on Barkla, it is assumed that ARPACK shared libraries are in the same directory as the executable
LIB += -Wl,-rpath,'$$ORIGIN'

all: main_dev static

main_dev: $(OBJ) ${HDR} ./main_dev.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJ) ./main_dev.cpp $(LIB) -o ./$@

static: $(OBJ) ${HDR} ./main_static.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJ) ./main_static.cpp $(LIB) -o ./$@

rt: $(OBJ) ${HDR} ./main_real_time.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJ) ./main_real_time.cpp $(LIB) -o ./$@

minmax_test: ./minmax.o ./minmax_test.cpp ./minmax.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) ./minmax.o ./minmax_test.cpp $(LIB) -o ./$@

%.o: %.cpp $(HDR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# -----------------------
# Test executable (gtest from source)
# -----------------------
TEST_EXE = ./test
TEST_SRC = ./test.cpp

TEST_DEV_EXE = ./test_dev
TEST_DEV_SRC = ./test_dev.cpp

GTEST_SRC = ./extern/googletest/googletest/src/gtest-all.cc
GTEST_INCLUDE = -I./extern/googletest/googletest/include -I./extern/googletest/googletest -pthread

# Create a list of application object files, EXCLUDING the main entry point
APP_OBJ = $(filter-out ./main_dev.o, $(OBJ))

# The test executable needs the test code, gtest, and all application code (minus main)
TEST_OBJ = test.o gtest-all.o $(APP_OBJ)
TEST_DEV_OBJ = test_dev.o gtest-all.o $(APP_OBJ)

test: $(TEST_EXE)
test_dev: $(TEST_DEV_EXE)

$(TEST_EXE): $(TEST_OBJ)
	$(CXX) $(CXXFLAGS) $^ $(LIB) -pthread -o $@

$(TEST_DEV_EXE): $(TEST_DEV_OBJ)
	$(CXX) $(CXXFLAGS) $^ $(LIB) -pthread -o $@

# Compile test.cpp with gtest include
test.o: $(TEST_SRC)
	$(CXX) $(CXXFLAGS) $(GTEST_INCLUDE) -c $< -o $@

test_dev.o: $(TEST_DEV_SRC)
	$(CXX) $(CXXFLAGS) $(GTEST_INCLUDE) -c $< -o $@

# Compile gtest-all.cc
gtest-all.o: $(GTEST_SRC)
	$(CXX) $(CXXFLAGS) $(GTEST_INCLUDE) -c $< -o $@

# -----------------------
# Clean
# -----------------------
clean:
	rm -f -v $(OBJ) test.o test_dev.o gtest-all.o $(TEST_EXE) $(TEST_DEV_EXE) ./main_dev ./static ./minmax_test
