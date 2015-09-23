##############################################################
# Makefile                                                   #
#                                                            #
# Author      : Cristobal Navarro <crinavar@dcc.uchile.cl>   #
# Version     : 1.0                                          #
# Date        : March 5 2014                                 #
# Discription : Works for Modular Host and Device programs	 #
##############################################################

# hierarchy, names
BIN		:= ./bin
OBJ		:= ./obj
SRC		:= ./src
DEP		:= ./dep
EXEC	:= mctrueke
DEBUG_EXEC := debug-mctrueke

# check if debug option was chosen
DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CXXFLAGS = -O3 -fopenmp -g -DDEBUG
else
	CXXFLAGS = -O3 -fopenmp -fopt-info-vec-optimized -ffast-math -funsafe-loop-optimizations -ftree-loop-if-convert-stores
	#CXXFLAGS = -O2 -fopenmp
endif

# compiler, include and lib parameters
INCD = 
LIBS = -lprofiler -lm -fopenmp

# source files
C_SOURCES	:= $(wildcard $(SRC)/*.c)
CPP_SOURCES	:= $(wildcard $(SRC)/*.cpp)
CPP_HEADERS	:= $(wildcard $(SRC)/*.h)

# object files
CPP_OBJS	:= $(patsubst $(SRC)/%.cpp, $(OBJ)/%.o, $(CPP_SOURCES))
C_OBJS		:= $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(C_SOURCES))

# dependencies
CPP_DEPS	:= $(patsubst $(OBJ)/%.o, $(DEP)/%.d,$(CPP_OBJS))
C_DEPS		:= $(patsubst $(OBJ)/%.o, $(DEP)/%.d,$(C_OBJS))

all: $(EXEC)

$(OBJ)/%.o: $(SRC)/%.cpp
	$(CXX) -c $(CXXFLAGS) $(INCD) -o $@ $<
	$(CXX) -M $(SRC)/$*.cpp -MT $(OBJ)/$*.o > $(DEP)/$*.d

$(OBJ)/%.o: $(SRC)/%.c
	$(CXX) -c $(CXXFLAGS) $(INCD) -o $@ $<
	$(CXX) -M $(SRC)/$*.c -MT $(OBJ)/$*.o > $(DEP)/$*.d

# rules 
$(EXEC): $(CPP_OBJS) $(C_OBJS)
	$(CXX) -o $(BIN)/$(EXEC) $(CPP_OBJS) $(C_OBJS) $(INCD) $(LIBS)

clean:
	rm -f $(BIN)/$(EXEC) $(OBJ)/*.o $(OBJ)/*.cuo $(DEP)/*.cud $(DEP)/*.d

# include dependency files
-include $(CPP_DEPS)
-include $(C_DEPS)
