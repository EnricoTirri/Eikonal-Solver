-include Makefile.inc
#CXX=clang++
VPATH=./include
CXX=g++
#OPTFLAGS=-g
OPTFLAGS=-O3 -funroll-all-loops
#OPTFLAGS=-O3 -funroll-all-loops -pg -g # Se si vuole il profiler

VERBOSITY=NOVERBOSE
EIGEN_DIR=/usr/local/include/eigen3/# Here where you have the Eigine library headers
#PYBIND11_DIR=#
#PYBIND11_DIR=../extern/pybind11/include/
STDFLAGS=-std=c++20
CXXFLAGS+=$(STDFLAGS) $(OPTFLAGS) -fPIC #-Wall -Wno-ignored-attributes -Wno-maybe-uninitialized
DIM=3

CPPFLAGS=-DDIMENSION=$(DIM) -I./include -I$(EIGEN_DIR) -D$(VERBOSITY) 
DEPEND=make.dep
LDFLAGS=$(OPTFLAGS) -Wl,-rpath=$(PWD)
CC=$(CXX)
#
# The general setting is as follows:
# mains are identified bt main_XX.cpp
# all other files are XX.cpp
#

# get all files *.cpp
SRCS=$(wildcard *.cpp)
# get the corresponding object file
OBJS= $(SRCS:.cpp=.o)
# take out python interface
# get all headers in the working directory
HEADERS=$(wildcard *.hpp)
#
exe_sources=$(filter main%.cpp,$(SRCS))
EXEC=$(exe_sources:.cpp=)
OBJS_NOEXEC=$(filter-out main%.o,$(OBJS))

#========================== ORA LA DEFINIZIONE DEGLI OBIETTIVI
.phony= 2d 3d objs all alldim clean distclean doc

.DEFAULT_GOAL = alldim

alldim:
	$(MAKE) 3d
	$(MAKE) 2d


#$(PY_OBJ) : $(PY_SRC)
#	$(CXX) $^  -c $(CPPFLAGS) $(PY_INCLUDES) $(CXXFLAGS) 

3d:
	$(MAKE) clean
	$(MAKE) all  DIM=3
	cp $(EXEC) main_eikonal3D

2d:
	$(MAKE) clean
	$(MAKE) all  DIM=2
	cp $(EXEC) main_eikonal2D

all: $(DEPEND) $(EXEC)

objs: $(DEPEND) $(OBJS_NOEXEC)

clean:
	$(RM) -f $(EXEC) $(OBJS)

distclean:
	$(MAKE) clean
	$(RM)  ./doc $(DEPEND)
	$(RM) *.out *.bak *~ *.so
doc:
	doxygen $(DOXYFILE)

$(EXEC): $(OBJS)


$(DEPEND): $(SRCS)
	$(RM) $(DEPEND)
	for f in $(SRCS); do \
	$(CXX) $(STDFLAGS) $(CPPFLAGS) -MM $$f >> $(DEPEND); \
	done

-include $(DEPEND)




