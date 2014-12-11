#---[ Paths ]-------------------------------------
oPath = obj
sPath = src
iPath = include
#=================================================


#---[ Compiler ]----------------------------------
compiler = g++
#=================================================


#---[ Files ]-------------------------------------
headers = $(wildcard $(iPath)/*.hpp) $(wildcard $(iPath)/*.tpp)
sources = $(wildcard $(sPath)/*.cpp)
objects  = $(subst $(sPath)/,$(oPath)/,$(sources:.cpp=.o))
#=================================================


#---[ Flags & Libs ]------------------------------
flags =  -I$(iPath) -m64 -fopenmp
libs   =  -L/usr/lib -llapack -lblas -lm

# Debug Option
ifeq ($(DEBUG), 1)
	flags += -g
else
	#flags += -O3
endif
#=================================================



#---[ Single | Double ]---------------------------
DOUBLE = 1
ifeq ($(DOUBLE), 1)
	flags += -Ddatafloat=double
else
	flags += -Ddatafloat=float
endif
#=================================================


#---[ COMPILATION ]-------------------------------

main: $(objects) $(headers)
	$(compiler) -o main $(flags) $(objects) $(libs)

$(oPath)/%.o:$(sPath)/%.cpp $(wildcard $(subst $(sPath)/,$(iPath)/,$(<:.cpp=.hpp))) $(wildcard $(subst $(sPath)/,$(iPath)/,$(<:.cpp=.tpp)))
	$(compiler) -o $@ $(flags) -c $(libs) $<

clean:
	rm -f $(oPath)/*;
	rm -f main;
#=================================================
