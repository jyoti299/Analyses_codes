EXECNAME = stats

# Which libraries to include
#USERLIBS= -lMathCore -lMathMore -lReflexDict -lCore -lCint -lReflex -lFoam -lMinuit
USERLIBS= -lMinuit -L${BATINSTALLDIR}/lib -lBATmodels -lBATmtf -lBAT

INCDIRS = -I${ROOTSYS}/include -I${BATINSTALLDIR}/include

# Include Root's architecture file that has options for all of it's
# supported compilers
include ${ROOTSYS}/etc/Makefile.arch

# Define compiler
CC = g++

# Here's where you'd want to define any #define with -DVariableName
CFLAGS += $(INCDIRS) -O2 -g

# Everything it needs to know about libraries
LIBS =  $(USERLIBDIRS) $(USERLIBS) $(shell root-config --glibs)

# Define what to compile as well as the necessary dependencies and
# object names.
SRCS =  stats.cc binneddata.cc fit.cc statistics.cc
DEPS =  $(patsubst %.cc, %.d, $(SRCS)) 
OBJS =  $(patsubst %.cc, %.o, $(SRCS))

# Make the executable
$(EXECNAME): $(OBJS)
	@echo "making stuff with $(OBJS) : $(SRCS)!"
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) $(EXPLLINKLIBS) $(OutPutOpt) $(EXECNAME) 

# Make the objects
%.o: %.cc
	$(CC) $(CXXFLAGS) $(INCDIRS) -c $(OutPutOpt) $@ $<

# Make the dependencies
%.d: %.cc
	@echo "Generating dependencies for $<"
	@set -e; $(CXX) -M $(CFLAGS) $< \
	| sed 's%\($*\)\.o[ :]*%\1.o $@ : %g' > $@; \
	[ -s $@ ] || rm -f $@

# For debugging purposes
echo:
	@echo For Debugging:
	@echo .
	@echo SRCS: $(SRCS)
	@echo .
	@echo INCLUDES: $(INCDIRS)
	@echo .
	@echo OBJECTS: $(OBJS)
	@echo .
	@echo DEPS: $(DEPS)
	@echo .
	@echo CFLAGS: $(CFLAGS)

# Clean everything
clean:
	@rm -f core* ${DEPS} ${OBJS}

clobber: clean
	@rm -f *~ *.gif $(EXECNAME)
