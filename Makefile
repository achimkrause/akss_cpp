#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

# define the C++ compiler to use
CPP = g++

# define any compile-time flags
CPPFLAGS = -std=c++0x -Wall -Werror -g -O0

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
LIBS = -lgmp -lgmpxx

SRCDIR = src
OBJDIR = bin

# define the C++ source files
SRCS = main.cpp

GTEST_DIR = ~/Documents/cpp/googletest/googletest


# define the C++ object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(addprefix $(OBJDIR)/,$(SRCS:.cpp=.o))

# define the executable file 
MAIN = main

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

vpath %.cpp $(SRCDIR)

all:    $(MAIN)

$(OBJS): | $(OBJDIR)

$(OBJDIR):
	@mkdir -p $@

$(MAIN): $(OBJS) 
	$(CPP) $(CPPFLAGS) $(INCLUDES) -o $(OBJDIR)/$(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJDIR)/%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c -o $@ $<

clean:
	$(RM) -r $(OBJDIR)

