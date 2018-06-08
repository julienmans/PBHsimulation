SHELL=/bin/bash
WORKDIR=$(shell pwd)

# Daniele
#GSL_DIR = /usr/local/Cellar/gsl/1.16/lib 
#GSL_INCLUDE = /usr/local/Cellar/gsl/1.16/include

# Francesca
#GSL_DIR = /Users/fcalore/Software/gsl-1.16/ 
#GSL_INCLUDE = /Users/fcalore/Software/gsl-1.16/

# Julien
GSL_DIR = /usr/lib
GSL_INCLUDE = /usr/include/gsl

CC=gcc
CXX=g++
CFLAGS = -g -Wall -O3 -L${GSL_DIR} -lgsl -lgslcblas -std=c++11 -fopenmp #Added the last flag for random lognormal usage
INCS=-I.  -I${GSL_INCLUDE}

${WORKDIR}/%.o : ${WORKDIR}/%.cpp
	@(echo "Compiling $(@F)")
	$(CC) $(CFLAGS) -c $< -o $@

NOBJS = ${WORKDIR}/main.o ${WORKDIR}/PBHpopulation.o ${WORKDIR}/randgen.o  ${WORKDIR}/input.o ${WORKDIR}/tinystr.o ${WORKDIR}/tinyxml.o ${WORKDIR}/tinyxmlerror.o ${WORKDIR}/tinyxmlparser.o 

all:    $(NOBJS) 
	$(CXX) -o ${WORKDIR}/PBH ${WORKDIR}/main.o ${WORKDIR}/PBHpopulation.o ${WORKDIR}/randgen.o ${WORKDIR}/input.o ${WORKDIR}/tinystr.o ${WORKDIR}/tinyxml.o ${WORKDIR}/tinyxmlerror.o ${WORKDIR}/tinyxmlparser.o $(CFLAGS) $(INCS) # I moved these two to the end (Julien)

clean:
	rm -f ${WORKDIR}/*.o ${WORKDIR}/PBH
