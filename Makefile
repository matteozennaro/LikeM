## Your favourite compiler
CC=mpicc

## Other options

CFLAGS=-lm -lgsl -lgslcblas -Wall
OFLAGS=-Ofast
ADD_HEADERS=-I./headers -I./likelihoods

## Do not touch things below this line

DEBUGFLAG=

COSMODIR=source
MODULEDIR=modules
LIKEDIR=likelihoods

COSMOFILES=$(COSMODIR)/LikeM.o
MODULEFILES=$(MODULEDIR)/runtime_messages.o $(MODULEDIR)/vectors.o
MODULEFILES+=$(MODULEDIR)/file_manipulation.o $(MODULEDIR)/linalg.o
MODULEFILES+=$(MODULEDIR)/read_ini_file.o $(MODULEDIR)/numerics.o
MODULEFILES+=$(MODULEDIR)/statistics.o $(MODULEDIR)/mcmc_functions.o

LIKEFILES=$(LIKEDIR)/TinkerHMF.o $(LIKEDIR)/Cole05bias.o $(LIKEDIR)/Cole05pk.o
LIKEFILES+=$(LIKEDIR)/NonLocalBias_RS_like.o

all: modules likes cosmo mcmc

mcmc: $(COSMOFILES) $(MODULEFILES) $(LIKEFILES)
	$(CC) -o LikeM $(COSMOFILES) $(MODULEFILES) $(LIKEFILES) $(CFLAGS) $(OFLAGS) $(ADD_HEADERS) $(DEBUGFLAG)

cosmo: $(COSMOFILES)

$(COSMOFILES): %.o : %.c
	$(CC) -c $< -o $@ $(CFLAGS) $(OFLAGS)  $(ADD_HEADERS) $(DEBUGFLAG)

modules: $(MODULEFILES)

$(MODULEFILES): %.o : %.c
	$(CC) -c $< -o $@ $(CFLAGS) $(OFLAGS) $(ADD_HEADERS) $(DEBUGFLAG)

likes: $(LIKEFILES)

$(LIKEFILES): %.o : %.c
	$(CC) -c $< -o $@ $(CFLAGS) $(OFLAGS) $(ADD_HEADERS) $(DEBUGFLAG)

debug: DEBUGFLAG=-DDEBUG
debug: all

onlyplanck: DEBUGFLAG=-DONLYPLANCK
onlyplanck: all

debugplanck: DEBUGFLAG=-DDEBUG -DONLYPLANCK
debugplanck: all

clean:
	rm -rf $(LIKEDIR)/*.o
	rm -rf $(MODULEDIR)/*.o
	rm -rf $(COSMODIR)/*.o
	rm -f LikeM
