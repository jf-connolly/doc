# Variables
CC     = g++
CFLAGS = -O3 -lm -fopenmp -ffast-math
SDIR   = src
ODIR   = obj

RM := rm -rf

#-- Header files
_DEPS = archive.h artmap.h configuration.h dbase.h ensemble.h \
        multiPso.h objFtn_fam.h particle.h result.h sorting.h
DEPS  = $(patsubst %,$(SDIR)/%,$(_DEPS))

#-- Objects
_OBJ = archive.o artmap.o configuration.o dbase.o ensemble.o \
       multiPso.o multisemble.o objFtn_fam.o particle.o result.o sorting.o
OBJ  = $(patsubst %,$(ODIR)/%,$(_OBJ))

#-- Compile the objects 
$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

#-- Multi-objective ensemble
multisemble: $(OBJ)
	echo "Making all"
	$(CC) -o $@ $^ $(CFLAGS)

#-- Secondary commands
clean:
	-$(RM) $(ODIR)/*.o multisemble

.PHONY: clean

