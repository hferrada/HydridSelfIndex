SDSL_DIR=/home/ferrada
RMQ_DIR=/home/ferrada/develop/RMQs/rmq

CPP=g++ -std=c++11
CPPFLAGS=-O9 -Wall -DVERBOSE
INCLUDES=-I$(SDSL_DIR)/include/ -I$(RMQ_DIR)/includes/
LIB=$(SDSL_DIR)/lib/libsdsl.a $(SDSL_DIR)/lib/libdivsufsort.a $(SDSL_DIR)/lib/libdivsufsort64.a $(RMQ_DIR)/rmqrmmBP.a
OBJECTS=HybridSelfIndex.o
BINS=build_index

all: clean stats $(OBJECTS) $(BINS)

HybridSelfIndex.o: HybridSelfIndex.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c HybridSelfIndex.cpp

stats:
	@echo
	@echo " COMPILING hsi.a"
	@echo " ###################"
	@echo "  * Compiler flags: $(CPPFLAGS)"
	@echo "  * Include dirs: $(INCLUDES)"
	@echo "  * Lib dirs: $(LIB)"
	@echo

clean:
	@echo " [CLN] Removing object files"
	@rm -f $(OBJECTS) $(BINS)

build_index: 
	@echo " [BLD] Building hsi.a"
	ar -rvcs hsi.a $(OBJECTS) $(LIB) 

build_binary:
	@echo " [BLD] Building binary build_hsi"
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o build_hsi build_hsi.cpp $(OBJECTS) $(LIB)

load_binary:
	@echo " [BLD] Building binary loadDL_hsi"
	@$(CPP) $(CPPFLAGS) $(INCLUDES) -o load_hsi load_hsi.cpp $(OBJECTS) $(LIB)


