CC = g++
OPT = -O3
CFLAGS = -Wall -Wno-deprecated -std=c++0x
CURRENT=/home/cczhong/Codes/RNAMotifScanX
BOOST_INC=/home/cczhong/Tools/boost_1_55_0
BOOST_LIB_STATIC=/home/cczhong/Tools/boost_1_55_0/stage/lib
THREADPOOL=$(CURRENT)/Libraries/threadpool
INCPATH = -I$(BOOST_INC) -I$(THREADPOOL)
LIBPATH = -L$(BOOST_LIB_STATIC)
LIB = -lboost_iostreams -lboost_filesystem -lboost_system -lboost_thread -lboost_program_options -lz -lrt -lm -lpthread
OBJS = structural_motif.o annotated_structure.o scan_structure.o motif_graph_matching.o simulation.o parameters.o map_PDB_index.o
OBJS_CUT = $(OBJS) main_cut.o
OBJS_SCAN = $(OBJS) main_scan.o
OBJS_ALIGN = $(OBJS) main_align.o
OBJS_ALL = structural_motif.o annotated_structure.o scan_structure.o motif_graph_matching.o simulation.o parameters.o map_PDB_index.o main_cut.o main_scan.o main_align.o

all: $(OBJS_ALL)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_CUT) -o cut $(LIB)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_SCAN) -o scan $(LIB)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_ALIGN) -o align $(LIB)
cut_motif: $(OBJS_CUT)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_CUT) -o cut $(LIB)
scan_structure: $(OBJS_SCAN)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_SCAN) -o scan $(LIB)
align_two_motifs: $(OBJS_ALIGN)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_ALIGN) -o align $(LIB)
structural_motif.o: structural_motif.cc structural_motif.h align_functor.h
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) -c structural_motif.cc
annotated_structure.o: annotated_structure.cc annotated_structure.h
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) -c annotated_structure.cc
scan_structure.o: scan_structure.cc scan_structure.h
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) -c scan_structure.cc
motif_graph_matching.o: motif_graph_matching.cc motif_graph_matching.h
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) -c motif_graph_matching.cc
simulation.o: simulation.cc simulation.h
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) -c simulation.cc
parameters.o: parameters.cc parameters.h
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) -c parameters.cc
map_PDB_index.o: map_PDB_index.cc map_PDB_index.h
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) -c map_PDB_index.cc
main_cut.o: main_cut.cc
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) -c main_cut.cc
main_scan.o: main_scan.cc
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) -c main_scan.cc
main_align.o: main_align.cc
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) -c main_align.cc

new:
	-touch *.cc
	-make all
delete:
	-rm cut_motif scan_structure align_two_motifs
	-rm $(OBJS_ALL)
