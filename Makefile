CC=g++
CFLAGS=-c -Wall -DNDEBUG -O3 -ffast-math -ffinite-math-only -I ./header/ -I ./../BBFMM2D/header/ -I ./../../ 
LDFLAGS=
SOURCES=./src/read_X_R_measurements.cpp ./../BBFMM2D/src/H2_2D_node.cpp ./../BBFMM2D/src/H2_2D_tree.cpp ./../BBFMM2D/src/kernel_base.cpp ./../BBFMM2D/src/kernel_types.cpp ./../BBFMM2D/src/read_location_H.cpp


SOURCESA=examples/FLIPACK2D_input_from_file_standard_kernel.cpp
SOURCESB=examples/FLIPACK2D_get_matrix_through_routine_standard_kernel.cpp 
SOURCESC=examples/FLIPACK2D_input_from_file_mykernel.cpp
SOURCESD=examples/FLIPACK2D_get_matrix_through_routine_mykernel.cpp

OBJECTSA=$(SOURCES:.cpp=.o) $(SOURCESA:.cpp=.o) 
OBJECTSB=$(SOURCES:.cpp=.o) $(SOURCESB:.cpp=.o) 
OBJECTSC=$(SOURCES:.cpp=.o) $(SOURCESC:.cpp=.o)
OBJECTSD=$(SOURCES:.cpp=.o) $(SOURCESD:.cpp=.o)


EXECUTABLEA=./exec/FLIPACK2D_input_from_file_standard_kernel
EXECUTABLEB=./exec/FLIPACK2D_get_matrix_through_routine_standard_kernel
EXECUTABLEC=./exec/FLIPACK2D_input_from_file_mykernel
EXECUTABLED=./exec/FLIPACK2D_get_matrix_through_routine_mykernel

input_from_file_standard_kernel: $(SOURCES) $(SOURCESA) $(EXECUTABLEA)
	
$(EXECUTABLEA): $(OBJECTSA)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSA) -o $@


get_matrix_through_routine_standard_kernel: $(SOURCES) $(SOURCESB) $(EXECUTABLEB)
	
$(EXECUTABLEB): $(OBJECTSB)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSB) -o $@


input_from_file_mykernel: $(SOURCES) $(SOURCESC) $(EXECUTABLEC)
	
$(EXECUTABLEC): $(OBJECTSC)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSC) -o $@


get_matrix_through_routine_mykernel: $(SOURCES) $(SOURCESD) $(EXECUTABLED)
	
$(EXECUTABLED): $(OBJECTSD)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSD) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(KERNEL) $(INDEX) $< -o $@

clean:
	rm -rf *.out examples/*.o exec/* src/*.o ./../BBFMM2D/src/*.o

tar:
	tar -zcvf FLIPACK2D.tar.gz ./exec ./src ./header ./examples ./Makefile ./README ./LICENSE
