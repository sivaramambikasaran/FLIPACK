CC=g++
CFLAGS=-c -Wall -DNDEBUG -O3 -ffast-math -ffinite-math-only -I ./header/ -I ./../BBFMM2D/header/ -I ./../../ 
LDFLAGS=
SOURCES= ./../BBFMM2D/src/H2_2D_Node.cpp ./../BBFMM2D/src/H2_2D_Tree.cpp ./../BBFMM2D/src/kernel_Base.cpp ./../BBFMM2D/src/kernel_Types.cpp ./src/read_metadata_FLIPACK.cpp ../BBFMM2D/src/write_Into_Binary_File.cpp

SOURCESTXT=./src/read_X_R_Measurements.cpp ./../BBFMM2D/src/read_Location_Charges.cpp

SOURCESBINARY= ./src/read_X_R_Measurements_binary.cpp ../BBFMM2D/src/read_Location_Charges_binary.cpp 

SOURCESA=examples/FLIPACK_textfile_standard_kernel.cpp
SOURCESB=examples/FLIPACK_get_matrix_through_routine_standard_kernel.cpp 
SOURCESC=examples/FLIPACK_textfile_mykernel.cpp
SOURCESD=examples/FLIPACK_get_matrix_through_routine_mykernel.cpp
SOURCESE=examples/FLIPACK_binary_file_standard_kernel.cpp
SOURCESF=examples/FLIPACK_binary_file_mykernel.cpp


OBJECTSA=$(SOURCES:.cpp=.o) $(SOURCESA:.cpp=.o) $(SOURCESTXT:.cpp=.o)
OBJECTSB=$(SOURCES:.cpp=.o) $(SOURCESB:.cpp=.o) $(SOURCESTXT:.cpp=.o)
OBJECTSC=$(SOURCES:.cpp=.o) $(SOURCESC:.cpp=.o) $(SOURCESTXT:.cpp=.o)
OBJECTSD=$(SOURCES:.cpp=.o) $(SOURCESD:.cpp=.o) $(SOURCESTXT:.cpp=.o)
OBJECTSE=$(SOURCES:.cpp=.o) $(SOURCESE:.cpp=.o) $(SOURCESBINARY:.cpp=.o)
OBJECTSF=$(SOURCES:.cpp=.o) $(SOURCESF:.cpp=.o) $(SOURCESBINARY:.cpp=.o)


EXECUTABLEA=./exec/FLIPACK_textfile_standard_kernel
EXECUTABLEB=./exec/FLIPACK_get_matrix_through_routine_standard_kernel
EXECUTABLEC=./exec/FLIPACK_textfile_mykernel
EXECUTABLED=./exec/FLIPACK_get_matrix_through_routine_mykernel
EXECUTABLEE=./exec/FLIPACK_binary_file_standard_kernel
EXECUTABLEF=./exec/FLIPACK_binary_file_mykernel

textfile_standard_kernel: $(SOURCES) $(SOURCESTXT) $(SOURCESA) $(EXECUTABLEA)

$(EXECUTABLEA): $(OBJECTSA)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSA) -o $@


get_matrix_through_routine_standard_kernel: $(SOURCES) $(SOURCESTXT) $(SOURCESB) $(EXECUTABLEB)

$(EXECUTABLEB): $(OBJECTSB)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSB) -o $@


textfile_mykernel: $(SOURCES) $(SOURCESTXT) $(SOURCESC) $(EXECUTABLEC)

$(EXECUTABLEC): $(OBJECTSC)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSC) -o $@


get_matrix_through_routine_mykernel: $(SOURCES) $(SOURCESTXT) $(SOURCESD) $(EXECUTABLED)

$(EXECUTABLED): $(OBJECTSD)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSD) -o $@

binary_file_standard_kernel: $(SOURCES) $(SOURCESBINARY) $(SOURCESE) $(EXECUTABLEE)

$(EXECUTABLEE): $(OBJECTSE)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSE) -o $@

binary_file_mykernel: $(SOURCES) $(SOURCESBINARY) $(SOURCESF) $(EXECUTABLEF)

$(EXECUTABLEF): $(OBJECTSF)
	$(CC) $(LDFLAGS) $(KERNEL) $(INDEX) $(OBJECTSF) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(KERNEL) $(INDEX) $< -o $@

clean:
	rm -rf *.out examples/*.o exec/* src/*.o ./../BBFMM2D/src/*.o

tar:
	tar -zcvf FLIPACK.tar.gz ./exec ./src ./header ./examples ./Makefile ./README.md ./LICENSE.md