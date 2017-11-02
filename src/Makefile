MAIN_DIR = ..

include $(MAIN_DIR)/make.inc

SRC = $(wildcard *.F90)
OBJ = $(patsubst %.F90,%.o,$(SRC))

lib:	$(OBJ)
	$(ARCH) $(ARFLAG) $(MAIN_DIR)/$(FSI_LIB) $(OBJ)
	$(RANLIB) $(MAIN_DIR)/$(FSI_LIB)

%.o:    %.F90
	$(FC) $(FC_FLAGS) -c $<

clean:
	rm -f $(MAIN_DIR)/libfsi.a *.o 
