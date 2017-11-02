MAIN_DIR = ..

include $(MAIN_DIR)/make.inc

LIB = $(MAIN_DIR)/$(FSI_LIB) $(LAPACK_LIB)

all: runFSI

%: %.o
	$(FC) $(FC_FLAGS) $< -o $@ $(LIB)
%.o: %.F90
	$(FC) $(FC_FLAGS) -c $<
clean:
	rm -f runFSI
