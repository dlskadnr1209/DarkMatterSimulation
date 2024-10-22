#Object files
OBJ = main.o fftcc.o grvit.o \
      init.o dens.o poten.o update.o out.o

# Compiler
FC     = ifort
#
# # Compiler options
OPTS   = -openmp -mcmodel=large
#
# # Library flags (if any)
LIB    =
#
# # Directory for executable
EXEDIR = ./
#
# # Pattern rule for compiling .f files to .o
%.o: %.f
	${FC} ${OPTS} -c $<
#
#         # Link object files to create the final executable
main.x: $(OBJ)
	$(FC) ${OPTS} -o $(EXEDIR)main.x $(OBJ) $(LIB)

# Ensure that com file is present for dependencies
$(OBJ): com
#
clean:
	rm -f $(OBJ) $(EXEC)
#
# Phony targets
.PHONY: all clean
#

