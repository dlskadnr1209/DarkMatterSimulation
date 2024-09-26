# Compiler
FC = ifort

# Compiler flags
FFLAGS = -O2 -fopenmp -mcmodel=large

# Source files
SRC = main.f initialize.f density.f potential.f update.f acceleration.f fft3ser.f fftsubs-fourt.f

# Object files
OBJ = $(SRC:.f=.o)

# Executable
EXEC = dark_matter_simulation

# Default target
all: $(EXEC)

# Link the executable
$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

# Compile each .f file to .o
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

# Clean up object files and executable
clean:
	rm -f $(OBJ) $(EXEC)

# Phony targets
.PHONY: all clean

