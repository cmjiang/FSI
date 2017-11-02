# ------------------------------------------
# Fortran compiler (choose one of the three)
# ------------------------------------------
# FC = ifort
# FC = gfortran
FC = ftn

# -----
# Flags
# -----
# FC_FLAGS = -O3 -qopenmp
# FC_FLAGS = -O3 -fopenmp
FC_FLAGS = -O3 -qopenmp -target-cpu=sandybridge

# -------------------------------------
# Linking to BLAS/LAPACK or MKL library
# -------------------------------------

# (1) Linking to MKL installed on the local machine
# -------------------------------------------------
# LAPACK_LIB = -L /opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 \
#              -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

# (2) Linking to MKL or LibScie installed on NERSC
# ------------------------------------------------
# LAPACK_LIB = -mkl
LAPACK_LIB =

# -----------
# FSI library
# -----------
FSI_LIB = libfsi.a

# ------------------
# Archiver and flags
# ------------------
ARCH       = ar
ARFLAG     = cr
RANLIB     = ranlib
