# PeleLMeX, 2D planar flame with Al Kassar multimode perturbation (Phase A, Lehmann 2024 method)
# Mechanism: Mei 2021 reduced to 31 species (no AR/HE/OH*/C species); PLOG -> one mechanism per pressure:
#   make PDIR=1_000atm | 5_000atm | 10_000atm | 0_987atm (1 bar) | 9_869atm (10 bar)
DIM = 2
COMP = gnu
PRECISION = DOUBLE
USE_EB = FALSE
USE_HYPRE = FALSE
TINY_PROFILE = FALSE
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_SYCL = FALSE
DEBUG = FALSE
USE_PLASMA = FALSE
USE_MANIFOLD = FALSE

PDIR ?= 1_000atm
Chemistry_Model = Mei2021NH3/$(PDIR)
Eos_Model = Fuego
Transport_Model = Simple
USERSuffix = .Mei2021NH3_$(PDIR)

PELE_HOME ?= $(HOME)/PeleLMeX
include $(PELE_HOME)/Exec/Make.PeleLMeX
