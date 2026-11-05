# Solver versions used for all simulations

| component | commit | repository |
|---|---|---|
| PeleLMeX | `111026b1eafbffd69d77f0e692c3ea27751e0f79` | https://github.com/AMReX-Combustion/PeleLMeX.git |
| PelePhysics | `10a59e5e18bd728f66f1f87ae1588318f8618dc7` | https://github.com/AMReX-Combustion/PelePhysics.git |
| AMReX | `5a08d94fc6ab8acc24ad62eedbe0f87b832ca9d6` | https://github.com/AMReX-Codes/amrex.git |
| SUNDIALS | `0eff39663606f2ff280c4059a947ed62ae38180a` | https://github.com/LLNL/sundials.git |
| AMReX-Hydro | `1c2b88c2648d8c5ec7ef96558cc41d05e7f9c0e9` | https://github.com/AMReX-Fluids/AMReX-Hydro.git |

PelePhysics carries the local patch `pelephysics_transport.patch` (transport parameters) and the CEPTR-generated mechanisms in `PelePhysics_Mechanisms/` (copy them to `$PELE_HOME/Submodules/PelePhysics/Mechanisms/`).
