# CMakeList.txt : CMake project for Lattice_Membrane_Monte_Carlo, include source and define
# project specific logic here.
#

# Add source to this project's executable.
add_executable (Lattice_Membrane_Monte_Carlo "Lattice_Membrane_Monte_Carlo.cpp" "Lattice_Membrane_Monte_Carlo.h" "evolution.h" "evolution.cpp" "membrane.h" "membrane.cpp" "lipid.h" "lipid.cpp" "math_functions.cpp" "math_functions.h" "qol_functions.h" "qol_functions.cpp" "energy.h" "energy.cpp" "periodicboundary.h" "periodicboundary.cpp" "parameters.h" "parameters.cpp" "output.h" "output.cpp" "system.h" "system.cpp" "mc_moves.cpp" "mc_moves.h" "initialize.cpp" "initialize.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET Lattice_Membrane_Monte_Carlo PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
