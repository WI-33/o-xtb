# This file is part of xtb.
# SPDX-Identifier: LGPL-3.0-or-later
#
# xtb is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# xtb is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with xtb.  If not, see <https://www.gnu.org/licenses/>.

# Modified by the o-xtb landscape task suite, 2026.

# Some user-configurable features
option(WITH_OpenMP "Enable support for shared memory parallelisation with OpenMP" TRUE)
option(WITH_TBLITE "Use tblite library as backend for xTB" TRUE)
option(WITH_CPCMX "Use CPCM-X solvation library for xTB" TRUE)
option(WITH_TESTS "Enable compilation of unit tests" TRUE)


# Set build type as CMake does not provide defaults
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(
    CMAKE_BUILD_TYPE "RelWithDebInfo"
    CACHE STRING "Build type to be used."
    FORCE
  )
  message(
    STATUS
    "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified."
  )
endif()

# Add modules to the CMake build 
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/modules")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" PARENT_SCOPE)

# specify module installation directory
install(
  DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/modules/"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# Compiler-specific configurations 
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(dialects "-fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fbacktrace")
  set(bounds "-fbounds-check -ffpe-trap=invalid,zero,overflow")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(dialects "-axAVX2 -r8 -traceback")
  set(bounds "-check all -fpe0")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  set(dialects "-Mbackslash -Mallocatable=03 -r8 -traceback")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang")
  set(dialects "-fdefault-real-8 -fdefault-double-8")
endif()

# This checkout supports -O2, not -O3, for GNU Fortran.
#
# gfortran 12.2.0 compiled the Hessian path of a numerically related source
# state (edcfbbe) **wrongly** at -O3 (frequencies -931.59/3329.78/3345.94 cm^-1
# instead of 1557.41/3677.84/3684.01), while -O0/-O1/-O2 were correct and
# bitwise identical to the reference binary; -O3 -fno-tree-vectorize fails as
# well, so vectorisation is not the mechanism.  The defect does **not**
# reproduce on this commit, but the mechanism is unknown, so the -O3 risk
# cannot be bounded.  CMake's `Release` build type is -O3 for GNU: configure
# with -DCMAKE_Fortran_FLAGS_RELEASE="-O2 -DNDEBUG" instead
# (docs/_simplify-recon/03-build-and-toolchain.md section 4.4, section 7).
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # Every configuration is inspected, not only the active build type: a
  # multi-config generator (Ninja Multi-Config, Visual Studio, Xcode) leaves
  # CMAKE_BUILD_TYPE empty and keeps the optimisation in
  # CMAKE_Fortran_FLAGS_<CONFIG>.
  set(xtb_opt_flags "${CMAKE_Fortran_FLAGS}")
  set(xtb_opt_configs "${CMAKE_BUILD_TYPE}")
  if(CMAKE_CONFIGURATION_TYPES)
    list(APPEND xtb_opt_configs ${CMAKE_CONFIGURATION_TYPES})
  endif()
  foreach(xtb_opt_config IN LISTS xtb_opt_configs)
    string(TOUPPER "${xtb_opt_config}" xtb_opt_config)
    if(DEFINED CMAKE_Fortran_FLAGS_${xtb_opt_config})
      string(APPEND xtb_opt_flags " ${CMAKE_Fortran_FLAGS_${xtb_opt_config}}")
    endif()
  endforeach()
  # CMake stores flag strings as `;` separated lists, so
  # -DCMAKE_Fortran_FLAGS="-O3;-funroll-loops" arrives as `;-O3;` and the
  # blank padded match below only sees it after the separators are normalised.
  string(REPLACE ";" " " xtb_opt_flags "${xtb_opt_flags}")
  if(xtb_opt_flags MATCHES "(^| )-O3( |$)")
    message(WARNING
      "Fortran optimisation level -O3 is not supported for this checkout: "
      "gfortran 12.2.0 gave wrong Hessian results at -O3 for a numerically "
      "related source state (mechanism unknown, "
      "docs/_simplify-recon/03-build-and-toolchain.md section 4.4). "
      "Use -DCMAKE_Fortran_FLAGS_RELEASE=\"-O2 -DNDEBUG\" with "
      "-DCMAKE_BUILD_TYPE=Release.")
  endif()
endif()

# Customize compiler flags
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}" PARENT_SCOPE)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialects}" PARENT_SCOPE)


# Populate xtb_version.fh with metadata
set(version ${PROJECT_VERSION})
set(oxtb_version ${OXTB_PRODUCT_VERSION})
# `WORKING_DIRECTORY` matters: without it the command runs in the build
# directory, so a build tree outside the repository (a scratch /tmp build, for
# instance) silently falls back to `unknown-commit` and the banner loses the
# revision it is supposed to show.
execute_process(COMMAND git rev-parse HEAD
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
  RESULT_VARIABLE git_return
  OUTPUT_VARIABLE commit)
if(git_return)
  set(commit "unknown-commit")
else()
  string(REGEX REPLACE "\n$" "" commit ${commit})
endif()
# `date` and `author` are the two build inputs a caller can pin
# (`SOURCE_DATE_EPOCH`, see that file); everything above comes from the commit.
include("${CMAKE_CURRENT_LIST_DIR}/version_metadata.cmake")
configure_file(
  "${PROJECT_SOURCE_DIR}/assets/templates/version.f90"
  "${PROJECT_BINARY_DIR}/include/xtb_version.fh"
  @ONLY
)
