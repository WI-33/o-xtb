# This file is part of xtb.
# SPDX-Identifier: LGPL-3.0-or-later
#
# xtb is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# xtb is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with xtb.  If not, see <https://www.gnu.org/licenses/>.

# Modified by the o-xtb landscape task suite, 2026.

cmake_minimum_required(VERSION 3.17)
option(WITH_OBJECT "To build using object library" TRUE)
option(INSTALL_MODULES "Install Fortran module files to include directory." FALSE)

# Buggy CMake versions
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.27.0 AND CMAKE_VERSION VERSION_LESS 3.28.0)
   set(WITH_OBJECT FALSE)
endif()

# Setup the o-xtb Project
project(
  "o-xtb"
  VERSION "6.7.1" 
  LANGUAGES "C" "Fortran"
)

# Product version of OxTB, i.e. of the command, its banner and `--version`.
# It is deliberately independent of `PROJECT_VERSION` above: the project
# version feeds the shared library `VERSION`/`SOVERSION` (`liboxtb.so.6.7.1`,
# soname `liboxtb.so.6`), the `libxtb.*` compatibility links and the CMake
# package version file, i.e. the ABI lineage inherited from upstream xTB, which
# must not move just because the command grew a product version of its own.
set(OXTB_PRODUCT_VERSION "0.3.1")

# Name of the command line program.  The project and its libraries are called
# `o-xtb`, the installed/registered command stays `oxtb` (the `xtb` module and
# symbol prefix of the sources is deliberately not renamed).
set(XTB_EXE_NAME "oxtb")

# Base name of the installed libraries (`liboxtb.a`, `liboxtb.so.6.7.1`).
# The project is still called `o-xtb` and the command still `oxtb`; only the
# *artifact* name is shortened.  It has to be set explicitly because CMake
# derives a library name from `${PROJECT_NAME}`, which would produce the
# awkward `libo-xtb.*`.  The library **file name** is a published interface
# (the `xtb-python` wheel opens `libxtb.so` by name, and distribution recipes
# reference the archive directly), so the install step additionally creates
# compatibility symlinks under the historical `libxtb.*` names - see the
# `XTB_COMPAT_LINKS` block below the installation rules.
set(XTB_LIB_NAME "oxtb")

# Follow GNU conventions for installing directories
include(GNUInstallDirs)


# Include CMake specific configurations
add_subdirectory("cmake")
include("${PROJECT_SOURCE_DIR}/cmake/modules/xtb-utils.cmake")

# Check a specific CMake targets for xtb build & execute corresponding Find scripts
if(NOT TARGET "mctc-lib::mctc-lib")
  find_package("mctc-lib" REQUIRED)
endif()

if(NOT TARGET "tblite::tblite" AND WITH_TBLITE)
   find_package("tblite" REQUIRED)
   add_compile_definitions(WITH_TBLITE)
endif()

foreach(_xtb_dep IN ITEMS toml-f test-drive mstore multicharge dftd4 s-dftd3 tblite)
   xtb_promote_target("${_xtb_dep}")
endforeach()
unset(_xtb_dep)

if(NOT TARGET "cpcmx::cpcmx" AND WITH_CPCMX)
   find_package("cpcmx" REQUIRED)
   add_compile_definitions(WITH_CPCMX)
endif()

foreach(_xtb_dep IN ITEMS cpcmx numsa toml-f test-drive mctc-lib)
   xtb_promote_target("${_xtb_dep}")
endforeach()
unset(_xtb_dep)

if(WITH_TESTS AND NOT TARGET "test-drive::test-drive")
  find_package("test-drive" REQUIRED)
endif()

# Sources: initialize program sources (prog) and library sources (srcs) empty
set(prog)
set(srcs)

# add filenames to the list variables
add_subdirectory("src")
add_subdirectory("symmetry")


# CMake modules for third-party software

# Windows (MinGW-w64) target: two things a Linux build never needs, both
# confined to this block so the Linux configuration is untouched.
if(WIN32)
  # gfortran preprocesses the `.F90` sources but defines **no** target macros
  # at all - not `_WIN32`, not `__MINGW32__`, not `_WIN64` - on either
  # mingw-w64 thread model, although gcc defines all three.  The two source
  # level `#ifdef _WIN32` branches of this tree (`src/mctc/systools.F90`: path
  # separator `\` and path delimiter `;`; `src/prog/main.F90`: `--ceasefiles`
  # writing `NUL`) are dead code without it, which costs the executable
  # relative parameter directory of `xtb_type_environment` (the `bin\oxtb.exe`
  # next to `share\o-xtb` lookup) and `--gfn 0` would need `XTBPATH` again.
  add_compile_definitions(_WIN32)

  # `-static` (see the executable rule) needs the static OpenMP runtime: the
  # imported `OpenMP::OpenMP_Fortran` hands over `libgomp.dll.a` and
  # `libpthread.dll.a`, which are import libraries, so the link then either
  # fails on `GOMP_parallel` or leaves `libwinpthread-1.dll` in the import
  # table.  Pointing find_package at the archives is the documented
  # `-DOpenMP_gomp_LIBRARY=...` override (here for all three libraries
  # `OpenMP_Fortran_LIB_NAMES` names), so the cross build needs no extra
  # configure flag.
  if(WITH_OpenMP AND CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    foreach(omp_lib IN ITEMS gomp pthread mingwthrd)
      execute_process(
        COMMAND ${CMAKE_Fortran_COMPILER} -print-file-name=lib${omp_lib}.a
        OUTPUT_VARIABLE xtb_omp_static
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
      )
      if(IS_ABSOLUTE "${xtb_omp_static}" AND EXISTS "${xtb_omp_static}")
        set(OpenMP_${omp_lib}_LIBRARY "${xtb_omp_static}" CACHE FILEPATH
          "Static ${omp_lib} of the Fortran compiler, for a -static Windows link"
          FORCE)
      else()
        message(WARNING "no static lib${omp_lib}.a found "
          "(${xtb_omp_static}): a -static Windows link will not be "
          "self-contained")
      endif()
    endforeach()
    unset(omp_lib)
    unset(xtb_omp_static)
  endif()
endif()

if(NOT TARGET "OpenMP::OpenMP_Fortran" AND WITH_OpenMP)
  find_package("OpenMP" REQUIRED)
endif()
find_package("LAPACK" REQUIRED)
find_package("BLAS" REQUIRED)

# `src/prog/main.F90` reaches the BLAS thread setter of the running process
# through `dlopen`/`dlsym` (`--blas-threads`).  glibc 2.34 folded libdl into
# libc, so the symbol is probed first and the extra library is only added where
# the dynamic loader is not part of the standard libraries.
include(CheckSymbolExists)
check_symbol_exists(dlopen "dlfcn.h" XTB_DLOPEN_IN_LIBC)
if(NOT XTB_DLOPEN_IN_LIBC)
  find_library(XTB_DL_LIBRARY NAMES dl)
  if(XTB_DL_LIBRARY)
    set(XTB_DL_LIBRARIES "${XTB_DL_LIBRARY}")
  else()
    message(WARNING "dlopen/dlsym were not found: --blas-threads will report "
      "that no OpenBLAS or MKL symbol is available")
  endif()
endif()


if(NOT EXISTS "${PROJECT_BINARY_DIR}/include")
   file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/include")
endif()

##################
# OBJECT LIBRARY #
##################
if(WITH_OBJECT)
   
   add_library(
      "${PROJECT_NAME}-object"
      OBJECT
      ${srcs}
   )


   # customize object library
   set_target_properties(
      "${PROJECT_NAME}-object"
      PROPERTIES
      Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
      POSITION_INDEPENDENT_CODE ON
   )

   # link object library against mctc-lib &
   # & conditionally against OpenMP, tblite, cpcmx
   target_link_libraries(
      "${PROJECT_NAME}-object"
      PUBLIC
      "mctc-lib::mctc-lib"
      "$<$<BOOL:${WITH_CPCMX}>:cpcmx::cpcmx>"
      "$<$<BOOL:${WITH_TBLITE}>:tblite::tblite>"
      "$<$<BOOL:${WITH_OpenMP}>:OpenMP::OpenMP_Fortran>"
   )

   # include directories
   target_include_directories(
      "${PROJECT_NAME}-object"
      PUBLIC
      ${xtb-config-dir}
      ${PROJECT_BINARY_DIR}
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
      $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
   )

endif()

##################
# Static Library #
##################
if(WITH_OBJECT)
   add_library(
      "lib-${PROJECT_NAME}-static"
      STATIC
      $<TARGET_OBJECTS:${PROJECT_NAME}-object>
   )
else()
   add_library(
      "lib-${PROJECT_NAME}-static"
      STATIC
      ${srcs}
   )
endif()


target_link_libraries(
   "lib-${PROJECT_NAME}-static"
   PUBLIC
   ${BLAS_LIBRARIES}
   ${LAPACK_LIBRARIES}
   "$<$<BOOL:${WITH_OpenMP}>:OpenMP::OpenMP_Fortran>"
   "mctc-lib::mctc-lib"
   "$<$<BOOL:${WITH_CPCMX}>:cpcmx::cpcmx>"
   $<$<BOOL:${WITH_TBLITE}>:tblite::tblite>
)

set_target_properties(
   "lib-${PROJECT_NAME}-static"
   PROPERTIES
   Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
   ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
   POSITION_INDEPENDENT_CODE ON
   OUTPUT_NAME "${XTB_LIB_NAME}"
 )


target_include_directories(
  "lib-${PROJECT_NAME}-static"
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
)


##################
# Shared Library #
##################
if (WITH_OBJECT) 
   add_library(
      "lib-${PROJECT_NAME}-shared"
      SHARED
      $<TARGET_OBJECTS:${PROJECT_NAME}-object>
   )

   target_link_libraries(
      "lib-${PROJECT_NAME}-shared"
      PUBLIC
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES}
      "$<$<BOOL:${WITH_OpenMP}>:OpenMP::OpenMP_Fortran>"
      "mctc-lib::mctc-lib"
      "$<$<BOOL:${WITH_CPCMX}>:cpcmx::cpcmx>"
      "$<$<BOOL:${WITH_TBLITE}>:tblite::tblite>"
   )

   set_target_properties(
      "lib-${PROJECT_NAME}-shared"
      PROPERTIES
      Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
      LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
      OUTPUT_NAME "${XTB_LIB_NAME}"
      VERSION "${PROJECT_VERSION}"
      SOVERSION "${PROJECT_VERSION_MAJOR}"
   )

   target_include_directories(
      "lib-${PROJECT_NAME}-shared"
      PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
      $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
   )
endif()

###############
# Executables #
###############

add_executable(
  ${XTB_EXE_NAME}
  ${prog}
)

target_link_libraries(
   ${XTB_EXE_NAME}
   PRIVATE
   "lib-${PROJECT_NAME}-static"
   ${XTB_DL_LIBRARIES}
)

set_target_properties(
  ${XTB_EXE_NAME}
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  OUTPUT_NAME "${XTB_EXE_NAME}"
)

# A Windows build ships one self-contained executable: `-static` keeps the
# libgfortran/libgomp/libquadmath/libwinpthread DLLs out of the import table
# (KERNEL32 and msvcrt remain, both part of Windows), and the 16 MiB stack is
# the conda-recipe value the upstream Windows build uses.  The OpenMP archive
# comes from the `OpenMP_gomp_LIBRARY` override above.
#
# `--no-insert-timestamp` removes the only input of the PE image that is not
# the build itself: `ld` stamps the current time into `TimeDateStamp` of the
# COFF header, and two clean builds of one commit differed in nothing but those
# bytes (measured, report 209).  The field is zero instead - no loader reads it
# and it carries no version information - so the image is a function of its
# inputs only.
if(WIN32)
  target_link_options(
    ${XTB_EXE_NAME}
    PRIVATE
    "-static"
    "-Wl,--stack=16777216"
    "-Wl,--no-insert-timestamp"
  )
endif()

target_include_directories(
   ${XTB_EXE_NAME}
   PRIVATE 
   ${PROJECT_SOURCE_DIR}/include
   $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

# The parameter files the executable needs at run time.  They are referenced
# twice below: copied into the build tree so the executable finds them from
# the build directory, and installed into `<prefix>/share/o-xtb`.
set(
  xtb_param_files
  "${PROJECT_SOURCE_DIR}/param_gfn0-xtb.txt"
  "${PROJECT_SOURCE_DIR}/param_gfn1-xtb.txt"
  "${PROJECT_SOURCE_DIR}/param_gfn1-si-xtb.txt"
  "${PROJECT_SOURCE_DIR}/param_gfn2-xtb.txt"
  "${PROJECT_SOURCE_DIR}/param_ipea-xtb.txt"
  "${PROJECT_SOURCE_DIR}/.param_gfnff.xtb"
)

# Build tree copy of the installed data directory.  The run time fallback of
# `xtb_type_environment` derives it from the executable path (`<bindir>/../
# share/o-xtb`, here `<builddir>/share/o-xtb` because the executable sits in
# the build root), so `--gfn 0` - the one method without a compiled-in
# parameter set - works from the build tree without `XTBPATH` as well.
file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/share/${PROJECT_NAME}")
file(COPY
  ${xtb_param_files}
  DESTINATION "${PROJECT_BINARY_DIR}/share/${PROJECT_NAME}"
)

######################
# Installation rules #
######################
# API for C/C++
install(
   FILES
   "${PROJECT_SOURCE_DIR}/include/xtb.h"
   DESTINATION
   "${CMAKE_INSTALL_INCLUDEDIR}"
)

if (INSTALL_MODULES)
   install(
      DIRECTORY
      "${PROJECT_BINARY_DIR}/include/"
      DESTINATION
      "${CMAKE_INSTALL_INCLUDEDIR}"
   )
endif()

# xtb-parameters 
install(
   FILES
   ${xtb_param_files}
   DESTINATION
   "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}"
)

# Manual pages.
# The manual pages are shipped as AsciiDoc sources and as the rendered roff
# pages.  `man/oxtb.1` and `man/xcontrol.7` are generated from the sources by
# `man/tools/adoc2man.py` (stdlib only) and committed, and
# `test/unit/landscape-man-render` fails when a source was edited without
# re-rendering.  The rendered page is installed under its real name (no
# `xtb.1` alias) because `man` cannot read `.adoc`, and the AsciiDoc source is
# installed under `share/doc` **not** under `share/man`: a source left in a man
# directory makes `man` print `ignoring bogus filename ...adoc` whenever that
# directory is in `MANPATH`.
install(
   FILES
   "${PROJECT_SOURCE_DIR}/man/oxtb.1"
   DESTINATION
   "${CMAKE_INSTALL_MANDIR}/man1"
)

install(
   FILES
   "${PROJECT_SOURCE_DIR}/man/xcontrol.7"
   DESTINATION
   "${CMAKE_INSTALL_MANDIR}/man7"
)

# The AsciiDoc sources, for readers and re-rendering.  Kept out of the man
# directories on purpose (see above).
install(
   FILES
   "${PROJECT_SOURCE_DIR}/man/oxtb.1.adoc"
   "${PROJECT_SOURCE_DIR}/man/xcontrol.7.adoc"
   DESTINATION
   "${CMAKE_INSTALL_DOCDIR}/man-src"
)

# Shell completions.
# `completions/oxtb.bash`, `completions/_oxtb` and `completions/oxtb.fish` are
# generated from the `--help` option table of `src/xhelp.f90` by
# `completions/tools/gen_completions.py` and committed;
# `test/unit/landscape-completions` fails when a source was edited without
# regenerating them.  They are installed into the shell's own data directory
# under the prefix - `share/bash-completion/completions`,
# `share/zsh/site-functions` and `share/fish/vendor_completions.d` - not into
# the user's home: the per-user locations are written by `install.sh`, and a
# `zsh` needs `fpath+=(.../share/zsh/site-functions)` before `compinit` to pick
# the file up.  The bash file is installed without its `.bash` suffix, because
# bash-completion looks up the file under the command name (`oxtb`), exactly
# like `install.sh` does it for the per-user location.
install(
   FILES
   "${PROJECT_SOURCE_DIR}/completions/oxtb.bash"
   DESTINATION
   "${CMAKE_INSTALL_DATADIR}/bash-completion/completions"
   RENAME
   "${XTB_EXE_NAME}"
)

install(
   FILES
   "${PROJECT_SOURCE_DIR}/completions/_oxtb"
   DESTINATION
   "${CMAKE_INSTALL_DATADIR}/zsh/site-functions"
)

install(
   FILES
   "${PROJECT_SOURCE_DIR}/completions/oxtb.fish"
   DESTINATION
   "${CMAKE_INSTALL_DATADIR}/fish/vendor_completions.d"
)

# Build output artifacts
if (WITH_OBJECT)
   install(
      TARGETS
      "lib-${PROJECT_NAME}-static"
      "lib-${PROJECT_NAME}-shared"
      "${XTB_EXE_NAME}"
      EXPORT "${PROJECT_NAME}-targets"
      LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
      INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"   
   )
else()
   install(
      TARGETS
      "lib-${PROJECT_NAME}-static"
      "${XTB_EXE_NAME}"
      EXPORT "${PROJECT_NAME}-targets"
      ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
      INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
   )
endif()

# Compatibility library names, see `XTB_LIB_NAME` above.  Consumers that
# resolve the library by **file name** expect the historical `libxtb.*`: the
# `xtb-python` wheel loads `libxtb.so`, and recipes test for or link against
# `libxtb.a`.  The C API (`include/xtb.h`, the `xtb_*` symbols) is unchanged,
# so a symlink is a drop-in and no consumer has to be rebuilt.  These links are
# a transition aid only: **remove them** once no supported consumer resolves
# the library by file name.  (`libxtb.*` was the artifact name for the whole
# upstream history; `libo-xtb.*` existed only between the `o-xtb` rename and
# this change, so that interim spelling gets no link.)
set(XTB_COMPAT_LINKS "libxtb.a=lib${XTB_LIB_NAME}.a")
if(WITH_OBJECT)
   list(APPEND XTB_COMPAT_LINKS "libxtb.so=lib${XTB_LIB_NAME}.so")
endif()
install(CODE "
   set(_xtb_libdir \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}\")
   foreach(_xtb_compat IN ITEMS ${XTB_COMPAT_LINKS})
      string(REPLACE \"=\" \";\" _xtb_pair \"\${_xtb_compat}\")
      list(GET _xtb_pair 0 _xtb_link)
      list(GET _xtb_pair 1 _xtb_target)
      message(STATUS \"Installing compatibility symlink: \${_xtb_link} -> \${_xtb_target}\")
      execute_process(
         COMMAND \"${CMAKE_COMMAND}\" -E create_symlink
                 \"\${_xtb_target}\" \"\${_xtb_libdir}/\${_xtb_link}\"
         RESULT_VARIABLE _xtb_link_result
      )
      if(NOT _xtb_link_result EQUAL 0)
         message(WARNING \"could not create the compatibility symlink \${_xtb_link}\")
      endif()
   endforeach()
")
unset(XTB_COMPAT_LINKS)

# CMake package files
include(CMakePackageConfigHelpers)

# -- Config version file
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY AnyNewerVersion
)

# -- Config file
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# -- Install config and configVersion
install(
  FILES
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# -- Targets file
# -- This makes the project importable from the build directory
export(
   EXPORT "${PROJECT_NAME}-targets"
   FILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-targets.cmake"
)

# -- This makes the project importable from the install directory
install(
   EXPORT "${PROJECT_NAME}-targets"
   FILE "${PROJECT_NAME}-targets.cmake"
   DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
   
# Tests
if(WITH_TESTS)
   enable_testing()
   add_subdirectory("test")
endif()

#################################
# Optional command registration #
#################################
# Off by default on purpose: a build must not silently write into the user's
# PATH.  With -DXTB_REGISTER_COMMAND=ON the freshly built command line program
# is linked into "~/.local/bin" at the end of the build.  The link target is
# the executable in the build tree, so a rebuild is picked up without
# reinstalling; use the installed copy (cmake --install) for a stable path.
# `install.sh` in the source root does the same registration explicitly and
# does not need this option.
option(
  XTB_REGISTER_COMMAND
  "Create or update a symlink to the built command in ~/.local/bin"
  OFF
)
if(XTB_REGISTER_COMMAND)
  # The executable target is named after the project in the pristine tree and
  # after ${XTB_EXE_NAME} in the renamed tree, so it is looked up by target
  # type instead of by name and this block stays independent of the rename.
  get_property(xtb_register_targets DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
  foreach(xtb_register_target IN LISTS xtb_register_targets)
    get_target_property(xtb_register_type "${xtb_register_target}" TYPE)
    if(xtb_register_type STREQUAL "EXECUTABLE")
      add_custom_command(
        TARGET "${xtb_register_target}"
        POST_BUILD
        COMMAND "${CMAKE_COMMAND}"
                "-DREGISTER_TARGET=$<TARGET_FILE:${xtb_register_target}>"
                -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/register_command.cmake"
        COMMENT "Registering ${xtb_register_target} in ~/.local/bin"
        VERBATIM
      )
    endif()
  endforeach()
  unset(xtb_register_targets)
endif()
