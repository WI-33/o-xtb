# This file is part of xtb.
#
# Copyright (C) 2026 xtb contributors
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

# Add all landscape sources with one build file.  Nested CMakeLists.txt files
# in types/, engine/ and optimize/ are deliberately avoided (PLAN section 4.1).

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/util/process.c"
  "${dir}/util/process.f90"
  "${dir}/types/settings.f90"
  "${dir}/types/results.f90"
  "${dir}/types/coords.f90"
  "${dir}/types/constraints.f90"
  "${dir}/types/log.f90"
  "${dir}/types/hessian.f90"
  "${dir}/types/json.f90"
  "${dir}/engine/engine.f90"
  "${dir}/optimize/hessian_update.f90"
  "${dir}/optimize/line_search.f90"
  "${dir}/optimize/sd.f90"
  "${dir}/optimize/cg.f90"
  "${dir}/optimize/lbfgs.f90"
  "${dir}/optimize/bfgs.f90"
  "${dir}/optimize/diis.f90"
  "${dir}/optimize/rfo.f90"
  "${dir}/optimize/driver.f90"
  "${dir}/saddle/eigenfollow.f90"
  "${dir}/saddle/bprfo.f90"
  "${dir}/saddle/prfo.f90"
  "${dir}/saddle/dimer.f90"
  "${dir}/tasks/sp.f90"
  "${dir}/tasks/opt.f90"
  "${dir}/scan/relax_scan.f90"
  "${dir}/scan/scan1d.f90"
  "${dir}/scan/scan2d.f90"
  "${dir}/tasks/scan.f90"
  "${dir}/vib/numhess.f90"
  "${dir}/vib/modes.f90"
  "${dir}/vib/thermo.f90"
  "${dir}/tasks/vib.f90"
  "${dir}/tasks/ts.f90"
  "${dir}/irc/corrector.f90"
  "${dir}/irc/gs.f90"
  "${dir}/irc/lqa.f90"
  "${dir}/irc/hpc.f90"
  "${dir}/irc/eulerpc.f90"
  "${dir}/tasks/irc.f90"
  "${dir}/path/idpp.f90"
  "${dir}/path/cineb.f90"
  "${dir}/path/neb.f90"
  "${dir}/path/autoneb.f90"
  "${dir}/path/string.f90"
  "${dir}/path/gsm.f90"
  "${dir}/landscape_driver.f90"
  "${dir}/landscape_api.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
