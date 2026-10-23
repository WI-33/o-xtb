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
#
# SPDX-Identifier: LGPL-3.0-or-later
#
# The two build inputs of `assets/templates/version.f90` a caller can pin: the
# build date and the author line.  They are compiled into the binary, so a
# reproducible build has to be able to set them; the other fields of the
# header (`version`, `oxtb_version`, `commit`) come from `PROJECT_VERSION` and
# `git` in `cmake/CMakeLists.txt` and are already determined by the commit.
#
# This file is pure computation with no side effect, so `cmake -P` can run it
# on its own: that is what `test/unit/landscape_repro_check.py` does
# (`xtb/landscape-repro-metadata`), which is why the two values are checked
# without a four minute double build.

# The build date.  `SOURCE_DATE_EPOCH` (seconds since the epoch, the
# reproducible-builds convention) wins when the environment sets it, today
# otherwise, so a normal developer build keeps the date it always showed.
#
# The epoch is formatted in **UTC**, not in CMake's default local time: local
# time makes the date depend on where the build runs as well, and the epoch
# 82800 (1970-01-01T23:00Z) formats as 1970/01/02 at UTC+8 (measured, report
# 209), i.e. a runner and a workstation would disagree about the date of the
# same epoch.
if(DEFINED ENV{SOURCE_DATE_EPOCH})
   string(TIMESTAMP date "%Y/%m/%d" UTC)
else()
   string(TIMESTAMP date "%Y/%m/%d")
endif()

# The author line.  It used to be `$ENV{USERNAME}` joined with the build
# host's system name (`@author@@@origin@` in the template), which put a
# workstation identity into every published binary and made the same commit
# build differently for two users.  It is a constant of the project now, and
# the printed line reads `compiled by o-xtb on <date>` (`src/prog/info.f90`
# and the four other sub-program headers).
set(author "o-xtb")
