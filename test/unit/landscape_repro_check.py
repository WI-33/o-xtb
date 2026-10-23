#!/usr/bin/env python3
#
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
"""The two pinned build inputs of the version header.

`cmake/version_metadata.cmake` computes `date` and `author` for
`assets/templates/version.f90`, and both are compiled into the binary, so both
decide whether two clean builds of one commit produce the same bytes (report
209).  Two of the three assertions below are the ones report 208 could only
*guess*:

  (1) a build without `SOURCE_DATE_EPOCH` keeps showing today's date, i.e. a
      developer build did not change;
  (2) with `SOURCE_DATE_EPOCH` set the date comes from the epoch and is
      formatted in **UTC**: the epoch 82800 is 1970-01-01T23:00Z, so any local
      time zone east of UTC would print 1970/01/02 and a runner and a
      workstation would date the same epoch differently;
  (3) the author is a constant of the project: `$ENV{USERNAME}` and the build
      host's system name must not reach the header (the template used to join
      them as `'user@host'`), so the rendered line stays the same for two users
      and carries no `@` at all.

The check renders the real template with the real metadata script in
`cmake -P` script mode, so it needs no configured build tree and no engine: a
wrong `date`/`author` is caught in milliseconds instead of by a double build.

Usage:

    landscape_repro_check.py <cmake> <source-dir> <rundir>

The exit status is non-zero on the first violated property.
"""
import datetime
import os
import re
import subprocess
import sys

# One line of the rendered header: `character(len=*),parameter :: name = "value"`.
HEADER_LINE = re.compile(
    r'^character\(len=\*\),parameter :: (\w+) = "(.*)"$')

# `SOURCE_DATE_EPOCH` 82800 = 1970-01-01T23:00Z, the small epoch that exposes
# local-time formatting: UTC is still 1970/01/01, UTC+8 is 1970/01/02.  The
# second value is far enough from it that a dead variable cannot pass both.
EPOCH_TZ_EDGE = "82800"
EPOCH_LATER = "1000000000"

# What `cmake/version_metadata.cmake` must produce whatever `USERNAME` says.
EXPECTED_AUTHOR = "o-xtb"


def fail(message):
    print("FAIL: {0}".format(message))
    raise SystemExit(1)


def render(cmake, source_dir, rundir, name, environment):
    """Render `xtb_version.fh` with the real metadata script and template.

    The wrapper script sets the two fields the metadata script does not own
    (the project and the commit line), so the rendered header is exactly the
    one the build produces apart from them."""
    work = os.path.join(rundir, name)
    os.makedirs(work, exist_ok=True)
    script = os.path.join(work, "render.cmake")
    header = os.path.join(work, "xtb_version.fh")
    with open(script, "w") as handle:
        handle.write(
            'set(version "6.7.1")\n'
            'set(commit "{0}")\n'
            'set(oxtb_version "0.0.0")\n'
            'include("{1}/cmake/version_metadata.cmake")\n'
            'configure_file("{1}/assets/templates/version.f90" "{2}" @ONLY)\n'
            .format("0" * 40, source_dir, header))

    env = dict(os.environ)
    env.pop("SOURCE_DATE_EPOCH", None)
    env.update(environment)
    result = subprocess.run(
        [cmake, "-P", script], env=env, capture_output=True, text=True)
    if result.returncode != 0:
        fail("`cmake -P {0}` failed with status {1}:\n{2}".format(
            script, result.returncode, result.stdout + result.stderr))

    values = {}
    with open(header) as handle:
        for line in handle:
            match = HEADER_LINE.match(line.strip())
            if match:
                values[match.group(1)] = match.group(2)
    for field in ("date", "author"):
        if field not in values:
            fail("the rendered header has no `{0}` line: {1!r}".format(
                field, values))
    return values


def main():
    if len(sys.argv) != 4:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2
    cmake, source_dir, rundir = sys.argv[1:]

    # (1) A build without `SOURCE_DATE_EPOCH` still shows today's date.
    plain = render(cmake, source_dir, rundir, "plain", {})
    today = datetime.date.today().strftime("%Y/%m/%d")
    if plain["date"] != today:
        fail("a build without SOURCE_DATE_EPOCH printed the date {0!r}, "
             "expected today ({1!r})".format(plain["date"], today))

    # (2) The epoch decides the date, and it is formatted in UTC.
    edge = render(cmake, source_dir, rundir, "epoch-edge",
                  {"SOURCE_DATE_EPOCH": EPOCH_TZ_EDGE})
    if edge["date"] != "1970/01/01":
        fail("SOURCE_DATE_EPOCH={0} formatted as {1!r}, expected 1970/01/01 "
             "(UTC; a local time zone east of UTC gives 1970/01/02)".format(
                 EPOCH_TZ_EDGE, edge["date"]))
    later = render(cmake, source_dir, rundir, "epoch-later",
                   {"SOURCE_DATE_EPOCH": EPOCH_LATER})
    if later["date"] != "2001/09/09":
        fail("SOURCE_DATE_EPOCH={0} formatted as {1!r}, expected "
             "2001/09/09".format(EPOCH_LATER, later["date"]))

    # (3) No workstation identity reaches the header.
    authors = {}
    for name, environment in (
            ("username-unset", {}),
            ("username-alice", {"USERNAME": "alice"}),
            ("username-bob", {"USERNAME": "bob"})):
        authors[name] = render(cmake, source_dir, rundir, name,
                               environment)["author"]
    if len(set(authors.values())) != 1:
        fail("the author depends on the environment: {0!r}".format(authors))
    author = authors["username-unset"]
    if author != EXPECTED_AUTHOR:
        fail("the author is {0!r}, expected {1!r}".format(
            author, EXPECTED_AUTHOR))
    if "@" in author or "@origin@" in author:
        fail("the author {0!r} still carries a user or host field".format(author))

    print("PASS: no SOURCE_DATE_EPOCH -> today ({0})".format(today))
    print("PASS: SOURCE_DATE_EPOCH={0} -> {1} (UTC), {2} -> {3}".format(
        EPOCH_TZ_EDGE, edge["date"], EPOCH_LATER, later["date"]))
    print("PASS: author is {0!r} with USERNAME unset, alice and bob".format(
        author))
    return 0


if __name__ == "__main__":
    sys.exit(main())
