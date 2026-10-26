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
"""Generate the command-line reference chapter of the o-xtb documentation site.

The single source of truth is the option table of ``src/xhelp.f90``: the list
of ``call add(spellings, category, summary, detail)`` records plus the two
17-element arrays ``cat_key`` / ``cat_title``.  Nothing here is a second option
table - the **extraction** of the spellings imports the ``--help`` coverage
gate, ``test/unit/landscape_help_check.py``:

  * ``HELP_LINE``           the gate's quoting regex, used to pull the quoted
                            strings of every record and of the category arrays;
  * ``SPLIT`` and ``TOKEN`` the gate's spelling tokeniser;
  * ``documented_flags``    the gate's documented spelling set (dialect A),
                            used as the zero-drift oracle.

The run writes one Markdown page per category into the output directory
(default ``docs/docs-site/docs/cli``), an ``index.md`` for the chapter, a
``nav-cli.yml`` fragment for the site ``nav`` and a ``filelist.txt``.  The
output is deterministic: no timestamp, no path and no set iteration order
reaches a file, so two runs are byte identical and ``--check`` is meaningful.

Usage::

    gen_cli_reference.py [--check] [--source PATH] [--outdir DIR]

``--check`` writes nothing and exits non-zero when the files on disk differ
from what this run would produce (a missing page, an extra page or a changed
page).  The local site gate ``tools/check_site.py`` and the CI job both use it.

The expected counts are asserted, so a change of the option table that is not
reviewed here fails loudly instead of silently shipping a stale reference.
"""

import argparse
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
REPO_ROOT = HERE.parents[3]
DEFAULT_SOURCE = REPO_ROOT / "src" / "xhelp.f90"
DEFAULT_OUTDIR = REPO_ROOT / "docs" / "docs-site" / "docs" / "cli"

# Import the gate unchanged rather than parsing the option table a second time.
sys.path.insert(0, str(REPO_ROOT / "test" / "unit"))
from landscape_help_check import (  # noqa: E402
    HELP_LINE,
    SPLIT,
    TOKEN,
    documented_flags,
)

#: The counts the gate is built on (`landscape_help_check.py` prints the same
#: numbers).  A record inserted or deleted in `src/xhelp.f90` must be
#: reflected here, so an unreviewed change of the table fails the generator.
EXPECTED_ENTRIES = 228
EXPECTED_SPELLINGS = 252
EXPECTED_CATEGORIES = 17

CAT_KEY_DECL = "cat_key(n_cat) = [character(len=16) ::"
CAT_TITLE_DECL = "cat_title(n_cat) = [character(len=64) ::"
#: One named legal-value array of the help table: `character(len=16),
#: parameter :: v_name(n) = [character(len=16) :: "a", "b"]`.  The `values=`
#: argument of a `call add` names one of these arrays.
VALUE_ARRAY_DECL = (
    r"character\(len=16\), parameter :: (v_[a-z]+)\(\d+\) = "
    r"\[character\(len=16\) ::(.*?)\]"
)
VALUE_REF = re.compile(r"values=(v_[a-z]+)")


def value_sets(source):
    """The named legal-value arrays of the help table, name -> list of values."""
    return {
        name: re.findall(r'"([^"]*)"', body)
        for name, body in re.findall(VALUE_ARRAY_DECL, source, re.S)
    }


def decode_string_array(source, decl):
    """The quoted elements of a `character(len=n), parameter :: decl = [...]`
    Fortran array, in source order."""
    start = source.index(decl)
    end = source.index("]", start)
    return [m.group(1) for m in HELP_LINE.finditer(source[start:end])]


def parse_records(source):
    """Split `src/xhelp.f90` at `call add(` and read the five fields.

    Uses the gate's `HELP_LINE` to pull the quoted strings, so the quoting rule
    is exactly the one `documented_flags` applies.  Source order is preserved.
    The legal value set is named by the optional `values=v_...` argument and is
    looked up in the named arrays by `value_sets`; a fix or a free value makes
    the list empty.
    """
    sets = value_sets(source)
    records = []
    for chunk in source.split("call add(")[1:]:
        strings = [m.group(1) for m in HELP_LINE.finditer(chunk)]
        spellings, category, summary = strings[0], strings[1], strings[2]
        detail = strings[3:]
        assert spellings.startswith("-"), spellings
        match = VALUE_REF.search(chunk)
        values = sets[match.group(1)] if match else []
        records.append((spellings, category, summary, detail, values))
    return records


def spelling_tokens(spellings):
    """The spellings of one record, split and `--[no]x`-expanded like the gate."""
    out = []
    for token in SPLIT.split(spellings.strip()):
        if not TOKEN.match(token):
            continue
        if "[no]" in token:
            out.append(token.replace("[no]", ""))
            out.append(token.replace("[no]", "no"))
        else:
            out.append(token)
    return out


def render_entry(spellings, summary, detail, values, level=3):
    """One option as Markdown: a heading of the spellings, the gloss, the
    full text and the uniform `Values:` line.  Only tokens that start with `-`
    are spellings; the argument placeholders (`[LEVEL]`, `FILE`, `REAL`, ...)
    go on a separate Arguments line so they never masquerade as options."""
    tokens = SPLIT.split(spellings.strip())
    spellings_only = [t for t in tokens if t.startswith("-")]
    placeholders = [t for t in tokens if not t.startswith("-")]
    heading = [
        "%s %s" % ("#" * level, " ".join("`%s`" % t for t in spellings_only))
    ]
    if placeholders:
        heading.append("")
        heading.append("*Arguments:* " + " ".join("`%s`" % p for p in placeholders))
    out = heading + ["", summary, ""]
    for line in detail:
        out.append(line[4:] if line.startswith("    ") else line)
    if values:
        out.append("Values: " + ", ".join(values))
    out.append("")
    return "\n".join(out).rstrip() + "\n"


def build_outputs(source):
    """Return `{relative path: text}` for the whole generated chapter.

    Pure function of `source`, so the same input always gives the same mapping
    and `--check` can compare it against the files on disk.
    """
    cat_keys = decode_string_array(source, CAT_KEY_DECL)
    cat_titles = decode_string_array(source, CAT_TITLE_DECL)
    assert len(cat_keys) == len(cat_titles) == EXPECTED_CATEGORIES, (
        cat_keys,
        cat_titles,
    )

    records = parse_records(source)
    outputs = {}

    # one page per category, in cat_key order
    for key, title in zip(cat_keys, cat_titles):
        entries = [r for r in records if r[1] == key]
        body = [
            "# %s" % title,
            "",
            "`oxtb --help %s` prints this group on the terminal." % key,
            "",
        ]
        for spellings, _, summary, detail, values in entries:
            body.append(render_entry(spellings, summary, detail, values))
        outputs[key + ".md"] = "\n".join(body)

    # the chapter index
    index = [
        "# Command line reference",
        "",
        "Generated from `src/xhelp.f90`; do not edit.  %d options across %d "
        "categories." % (len(records), len(cat_keys)),
        "",
        "## Categories",
        "",
    ]
    for key, title in zip(cat_keys, cat_titles):
        n = sum(1 for r in records if r[1] == key)
        index.append("- [%s](%s.md)  (%d options)" % (title, key, n))
    for key, title in zip(cat_keys, cat_titles):
        index += ["", "## " + title, ""]
        for spellings, _, summary, _, _ in (r for r in records if r[1] == key):
            index.append("- [`%s`](%s.md) -- %s" % (spellings.strip(), key, summary))
    outputs["index.md"] = "\n".join(index) + "\n"

    # the nav fragment and the deterministic file list
    nav = ["- Command line reference:", "    - Overview: index.md"]
    for key, title in zip(cat_keys, cat_titles):
        nav.append("    - %s: %s.md" % (title, key))
    outputs["nav-cli.yml"] = "\n".join(nav) + "\n"
    outputs["filelist.txt"] = (
        "\n".join(["index.md"] + [k + ".md" for k in cat_keys]) + "\n"
    )
    return outputs, records, cat_keys


def self_check(source, records, cat_keys, outputs):
    """Assert the gate's own invariants and print the counts.

    The spelling set built from the records must be exactly the gate's
    documented spelling set - the generator and the gate share one rule, so a
    divergence is a defect here, not a difference of opinion.
    """
    by_record = {
        t for spellings, _, _, _, _ in records for t in spelling_tokens(spellings)
    }
    by_gate = documented_flags(source, expand_no=True)
    drift = by_record ^ by_gate
    written = sum(
        outputs[k + ".md"].count("\n### ") for k in cat_keys
    )
    print("records (call add):        %d  (expected %d)" % (len(records), EXPECTED_ENTRIES))
    print("category pages:            %d  (expected %d)" % (len(cat_keys), EXPECTED_CATEGORIES))
    print(
        "generated entry headings:  %d  (== records: %s)"
        % (written, written == len(records))
    )
    print(
        "distinct spellings:        %d  (expected %d)"
        % (len(by_record), EXPECTED_SPELLINGS)
    )
    print(
        "documented_flags() == set: %s  (symmetric difference %d)"
        % (not drift, len(drift))
    )
    ok = (
        len(records) == EXPECTED_ENTRIES
        and len(cat_keys) == EXPECTED_CATEGORIES
        and written == len(records)
        and len(by_record) == EXPECTED_SPELLINGS
        and not drift
    )
    return ok, drift


def check(outdir, outputs):
    """Compare `outputs` with `outdir` on disk; print and return the drift."""
    problems = []
    on_disk = set()
    if outdir.is_dir():
        on_disk = {p.name for p in outdir.iterdir() if p.is_file()}
    for name, text in sorted(outputs.items()):
        path = outdir / name
        if not path.is_file():
            problems.append("%s is missing" % name)
        elif path.read_text() != text:
            problems.append("%s differs from the generated output" % name)
    for name in sorted(on_disk - set(outputs)):
        problems.append("%s is not a generated file (stale)" % name)
    return problems


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="write nothing; exit non-zero when the files on disk differ",
    )
    parser.add_argument("--source", default=str(DEFAULT_SOURCE), help="path to src/xhelp.f90")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="output directory")
    args = parser.parse_args(argv)

    source = Path(args.source).read_text()
    outdir = Path(args.outdir)

    outputs, records, cat_keys = build_outputs(source)
    ok, _ = self_check(source, records, cat_keys, outputs)

    if args.check:
        problems = check(outdir, outputs)
        for problem in problems:
            print("DRIFT: " + problem)
        if not ok:
            print("FAIL: the generator's own counts do not match the source")
        if problems or not ok:
            return 1
        print("check: %d files match, no drift" % len(outputs))
        return 0

    outdir.mkdir(parents=True, exist_ok=True)
    stale = {
        p.name for p in outdir.iterdir() if p.is_file()
    } - set(outputs)
    for name, text in sorted(outputs.items()):
        (outdir / name).write_text(text)
    print("wrote %d files to %s" % (len(outputs), outdir))
    for name in sorted(stale):
        print("note: %s is on disk but is no longer generated" % name)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
