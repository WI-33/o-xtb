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
"""Generate the bash, zsh and fish completions of ``oxtb``.

The single source of truth is the option table of ``src/xhelp.f90`` - the
``call add(...)`` records and the named legal-value arrays they name with
``values=v_...``.  Nothing here is a second option table:

  * ``test/unit/landscape_help_check.py`` is the ``--help`` coverage gate.  Its
    ``HELP_LINE`` quoting rule, ``SPLIT``/``TOKEN`` tokeniser and
    ``documented_flags`` spelling set are the zero-drift oracle this generator
    is checked against, so the completions, the help page and the man page
    cannot disagree about what a spelling is;
  * ``docs/docs-site/tools/gen_cli_reference.py``'s ``parse_records`` is the
    same reader the documentation site uses (that module imports the gate
    too), so the completions and the site reference cannot disagree about what
    a legal value is either.

Three completion classes are derived per record:

  values    the record names a value set (``values=v_...``) - complete those
            words, nothing else;
  files     the record's placeholder is ``FILE``/``[FILE]``/``[EXECUTABLE]``
            or a spelling ends in ``-file`` - complete file names;
  argument  any other argument (``INT``, ``REAL``, ``VALUE``, ``N``, ...) -
            complete nothing, so a numeric option does not offer files;
  flag      no argument at all - complete nothing.

A spelling is offered as an option only when the parser would see it as a
flag (``src/prog/argparser.f90``, ``initArgument``: starts with a dash,
nothing but ``[A-Za-z0-9_-]`` after it, and not number-like).  The numeric
level aliases ``-1``/``-2``/``-3`` of the ``levels`` help text are therefore
values of ``--opt`` and never completion candidates.

Only the space separated value syntax is offered: measured, ``--gfn 1``
selects GFN1-xTB while ``--gfn=1`` is neither a flag nor a file name for
``initArgument`` and is silently ignored (GFN2 default) - see report
``167-shell-completions.md``.

Usage::

    gen_completions.py [--check] [--source PATH] [--outdir DIR]

``--check`` writes nothing and exits non-zero when an artifact on disk differs
from what this run would produce (the drift gate ``xtb/landscape-completions``).
The output is deterministic: no timestamp, no path and no set iteration order
reaches a file, so two runs are byte identical.
"""

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
REPO_ROOT = HERE.parents[2]
DEFAULT_SOURCE = REPO_ROOT / "src" / "xhelp.f90"
DEFAULT_OUTDIR = REPO_ROOT / "completions"

sys.path.insert(0, str(REPO_ROOT / "test" / "unit"))
sys.path.insert(0, str(REPO_ROOT / "docs" / "docs-site" / "tools"))
from landscape_help_check import SPLIT, TOKEN, documented_flags  # noqa: E402
from gen_cli_reference import parse_records, spelling_tokens  # noqa: E402

#: The counts this generator is built on.  The gate prints the same numbers, so
#: an option added to or removed from `src/xhelp.f90` has to be reviewed here as
#: well: the counts fail loudly instead of shipping a silently stale script.
#: `EXPECTED_SPELLINGS` is one below the gate's documented set (247) by design:
#: the numeric level aliases `-1`/`-2`/`-3` are argument values, not flags.
EXPECTED_RECORDS = 228
EXPECTED_SPELLINGS = 249

#: Placeholders that mean "a path or a program".
FILE_PLACEHOLDERS = ("FILE", "[FILE]", "EXECUTABLE", "[EXECUTABLE]")
#: Characters the parser accepts after the leading dashes of a flag.
FLAG_CHARS = set(
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-"
)
NUMBER_CHARS = set("0123456789.eE+-")

GENERATED = (
    "generated from src/xhelp.f90 by completions/tools/gen_completions.py, "
    "do not edit"
)

CLASS_VALUES = "values"
CLASS_FILES = "files"
CLASS_ARGUMENT = "argument"
CLASS_FLAG = "flag"


def is_parser_flag(token):
    """True when `src/prog/argparser.f90`'s `initArgument` would read it as a
    flag: a leading dash, only `[A-Za-z0-9_-]` afterwards and not number-like."""
    if TOKEN.match(token) is None or not token.startswith("-"):
        return False
    if set(token) - FLAG_CHARS:
        return False
    return bool(set(token) - NUMBER_CHARS)


def option_entries(source):
    """Read the option table into `[(category, spellings, summary, values, kind,
    has_optional_argument)]` in `--help` order.

    A record's spellings are expanded exactly like the gate does it (a merged
    `--[no]x` becomes `--x` and `--nox`); a record the parser would not see as
    a flag is dropped by `is_parser_flag`.
    """
    entries = []
    dropped = []
    for spellings, category, summary, _detail, values in parse_records(source):
        tokens = [t for t in SPLIT.split(spellings.strip()) if t]
        placeholders = [t for t in tokens if not TOKEN.match(t)]
        flags = [w for w in spelling_tokens(spellings) if is_parser_flag(w)]
        if not flags:
            dropped.append(spellings.strip())
            continue
        kind = CLASS_FLAG
        if values:
            kind = CLASS_VALUES
        elif any(p in FILE_PLACEHOLDERS for p in placeholders) or any(
            f.endswith("-file") for f in flags
        ):
            kind = CLASS_FILES
        elif placeholders:
            kind = CLASS_ARGUMENT
        entries.append(
            (
                category,
                flags,
                summary,
                list(values),
                kind,
                any(p.startswith("[") and p.endswith("]") for p in placeholders),
            )
        )
    return entries, dropped


def all_spellings(entries):
    """The candidate option spellings, in `--help` order, without duplicates."""
    seen = set()
    ordered = []
    for _category, flags, _summary, _values, _kind, _optional in entries:
        for flag in flags:
            if flag not in seen:
                seen.add(flag)
                ordered.append(flag)
    return ordered


def self_check(source, entries, dropped):
    """Assert the invariants this generator is built on and print the counts."""
    offered = all_spellings(entries)
    by_gate = documented_flags(source, expand_no=True)
    excluded = sorted(by_gate - set(offered))
    # Every spelling the gate sees but this generator does not offer has to be
    # one the parser does not read as a flag - the numeric level aliases.
    bad = [token for token in excluded if is_parser_flag(token)]
    unknown = sorted(set(offered) - by_gate)
    print(
        "records (call add):        %d  (expected %d, %d without a flag)"
        % (len(entries) + len(dropped), EXPECTED_RECORDS, len(dropped))
    )
    print(
        "distinct spellings:        %d  (expected %d)"
        % (len(offered), EXPECTED_SPELLINGS)
    )
    print(
        "valued spellings:          %d"
        % sum(len(f) for _c, f, _s, _v, k, _o in entries if k == CLASS_VALUES)
    )
    print(
        "file spellings:            %d"
        % sum(len(f) for _c, f, _s, _v, k, _o in entries if k == CLASS_FILES)
    )
    print(
        "records without a flag:    %s"
        % ("; ".join(dropped) if dropped else "(none)")
    )
    print(
        "not offered (not a flag):  %s"
        % (" ".join(excluded) if excluded else "(none)")
    )
    print(
        "offered but not in the gate's documented set: %s"
        % (" ".join(unknown) if unknown else "(none)")
    )
    ok = (
        len(entries) + len(dropped) == EXPECTED_RECORDS
        and len(offered) == EXPECTED_SPELLINGS
        and not bad
        and not unknown
    )
    return ok


# --- rendering ------------------------------------------------------------


def zsh_quote(text):
    """A single-segment zsh spec from `text`: quote for the shell and escape the
    characters `_arguments` reads (`[`, `]`, `:`) inside a description."""
    escaped = ""
    for char in text:
        if char in "\\[]:":
            escaped += "\\" + char
        else:
            escaped += char
    return escaped.replace("'", "'\\''")


def fish_quote(text):
    """A single-quoted fish word."""
    return "'" + text.replace("\\", "\\\\").replace("'", "\\'") + "'"


def bash_arms(entries, kind):
    """`pattern|pattern)` arms of one completion class, in `--help` order."""
    arms = []
    for _category, flags, _summary, values, entry_kind, _optional in entries:
        if entry_kind != kind:
            continue
        arms.append((flags, values))
    return arms


def render_bash(entries, spellings):
    out = [
        "# bash completion for oxtb(1).",
        "# " + GENERATED + ".",
        "#",
        "# Source of truth: src/xhelp.f90.  Install per user:",
        "#   source ~/.local/share/bash-completion/completions/oxtb",
        "# The file is self contained - it uses only the `compgen` and `compopt`",
        "# builtins, so the bash-completion package (and its _init_completion) is",
        "# not required.  Only the space separated value syntax is offered:",
        "# `--gfn 1` is parsed, `--gfn=1` is silently ignored by the program.",
        "",
        "_oxtb_spellings='",
    ]
    out.extend("  " + spelling for spelling in spellings)
    out.append("'")
    out.append("")
    out.append("# file names: the geometry argument and the file taking options")
    out.append("_oxtb_files()")
    out.append("{")
    out.append("    COMPREPLY=( $(compgen -f -- \"$1\") )")
    out.append("    if [ ${#COMPREPLY[@]} -gt 0 ] && type compopt >/dev/null 2>&1; then")
    out.append("        compopt -o filenames 2>/dev/null")
    out.append("    fi")
    out.append("}")
    out.append("")
    out.append("_oxtb()")
    out.append("{")
    out.append("    local cur prev")
    out.append("    cur=${COMP_WORDS[COMP_CWORD]}")
    out.append("    prev=''")
    out.append("    if [ \"${COMP_CWORD}\" -gt 0 ]; then")
    out.append("        prev=${COMP_WORDS[COMP_CWORD-1]}")
    out.append("    fi")
    out.append("    case $prev in")
    for flags, values in bash_arms(entries, CLASS_VALUES):
        out.append("    %s)" % "|".join(flags))
        out.append("        COMPREPLY=( $(compgen -W '%s' -- \"$cur\") )" % " ".join(values))
        out.append("        return 0")
        out.append("        ;;")
    for flags, _values in bash_arms(entries, CLASS_FILES):
        out.append("    %s)" % "|".join(flags))
        out.append("        _oxtb_files \"$cur\"")
        out.append("        return 0")
        out.append("        ;;")
    argument_arms = [flags for flags, _v in bash_arms(entries, CLASS_ARGUMENT)]
    if argument_arms:
        out.append("    %s)" % "|".join(flag for flags in argument_arms for flag in flags))
        out.append("        # an argument without a value set: no candidates")
        out.append("        COMPREPLY=()")
        out.append("        return 0")
        out.append("        ;;")
    out.append("    esac")
    out.append("    if [[ $cur == -* ]]; then")
    out.append("        COMPREPLY=( $(compgen -W \"$_oxtb_spellings\" -- \"$cur\") )")
    out.append("        return 0")
    out.append("    fi")
    out.append("    _oxtb_files \"$cur\"")
    out.append("}")
    out.append("")
    out.append("complete -F _oxtb oxtb")
    out.append("")
    return "\n".join(out)


def zsh_specs(entries):
    """One `_arguments` spec per spelling, in `--help` order."""
    specs = []
    for _category, flags, summary, values, kind, optional in entries:
        description = zsh_quote(summary)
        for flag in flags:
            if kind == CLASS_VALUES:
                action = ("::" if optional else ":") + "values:(%s)" % (
                    " ".join(values),
                )
            elif kind == CLASS_FILES:
                action = "::file:_files" if optional else ":file:_files"
            elif kind == CLASS_ARGUMENT:
                action = "::argument:" if optional else ":argument:"
            else:
                action = ""
            specs.append("        '%s[%s]%s'" % (flag, description, action))
    return specs


def render_zsh(entries):
    out = [
        "#compdef oxtb",
        "# zsh completion for oxtb(1).",
        "# " + GENERATED + ".",
        "#",
        "# Source of truth: src/xhelp.f90.  Install per user by putting",
        "#   fpath+=(~/.local/share/zsh/site-functions)",
        "# into ~/.zshrc *before* the compinit call.",
        "# Only the space separated value syntax is offered: `--gfn 1` is parsed,",
        "# `--gfn=1` is silently ignored by the program.",
        "",
        "_oxtb()",
        "{",
        "    local -a _oxtb_specs",
        "    _oxtb_specs=(",
    ]
    out.extend(zsh_specs(entries))
    out.append("    )")
    out.append("    _arguments -s -S \"${_oxtb_specs[@]}\" '*:file:_files'")
    out.append("}")
    out.append("")
    out.append("_oxtb \"$@\"")
    out.append("")
    return "\n".join(out)


def render_fish(entries):
    out = [
        "# fish completion for oxtb(1).",
        "# " + GENERATED + ".",
        "#",
        "# Source of truth: src/xhelp.f90.  Install per user:",
        "#   ~/.config/fish/completions/oxtb.fish",
        "# Only the space separated value syntax is offered: `--gfn 1` is parsed,",
        "# `--gfn=1` is silently ignored by the program.",
        "",
    ]
    for _category, flags, summary, values, kind, _optional in entries:
        for flag in flags:
            if flag.startswith("--"):
                head = "-l " + flag[2:]
            else:
                head = "-s " + flag[1:]
            line = "complete -c oxtb " + head
            if kind == CLASS_VALUES:
                line += " -x -a '%s'" % " ".join(values)
            elif kind == CLASS_FILES:
                line += " -r -F"
            elif kind == CLASS_ARGUMENT:
                line += " -x"
            line += " -d " + fish_quote(summary)
            out.append(line)
    out.append("")
    return "\n".join(out)


def build_outputs(source):
    """`{file name: text}` for the three scripts - a pure function of `source`."""
    entries, dropped = option_entries(source)
    spellings = all_spellings(entries)
    return (
        {
            "oxtb.bash": render_bash(entries, spellings),
            "_oxtb": render_zsh(entries),
            "oxtb.fish": render_fish(entries),
        },
        entries,
        dropped,
    )


def check(outdir, outputs):
    """Compare `outputs` with `outdir` on disk; print and return the drift."""
    problems = []
    for name, text in sorted(outputs.items()):
        path = outdir / name
        if not path.is_file():
            problems.append("%s is missing" % path)
        elif path.read_text() != text:
            problems.append("%s differs from the generated output" % path)
    return problems


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="write nothing; exit non-zero when an artifact is stale",
    )
    parser.add_argument("--source", default=str(DEFAULT_SOURCE))
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR))
    args = parser.parse_args(argv)

    source = Path(args.source).read_text()
    outdir = Path(args.outdir)

    outputs, entries, dropped = build_outputs(source)
    ok = self_check(source, entries, dropped)

    if args.check:
        problems = check(outdir, outputs)
        for problem in problems:
            print("DRIFT: " + problem)
        if problems or not ok:
            if not ok:
                print("FAIL: the generator's own counts do not match the source")
            return 1
        print("check: %d files match, no drift" % len(outputs))
        return 0

    outdir.mkdir(parents=True, exist_ok=True)
    for name, text in sorted(outputs.items()):
        (outdir / name).write_text(text)
    print("wrote %d files to %s" % (len(outputs), outdir))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
