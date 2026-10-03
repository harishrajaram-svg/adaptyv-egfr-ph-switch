#!/usr/bin/env python3
"""ESMFold2 CIF -> PDB, with the atom-name column rule enforced.

Why this file exists (2026-10-02): the converter was written inline, from memory,
four separate times. The fourth copy left-justified atom names into column 13.
Coordinates were untouched, so contact counts matched to the atom and the output
looked right -- but propka types atoms from columns 13-16, so every pKa it
returned was garbage. That produced a clean, plausible, completely false result
(a design "passing" the pH gate 10/10) that was only caught by byte-diffing
against an earlier known-good file.

Two rules, both enforced here rather than remembered:

1. PDB columns 13-16 hold the atom name. A one-character element symbol is
   right-justified in 13-14, so the name starts in column 14 (" CA "). A
   two-character element (FE, CL) starts in column 13 ("FE  "). Proteins are
   all C/N/O/S, so in practice every name needs a leading space. Getting this
   wrong does not fail loudly -- it silently changes atom typing downstream.

2. Column indices in the CIF are read from the _atom_site loop header by tag
   name, never hardcoded. The wrapper's column order is stable today; it is not
   a promise.

Usage:
    cif2pdb.py <in.cif> <out.pdb>
    cif2pdb.py --self-test
"""
import sys
from pathlib import Path

NEEDED = ("type_symbol", "label_atom_id", "label_comp_id", "auth_seq_id",
          "auth_asym_id", "B_iso_or_equiv", "occupancy",
          "Cartn_x", "Cartn_y", "Cartn_z")


def _atom_site_columns(lines):
    """Map _atom_site tag -> field index, read from the loop header."""
    cols = {}
    for line in lines:
        s = line.strip()
        if s.startswith("_atom_site."):
            cols[s.split(".", 1)[1].split()[0]] = len(cols)
        elif cols and (s.startswith("ATOM") or s.startswith("HETATM")):
            break
    missing = [t for t in NEEDED if t not in cols]
    if missing:
        raise SystemExit(f"CIF _atom_site loop is missing: {', '.join(missing)}")
    return cols


def pdb_atom_name(name, element):
    """Place an atom name in PDB columns 13-16. See rule 1 in the module docstring."""
    name = name.strip().strip('"').strip("'")
    if len(element.strip()) == 2:       # FE, CL, ZN ... start at column 13
        return f"{name:<4s}"
    return f" {name:<3s}"               # C, N, O, S ... start at column 14


def convert(src, dst):
    lines = Path(src).read_text().splitlines()
    c = _atom_site_columns(lines)
    out, n = [], 0
    for line in lines:
        if not line.startswith("ATOM"):
            continue
        f = line.split()
        if len(f) <= max(c.values()):
            continue
        el = f[c["type_symbol"]]
        n += 1
        out.append(
            f"ATOM  {n:5d} {pdb_atom_name(f[c['label_atom_id']], el)}"
            f"{f[c['label_comp_id']]:>4s} {f[c['auth_asym_id']]}"
            f"{int(f[c['auth_seq_id']]):4d}    "
            f"{float(f[c['Cartn_x']]):8.3f}{float(f[c['Cartn_y']]):8.3f}"
            f"{float(f[c['Cartn_z']]):8.3f}"
            f"{float(f[c['occupancy']]):6.2f}{float(f[c['B_iso_or_equiv']]):6.2f}"
            f"          {el:>2s}\n"
        )
    if not out:
        raise SystemExit(f"no ATOM records produced from {src}")
    validate(out, src)
    out.append("END\n")
    Path(dst).write_text("".join(out))
    return n


def validate(atom_lines, src=""):
    """Reject the exact failure that made this file necessary."""
    for line in atom_lines:
        el = line[76:78].strip()
        if len(el) == 1 and line[12] != " ":
            raise SystemExit(
                f"REFUSING TO WRITE {src}: atom name occupies column 13 with a "
                f"one-character element ({el!r}). propka would mis-type every "
                f"atom and return plausible, wrong pKa values.\n  {line.rstrip()}"
            )


def self_test():
    # The bug: a 1-char element whose name was left-justified into column 13.
    bad = ("ATOM      1 N    LEU A   1      35.453  -2.419   0.407"
           "  1.00 43.77           N\n")
    good = ("ATOM      1  N   LEU A   1      35.453  -2.419   0.407"
            "  1.00 43.77           N\n")
    assert pdb_atom_name("N", "N") == " N  ", repr(pdb_atom_name("N", "N"))
    assert pdb_atom_name("CA", "C") == " CA ", repr(pdb_atom_name("CA", "C"))
    assert pdb_atom_name("CD1", "C") == " CD1", repr(pdb_atom_name("CD1", "C"))
    assert pdb_atom_name("FE", "FE") == "FE  ", repr(pdb_atom_name("FE", "FE"))
    validate([good])                                  # must pass
    try:
        validate([bad])
    except SystemExit:
        pass
    else:
        raise AssertionError("validate() failed to catch the column-13 bug")
    assert good[12:16] == " N  " and good[17:20] == "LEU" and good[21] == "A"
    print("self-test OK: column rule enforced, column-13 bug is caught")


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        self_test()
    elif len(sys.argv) == 3:
        print(convert(sys.argv[1], sys.argv[2]), sys.argv[2])
    else:
        raise SystemExit(__doc__)
