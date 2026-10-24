#!/usr/bin/env python3
"""Tests of the annotation pipeline's phase_mapper.py structure rules (2026-10-24, James Saal review): cF8 alone no longer
names DIAMOND_A4; rock-salt (B1) carbides map to the TDB's own NaCl-type phase or to FCC_A1 (M)(C,VA); sphalerite SiC to
SIC_B3; FeSi-type silicides never to DIAMOND_A4; a Pearson-only standard name never overrides a printed Strukturbericht.
Reads the SGTE TDB files of open-tdb (read only).
usage: ~/orchid-py313/bin/python tools/tests/test_phase_mapper_structures.py [--tools DIR] [--tdb-dir DIR]"""
import argparse
import os
import sys

sys.dont_write_bytecode = True
HOME = os.path.expanduser("~")
ap = argparse.ArgumentParser()
ap.add_argument("--tools", default=os.path.join(HOME, "projects", "open-tdb", "build", "okamoto-pilot", "tools"))
ap.add_argument("--tdb-dir", default=os.path.join(HOME, "projects", "open-tdb", "databases", "sgte-binary"))
A = ap.parse_args()
sys.path.insert(0, A.tools)
import phase_mapper as PM  # noqa: E402

fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


def m(tdb_file, label, info, axis):
    p = os.path.join(A.tdb_dir, tdb_file)
    tdb = PM.tdb_phases(p)
    comps = PM.tdb_compositions(p, list(axis))
    return PM.map_phase(label, info, tdb, comps, list(axis))[0]


B1 = {"pearson_symbol": "cF8", "strukturbericht": "B1", "prototype": "NaCl"}
B3 = {"pearson_symbol": "cF8", "strukturbericht": "B3", "prototype": "ZnS (sphalerite)"}
B20 = {"pearson_symbol": "cF8", "strukturbericht": "B20", "prototype": "FeSi"}   # as the handbook table was read
A4 = {"pearson_symbol": "cF8", "strukturbericht": "A4", "prototype": "C (diamond)"}

check("cF8 no longer a standard compound key", "cF8" not in PM.STD_COMPOUND)
for f, lab, ax in (("HfC-16Pan.tdb", "HfC", ("HF", "C")), ("NbC-01Lee.tdb", "NbC", ("NB", "C")),
                   ("TiC-03Fri.tdb", "CTi", ("TI", "C")), ("ZrC-95Fer-LB.tdb", "ZrC", ("ZR", "C")),
                   ("TaC-16Pan.tdb", "TaC", ("TA", "C")), ("VC-06Hu.tdb", "VC", ("V", "C")),
                   ("WC-86Gus.tdb", "βWC", ("W", "C"))):
    check("B1 carbide %s in %s -> FCC_A1" % (lab, f), m(f, lab, B1, ax) == "FCC_A1")
check("B1 PuC with its own PUC_B1 phase -> PUC_B1", m("PuC-08Fis.tdb", "PuC", B1, ("PU", "C")) == "PUC_B1")
check("B1 SnTe -> SNTE_B1", m("SnTe-10Liu.tdb", "SnTe", B1, ("SN", "TE")) == "SNTE_B1")
check("B1 AsSn -> SNAS_B1", m("SnAs-21Wu.tdb", "AsSn", B1, ("AS", "SN")) == "SNAS_B1")
check("B1 letter-named UC (delta, 47-50 at.% C) -> FCC_B1",
      m("UC-01Che.tdb", "δ", dict(B1, composition_at_pct=[47.0, 50.0]), ("U", "C")) == "FCC_B1")
check("B3 SiC -> SIC_B3", m("SiC-17Lia.tdb", "SiC", B3, ("C", "SI")) == "SIC_B3")
check("B20 CoSi -> COSI_B20", m("CoSi-06Zha.tdb", "CoSi", B20, ("CO", "SI")) == "COSI_B20")
check("B20 CrSi never DIAMOND_A4", m("CrSi-00Du.tdb", "CrSi", B20, ("CR", "SI")) != "DIAMOND_A4")
check("B20 MnSi never DIAMOND_A4", m("MnSi-14Ber.tdb", "MnSi", B20, ("MN", "SI")) != "DIAMOND_A4")
check("solution (Si) A4 still DIAMOND_A4", m("CuSi-00Yan.tdb", "(Si)", A4, ("CU", "SI")) == "DIAMOND_A4")
check("solution (Ge) A4 still DIAMOND_A4", m("AgGe-11Wan.tdb", "(Ge)", A4, ("AG", "GE")) == "DIAMOND_A4")

# synthetic TDB dictionaries (no file)
ox = {"FCC_A1": {"pearson": "cF4", "constituents": [["NI"], ["O", "VA"]]}, "LIQUID": {}}
check("B1 oxide NiO is not FCC_A1 (only C/N interstitial carbides/nitrides)",
      PM.map_phase("NiO", B1, ox, None, ["NI", "O"])[0] is None)
hal = dict(ox, HALITE={"pearson": "cF8", "prototype": "NaCl"})
check("B1 oxide with a HALITE phase -> HALITE", PM.map_phase("NiO", B1, hal, None, ["NI", "O"])[0] == "HALITE")
nit = {"FCC_A1": {"constituents": [["TI"], ["N", "VA"]]}}
check("B1 nitride TiN -> FCC_A1 (Ti)(N,VA)", PM.map_phase("TiN", B1, nit, None, ["TI", "N"])[0] == "FCC_A1")
dh = {"DHCP": {"pearson": "hP4"}}
check("Pearson-only hP4 does not override a printed B8_1 (NiAs type)",
      PM.map_phase("NiSb", {"pearson_symbol": "hP4", "strukturbericht": "B81"}, dh, None, ["NI", "SB"])[0] != "DHCP")
check("printed A3' with hP4 still DHCP", PM.map_phase("αLa", {"pearson_symbol": "hP4", "strukturbericht": "A3'"},
                                                       dh, None, ["LA", "X"])[0] == "DHCP")

print("%d failed" % len(fails))
sys.exit(1 if fails else 0)
