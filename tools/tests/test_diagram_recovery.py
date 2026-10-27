import sys
from pathlib import Path
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compute_binary_diagrams import prepare_source, validate_geometry, source_phase_selection


def test_legacy_encoding_preserves_parameter_text():
    raw = b"$ H\xe4m\xe4l\xe4inen\nPARAMETER G(LIQUID,CU;0) 298.15 +123.45; 6000 N !"
    text, encoding, removed = prepare_source(raw)
    assert text.encode("latin-1") == raw
    assert encoding == "latin-1" and not removed


def test_only_unused_dangling_order_directive_removed():
    directive = "TYPE-DEF O GES AMEND_PHASE_DESCRIPTION B2_BCC DIS_PART BCC_A2 !"
    base = "PHASE BCC_A2 %B 2 1 3 !\n" + directive
    text, _, changes = prepare_source(base.encode())
    assert directive not in text and len(changes) == 1
    for protected in (base.replace("%B", "%BO"), base + "\nPHASE B2_BCC % 2 1 1 !"):
        assert prepare_source(protected.encode())[0] == protected
    assert prepare_source(("$ " + directive).encode())[0] == "$ " + directive


@pytest.mark.parametrize("x,t", [(float("nan"), 500), (.2,float("inf")), (1.1,500), (.2,100)])
def test_invalid_geometry_rejected(x,t):
    with pytest.raises(ValueError):
        validate_geometry({"T_range":[300,2500],"regions":[[{"x":[x],"T":[t]}]],"invariants":[]})


def test_valid_geometry_and_mismatched_lengths():
    rec = {"T_range":[300,2500],"regions":[[{"x":[0,.2],"T":[300,500]}]],"invariants":[]}
    validate_geometry(rec)
    rec["regions"][0][0]["T"].pop()
    with pytest.raises(ValueError): validate_geometry(rec)


@pytest.mark.parametrize("regions", [[], [[{"x": [], "T": []}]]])
def test_empty_mapping_is_not_success(regions):
    with pytest.raises(ValueError, match="No mapped phase boundaries"):
        validate_geometry({"T_range": [300, 2500], "regions": regions, "invariants": []})


def test_source_phase_selection_uses_explicit_defaults_only():
    source = "$ DEFAULT_COM REJECT_PHASE O2GAS !\nDEFAULT-COM REJECT_PHASE GAS !"
    phases, evidence = source_phase_selection(source, ["GAS", "O2GAS", "LIQUID"])
    assert phases == ["LIQUID", "O2GAS"]
    assert evidence["excluded"] == ["GAS"]
    assert source_phase_selection("", ["GAS", "O2GAS"])[0] == ["GAS", "O2GAS"]


def test_unknown_source_phase_fails_closed():
    with pytest.raises(ValueError, match="Unknown phase"):
        source_phase_selection("DEFAULT_COM REJECT_PHASE GSA !", ["GAS"])


def test_mgy_repair_preserves_every_other_byte_and_is_hash_scoped():
    import hashlib
    source = Path(__file__).resolve().parents[3] / "open-tdb/databases/sgte-binary/MgY-98Luk.tdb"
    raw = source.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == "68740131a132fee3651ec3eb1fd044d3a56428ed2828054a4818d2079e1115cc"
    text, encoding, changes = prepare_source(raw)
    old = " CONST MGY_B2 : MG : MG Y% : VA : !"
    new = " CONST MGY_B2 : MG : MG Y% : !"
    assert text.encode(encoding) == raw.replace(old.encode(), new.encode(), 1)
    assert changes == [{"rule": "verified_mgy_extra_vacancy_sublattice", "statement": old.strip(),
                        "replacement": new.strip(), "evidence": "review/mgy-1998-20261026/provenance.json"}]
    modified = raw + b"\n$ Different source version\n"
    assert old in prepare_source(modified)[0]
