"""export_literature_data.apply_contributions: a contributed correction is applied only after a confirmed page check."""
import json
import os
import shutil
import sys

import pytest

TOOLS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, TOOLS)
DATA = os.path.join(os.path.dirname(os.path.dirname(TOOLS)), "orchid-data")
CONTRIB = os.path.join(DATA, "contributions", "tools", "orchid_contrib.py")


@pytest.mark.skipif(not os.path.exists(CONTRIB), reason="orchid-data with contributions/ not found next to this repository")
def test_apply_contributions(tmp_path, monkeypatch):
    import export_literature_data as ex
    (tmp_path / "contributions" / "tools").mkdir(parents=True)
    shutil.copy(CONTRIB, tmp_path / "contributions" / "tools")
    (tmp_path / "contributions" / "corrections").mkdir()
    corr = {"record_id": "ORC-X1", "field": "temperature_K", "old": 590.15, "new": 591.15,
            "evidence": {"page": 2, "quote": "318 °C"}, "contributor": {"name": "A B"}, "licence_consent": "CC-BY-4.0"}
    (tmp_path / "contributions" / "corrections" / "a.jsonl").write_text(json.dumps(corr) + "\n")
    recs = [{"id": "ORC-X1", "temperature_K": 590.15, "flags": {}}]
    monkeypatch.setattr(ex.gates, "load_visual", lambda: {})
    assert ex.apply_contributions(recs, str(tmp_path)) == [] and recs[0]["temperature_K"] == 590.15
    monkeypatch.setattr(ex.gates, "load_visual", lambda: {"contributed_correction:ORC-X1:temperature_K": {"verdict": "confirmed", "page": 2}})
    assert ex.apply_contributions(recs, str(tmp_path)) == ["contributed_correction:ORC-X1:temperature_K"]
    assert recs[0]["temperature_K"] == 591.15 and recs[0]["corrections"][0]["before"] == 590.15
    assert recs[0]["contributors"][0]["name"] == "A B"
    assert ex.apply_contributions(recs, str(tmp_path / "nothing")) == []
