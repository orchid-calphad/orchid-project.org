"""report_issue_url: pre-filled issue links carry only the record's non-personal fields (same vectors as orchid-mcp)."""
import os
import sys
from urllib.parse import parse_qs, urlsplit

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from orchid_report import report_issue_url  # noqa: E402

# shared test vectors: orchid-mcp tests/test_report.py has the same three, with the same expected URLs
VECTORS = [
    (("ORC-RZXYDDJGWB", "AG-AL", "10.1016/j.jallcom.2004.04.126", 2, "590.15 K"),
     "https://github.com/orchid-calphad/orchid-data/issues/new?template=report-value-error.yml"
     "&title=%5Bvalue%5D%20ORC-RZXYDDJGWB%20%28AG-AL%29&record_id=ORC-RZXYDDJGWB&system=AG-AL"
     "&source=doi%3A10.1016%2Fj.jallcom.2004.04.126%2C%20p.%202&value_shown=590.15%20K"),
    (("ORC-5D324PIVW6.310d73", "AL-CR", "https://doi.org/10.1007/bf02735058", "5", "1.6386 nm"),
     "https://github.com/orchid-calphad/orchid-data/issues/new?template=report-value-error.yml"
     "&title=%5Bvalue%5D%20ORC-5D324PIVW6.310d73%20%28AL-CR%29&record_id=ORC-5D324PIVW6.310d73&system=AL-CR"
     "&source=doi%3A10.1007%2Fbf02735058%2C%20p.%205&value_shown=1.6386%20nm"),
    (("ORC-ABCDEFGHIJ", None, None, None, None),
     "https://github.com/orchid-calphad/orchid-data/issues/new?template=report-value-error.yml"
     "&title=%5Bvalue%5D%20ORC-ABCDEFGHIJ&record_id=ORC-ABCDEFGHIJ"),
]


def test_vectors():
    for args, url in VECTORS:
        assert report_issue_url(*args) == url


def test_only_whitelisted_fields_and_validation():
    u = report_issue_url("ORC-RZXYDDJGWB", "ag-al; drop", "not a doi", "p2<script>", "x" * 500)
    q = parse_qs(urlsplit(u).query)
    assert set(q) == {"template", "title", "record_id", "value_shown"}
    assert len(q["value_shown"][0]) == 120 and q["value_shown"][0].endswith("…")
    assert report_issue_url("not-an-id") is None and report_issue_url(None) is None
    assert report_issue_url("ORC-<b>bad</b>") is None


def test_round_trip_and_no_personal_parameters():
    u = report_issue_url("ORC-RZXYDDJGWB", "AG-AL", "10.1016/j.jallcom.2004.04.126", 2, "−5.2 ± 0.3 kJ/mol-atom")
    q = {k: v[0] for k, v in parse_qs(urlsplit(u).query).items()}
    assert q["value_shown"] == "−5.2 ± 0.3 kJ/mol-atom"
    assert q["source"] == "doi:10.1016/j.jallcom.2004.04.126, p. 2"
    assert not set(q) & {"name", "email", "orcid", "credit_name", "user", "assignees"}
    assert " " not in u and "&amp;" not in u
