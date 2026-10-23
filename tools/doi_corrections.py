#!/usr/bin/env python3
"""Shared, versioned DOI correction table (tools/doi_corrections.json) and the lookups every DOI-emitting exporter uses.

The table lists printed references that an upstream matcher resolved to a DOI whose Crossref record describes a
different work (seeded from bibliography/DOI_MISMATCHES.md). Each entry names the context it applies to:

  literature_cited_source  match.literature_paper   the literature paper keyed by the DOI and every value pointing to it
  tdb_header               match.tdb                the link from one TDB file's header citation to the DOI
  handbook_source_line     match.handbook_system    the source line printed under one handbook diagram

`scope: context_only` drops the DOI only in that context (e.g. 10.2320/jinstmet.73.180 and 10.1016/j.tca.2014.01.002
are correct for their own papers and wrong only for one TDB header); `scope: everywhere` (a journal-level DOI such
as 10.1557/jmr) drops it in every context. Where a correction applies the DOI becomes null, `doi_status` =
"dropped_mismatch", `doi_correction` = the entry id, and the printed citation is kept.

The table is frozen and hand-maintained: no builder regenerates it.

  python3 tools/doi_corrections.py            # summary
  python3 tools/doi_corrections.py --check    # validate the table (exit 1 on a malformed entry)
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import sys

TABLE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "doi_corrections.json")
STATUS = "dropped_mismatch"
CONTEXTS = {"literature_cited_source": "literature_paper", "tdb_header": "tdb", "handbook_source_line": "handbook_system"}


def norm_doi(d):
    if not d or not isinstance(d, str):
        return None
    d = re.sub(r"^(https?://(dx\.)?doi\.org/|doi:)", "", d.strip(), flags=re.I)
    return d.lower() or None


def norm_system(s):
    if isinstance(s, (list, tuple)):
        s = "-".join(s)
    return "-".join(sorted(x.strip().upper() for x in str(s or "").replace("_", "-").split("-") if x.strip()))


class Corrections:
    def __init__(self, path=TABLE_PATH):
        self.path = path
        raw = open(path, "rb").read()
        self.sha256 = hashlib.sha256(raw).hexdigest()
        t = json.loads(raw)
        self.version = t["version"]
        self.entries = t["corrections"]
        self.by_id = {e["id"]: e for e in self.entries}
        self._by_doi = {}
        for e in self.entries:
            self._by_doi.setdefault(norm_doi(e["doi"]), []).append(e)

    # ---- lookups: the entry that drops `doi` in the given context, or None
    def _find(self, doi, context, key, value, norm=lambda v: v):
        es = self._by_doi.get(norm_doi(doi)) or []
        for e in es:
            if e["context"] == context and norm(e["match"].get(key)) == norm(value):
                return e
        for e in es:   # a DOI that names no article is dropped wherever it appears
            if e["scope"] == "everywhere":
                return e
        return None

    def literature(self, doi):
        """The literature paper keyed by this DOI (and every literature field pointing to it)."""
        return self._find(doi, "literature_cited_source", "literature_paper", "doi:" + (norm_doi(doi) or ""))

    def tdb(self, tdb, doi):
        """The DOI linked to TDB file `tdb` through its header citation."""
        return self._find(doi, "tdb_header", "tdb", os.path.basename(tdb or ""), norm=lambda v: (v or "").lower())

    def handbook(self, system, doi):
        """The DOI matched to the source line printed under the handbook diagram of `system`."""
        return self._find(doi, "handbook_source_line", "handbook_system", system, norm=norm_system)

    def anywhere(self, doi):
        """Only `everywhere` entries (e.g. a PDF fetched under a journal-level DOI)."""
        for e in self._by_doi.get(norm_doi(doi)) or []:
            if e["scope"] == "everywhere":
                return e
        return None

    def mark(self, obj, entry, doi_key="doi", status_key="doi_status"):
        """Null obj[doi_key] and record the correction (in place; returns obj)."""
        obj[doi_key] = None
        obj[status_key] = STATUS
        obj["doi_correction"] = entry["id"]
        return obj

    def provenance(self):
        return {"table": "tools/doi_corrections.json", "version": self.version, "sha256": self.sha256,
                "entries": len(self.entries)}


_DEFAULT = None


def load():
    global _DEFAULT
    if _DEFAULT is None:
        _DEFAULT = Corrections()
    return _DEFAULT


def validate(c):
    bad = []
    ids = set()
    for e in c.entries:
        if e["id"] in ids:
            bad.append("duplicate id " + e["id"])
        ids.add(e["id"])
        if e.get("context") not in CONTEXTS:
            bad.append("%s: unknown context %r" % (e["id"], e.get("context")))
        elif not (e.get("match") or {}).get(CONTEXTS[e["context"]]):
            bad.append("%s: match.%s missing" % (e["id"], CONTEXTS[e["context"]]))
        if e.get("scope") not in ("context_only", "everywhere"):
            bad.append("%s: scope %r" % (e["id"], e.get("scope")))
        if norm_doi(e.get("doi")) != e.get("doi"):
            bad.append("%s: DOI not normalised" % e["id"])
        if not e.get("printed_citation"):
            bad.append("%s: printed_citation missing" % e["id"])
        if e.get("status") != STATUS:
            bad.append("%s: status %r" % (e["id"], e.get("status")))
    return bad


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    c = load()
    bad = validate(c)
    if bad:
        print("\n".join(bad), file=sys.stderr)
        return 1
    if "--check" in argv:
        print("doi_corrections.json v%s OK (%d entries, sha256 %s)" % (c.version, len(c.entries), c.sha256[:12]))
        return 0
    for e in c.entries:
        print("%s  %-34s %-24s %-14s %s" % (e["id"], e["doi"], e["context"], e["scope"],
                                           e["match"][CONTEXTS[e["context"]]]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
