#!/usr/bin/env python3
"""Pre-filled "Report an error in a value" GitHub issue URLs for ORCHID records.

The URL opens the issue form `.github/ISSUE_TEMPLATE/report-value-error.yml` of orchid-calphad/orchid-data with
the record's non-personal fields filled in: record id, chemical system, source (DOI and page) and the value as
shown. Nothing about the person reporting is ever put in the URL, and nothing is posted: the reader submits the
form on GitHub. The same function (same output for the same input) is in orchid-mcp (`orchid_mcp.report`).

Standard library only; deterministic.
"""
from __future__ import annotations

import re
from urllib.parse import quote, urlencode

ISSUE_URL = "https://github.com/orchid-calphad/orchid-data/issues/new"
TEMPLATE = "report-value-error.yml"
VALUE_MAX = 120
_ID = re.compile(r"^ORC-[A-Za-z0-9.\-]{4,40}$")
_SYS = re.compile(r"^[A-Z]{1,2}(?:-[A-Z]{1,2}){0,7}$")
_DOI = re.compile(r"^10\.\d{4,9}/\S{1,200}$")
_PAGE = re.compile(r"^[A-Za-z]{0,3}\d{1,5}$")


def _clean(s, n):
    s = " ".join(str(s).split())
    return s if len(s) <= n else s[: n - 1].rstrip() + "…"


def report_issue_url(record_id, system=None, doi=None, page=None, value=None):
    """URL of the pre-filled issue form for one value, or None when record_id is not an ORCHID id.

    Only these fields go into the URL, each validated: record_id (ORC-...), system (CU-ZR), doi (10.xxxx/...),
    page (12, S3) and value (the value as shown, at most 120 characters).
    """
    rid = str(record_id or "").strip()
    if not _ID.match(rid):
        return None
    sysid = str(system or "").strip().upper()
    sysid = sysid if _SYS.match(sysid) else None
    d = re.sub(r"(?i)^(https?://(dx\.)?doi\.org/|doi:)", "", str(doi or "").strip())
    d = d if _DOI.match(d) else None
    pg = str(page).strip() if page is not None and _PAGE.match(str(page).strip()) else None
    params = [("template", TEMPLATE), ("title", "[value] %s%s" % (rid, " (%s)" % sysid if sysid else "")),
              ("record_id", rid)]
    if sysid:
        params.append(("system", sysid))
    if d or pg:
        params.append(("source", ", ".join(x for x in ("doi:%s" % d if d else None, "p. %s" % pg if pg else None) if x)))
    if value not in (None, ""):
        params.append(("value_shown", _clean(value, VALUE_MAX)))
    return ISSUE_URL + "?" + urlencode(params, quote_via=quote, safe="")
