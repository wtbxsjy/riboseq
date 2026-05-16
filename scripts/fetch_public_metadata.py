#!/usr/bin/env python3
"""Fetch public sequencing metadata from ENA and NCBI.

This helper expands public accessions to run-level metadata, infers likely
assay type (`riboseq`, `rnaseq`, `tiseq`, `unknown`), suggests grouping hints,
and optionally emits a download manifest and a samplesheet template.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


ENA_FIELDS = [
    "run_accession",
    "study_accession",
    "secondary_study_accession",
    "experiment_accession",
    "sample_accession",
    "secondary_sample_accession",
    "scientific_name",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_layout",
    "fastq_ftp",
    "fastq_md5",
    "submitted_ftp",
    "submitted_format",
    "study_title",
    "experiment_title",
    "sample_title",
    "instrument_platform",
    "instrument_model",
]

CURATED_COLUMNS = [
    "input_accession",
    "run_accession",
    "study_accession",
    "project_accession",
    "experiment_accession",
    "sample_accession",
    "biosample_accession",
    "sample_title",
    "experiment_title",
    "study_title",
    "scientific_name",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_layout",
    "platform",
    "model",
    "fastq_ftp_1",
    "fastq_ftp_2",
    "fastq_md5_1",
    "fastq_md5_2",
    "sra_download_accession",
    "inferred_type",
    "type_confidence",
    "type_evidence",
    "suggested_group",
    "group_evidence",
    "replicate_hint",
    "needs_manual_review",
]

DOWNLOAD_COLUMNS = [
    "input_accession",
    "run_accession",
    "file_role",
    "file_name",
    "ftp_url",
    "md5",
    "recommended_method",
    "fallback_sra_accession",
]

SAMPLESHEET_COLUMNS = [
    "sample",
    "fastq_1",
    "fastq_2",
    "strandedness",
    "type",
    "group",
    "run_accession",
    "input_accession",
]

TYPE_PATTERNS = {
    "riboseq": [
        r"\bribo[\s\-_]?seq\b",
        r"\bribosome profiling\b",
        r"\bribosome protected fragments?\b",
        r"\brpf\b",
        r"\bfootprint(s|ing)?\b",
    ],
    "rnaseq": [
        r"\brna[\s\-_]?seq\b",
        r"\btranscriptome\b",
        r"\bmrna\b",
    ],
    "tiseq": [
        r"\bti[\s\-_]?seq\b",
        r"\btranslation initiation\b",
    ],
}

GROUP_HINT_PATTERNS = [
    ("control", [r"\bcontrol\b", r"\buntreated\b", r"\bvehicle\b", r"\bmock\b"]),
    ("treated", [r"\btreated\b", r"\btrt\b", r"\bdrug\b"]),
    ("wt", [r"\bwt\b", r"\bwild[\s-]?type\b"]),
    ("ko", [r"\bko\b", r"\bknock[\s-]?out\b"]),
    ("kd", [r"\bkd\b", r"\bknock[\s-]?down\b"]),
    ("oe", [r"\boe\b", r"\boverexpression\b"]),
]

REPLICATE_PATTERNS = [
    re.compile(r"\brep(?:licate)?[_\-\s]?(\d+)\b", re.IGNORECASE),
]

SOURCE_ORDERS = {
    "ena-first": ["ena", "ncbi"],
    "ncbi-first": ["ncbi", "ena"],
    "ena-only": ["ena"],
    "ncbi-only": ["ncbi"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch ENA/NCBI public sequencing metadata and generate helper tables."
    )
    parser.add_argument(
        "--accession",
        action="append",
        default=[],
        help="Public accession to query. Can be supplied multiple times.",
    )
    parser.add_argument(
        "--accession-file",
        help="Text file containing one accession per line. Empty lines and # comments are ignored.",
    )
    parser.add_argument(
        "--source-strategy",
        choices=sorted(SOURCE_ORDERS),
        default="ena-first",
        help="Metadata source priority (default: ena-first).",
    )
    parser.add_argument(
        "--email",
        default=None,
        help="Optional email for NCBI requests.",
    )
    parser.add_argument(
        "--api-key",
        default=None,
        help="Optional NCBI API key.",
    )
    parser.add_argument(
        "--output-prefix",
        default="public_metadata",
        help="Prefix for generated output files (default: public_metadata).",
    )
    parser.add_argument(
        "--format",
        default="tsv,json",
        help="Comma-separated raw output formats to emit. Supported: tsv,json (default: tsv,json).",
    )
    parser.add_argument(
        "--emit-download-manifest",
        action="store_true",
        help="Write a download manifest TSV.",
    )
    parser.add_argument(
        "--emit-samplesheet-template",
        action="store_true",
        help="Write a candidate samplesheet CSV.",
    )
    parser.add_argument(
        "--strandedness-default",
        default="auto",
        help="Default strandedness in the generated samplesheet (default: auto).",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=30,
        help="HTTP timeout in seconds (default: 30).",
    )
    return parser.parse_args()


def read_accessions(args: argparse.Namespace) -> List[str]:
    accessions = [value.strip() for value in args.accession if value and value.strip()]
    if args.accession_file:
        for line in Path(args.accession_file).read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            accessions.append(line)
    seen = set()
    ordered: List[str] = []
    for accession in accessions:
        upper = accession.strip().upper()
        if upper and upper not in seen:
            seen.add(upper)
            ordered.append(upper)
    if not ordered:
        raise SystemExit("No accession was provided. Use --accession and/or --accession-file.")
    return ordered


def fetch_json(url: str, timeout: int) -> List[dict]:
    request = urllib.request.Request(url, headers={"User-Agent": "riboseq-public-metadata/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8", "replace"))


def fetch_text(url: str, timeout: int) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "riboseq-public-metadata/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", "replace")


GEO_ACCESSION_RE = re.compile(r"^(GSE|GDS)\d+$", re.IGNORECASE)


def is_geo_accession(accession: str) -> bool:
    """Return True if *accession* looks like a GEO Series or DataSet identifier."""
    return bool(GEO_ACCESSION_RE.match(accession))


def resolve_geo_to_sra(accession: str, timeout: int) -> List[str]:
    """Resolve a GEO accession (GSE/GDS) to SRA study accessions via NCBI E-utilities.

    Uses three-step resolution: esearch → elink → esummary, extracting
    ``studyacc`` (SRP/ERP/DRP) entries so the caller can re-query ENA for
    comprehensive run-level metadata including FTP URLs.
    """
    # Step 1 – esearch: GEO term → GEO UID list
    search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=gds&term={accession}&retmode=json"
    )
    try:
        search_result = fetch_json(search_url, timeout)
    except Exception:
        return []
    geo_ids = search_result.get("esearchresult", {}).get("idlist", [])
    if not geo_ids:
        return []

    # Step 2 – elink: GEO UID → linked SRA UIDs
    link_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
        f"?dbfrom=gds&db=sra&id={','.join(geo_ids)}&retmode=json"
    )
    try:
        link_result = fetch_json(link_url, timeout)
    except Exception:
        return []
    sra_ids: List[str] = []
    for linkset in link_result.get("linksets", []):
        for linksetdb in linkset.get("linksetdbs", []):
            if linksetdb.get("linkname") == "gds_sra":
                sra_ids.extend(linksetdb.get("links", []))

    if not sra_ids:
        return []

    # Step 3 – esummary (db=sra): extract study accessions
    summary_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        f"?db=sra&id={','.join(sra_ids)}&retmode=json"
    )
    try:
        summary_result = fetch_json(summary_url, timeout)
    except Exception:
        return []

    study_accessions: List[str] = []
    seen_studies: set = set()
    for uid, summary in summary_result.get("result", {}).items():
        if uid == "uids":
            continue
        expxml = summary.get("expxml", "")
        if not expxml:
            continue
        # expxml is an XML fragment; extract Study/@acc via regex
        match = re.search(r'<Study\s+acc="([^"]+)"', expxml)
        if match:
            study_acc = match.group(1)
            if study_acc not in seen_studies:
                seen_studies.add(study_acc)
                study_accessions.append(study_acc)

    return study_accessions


def accession_query_field(accession: str) -> Optional[Tuple[str, str]]:
    if re.match(r"^(SRR|ERR|DRR)\d+$", accession):
        return "run_accession", accession
    if re.match(r"^(SRX|ERX|DRX)\d+$", accession):
        return "experiment_accession", accession
    if re.match(r"^(SRP|ERP|DRP)\d+$", accession):
        return "secondary_study_accession", accession
    if re.match(r"^(PRJNA|PRJEB|PRJDB|PRJDA)\d+$", accession):
        return "study_accession", accession
    if re.match(r"^(SAMN|SAMEA|SAMD)\d+$", accession):
        return "sample_accession", accession
    if re.match(r"^(SRS|ERS|DRS)\d+$", accession):
        return "secondary_sample_accession", accession
    return None


def query_ena(accession: str, timeout: int) -> List[dict]:
    field_query = accession_query_field(accession)
    if not field_query:
        return []
    field, value = field_query
    params = {
        "result": "read_run",
        "fields": ",".join(ENA_FIELDS),
        "query": f'{field}="{value}"',
        "format": "json",
    }
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + urllib.parse.urlencode(params)
    return fetch_json(url, timeout)


def query_ncbi_runinfo(accession: str, timeout: int, email: Optional[str], api_key: Optional[str]) -> List[dict]:
    params = {"acc": accession}
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    url = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?" + urllib.parse.urlencode(params)
    text = fetch_text(url, timeout)
    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) < 2:
        return []
    return list(csv.DictReader(lines))


def non_empty(value: Optional[str]) -> str:
    return value.strip() if isinstance(value, str) else ""


def split_semicolon(value: str) -> List[str]:
    if not value:
        return []
    return [item.strip() for item in value.split(";") if item.strip()]


def normalize_ena_row(row: dict, accession: str) -> dict:
    ftp_values = split_semicolon(non_empty(row.get("fastq_ftp")))
    md5_values = split_semicolon(non_empty(row.get("fastq_md5")))
    biosample_accession = non_empty(row.get("sample_accession"))
    secondary_sample = non_empty(row.get("secondary_sample_accession"))
    project_accession = non_empty(row.get("study_accession")) or non_empty(row.get("secondary_study_accession"))
    return {
        "input_accession": accession,
        "run_accession": non_empty(row.get("run_accession")),
        "study_accession": non_empty(row.get("secondary_study_accession")) or non_empty(row.get("study_accession")),
        "project_accession": project_accession,
        "experiment_accession": non_empty(row.get("experiment_accession")),
        "sample_accession": secondary_sample or biosample_accession,
        "biosample_accession": biosample_accession,
        "sample_title": non_empty(row.get("sample_title")),
        "experiment_title": non_empty(row.get("experiment_title")),
        "study_title": non_empty(row.get("study_title")),
        "scientific_name": non_empty(row.get("scientific_name")),
        "library_strategy": non_empty(row.get("library_strategy")),
        "library_source": non_empty(row.get("library_source")),
        "library_selection": non_empty(row.get("library_selection")),
        "library_layout": non_empty(row.get("library_layout")),
        "platform": non_empty(row.get("instrument_platform")),
        "model": non_empty(row.get("instrument_model")),
        "fastq_ftp_1": ftp_values[0] if len(ftp_values) >= 1 else "",
        "fastq_ftp_2": ftp_values[1] if len(ftp_values) >= 2 else "",
        "fastq_md5_1": md5_values[0] if len(md5_values) >= 1 else "",
        "fastq_md5_2": md5_values[1] if len(md5_values) >= 2 else "",
        "sra_download_accession": non_empty(row.get("run_accession")),
        "_source": "ena",
        "_source_row": row,
    }


def normalize_ncbi_row(row: dict, accession: str) -> dict:
    download_path = non_empty(row.get("download_path"))
    biosample = non_empty(row.get("BioSample"))
    sample_accession = non_empty(row.get("Sample"))
    study_accession = non_empty(row.get("SRAStudy"))
    project_accession = non_empty(row.get("BioProject")) or study_accession
    return {
        "input_accession": accession,
        "run_accession": non_empty(row.get("Run")),
        "study_accession": study_accession,
        "project_accession": project_accession,
        "experiment_accession": non_empty(row.get("Experiment")),
        "sample_accession": sample_accession or biosample,
        "biosample_accession": biosample or sample_accession,
        "sample_title": non_empty(row.get("SampleName")),
        "experiment_title": "",
        "study_title": "",
        "scientific_name": non_empty(row.get("ScientificName")),
        "library_strategy": non_empty(row.get("LibraryStrategy")),
        "library_source": non_empty(row.get("LibrarySource")),
        "library_selection": non_empty(row.get("LibrarySelection")),
        "library_layout": non_empty(row.get("LibraryLayout")),
        "platform": non_empty(row.get("Platform")),
        "model": non_empty(row.get("Model")),
        "fastq_ftp_1": "",
        "fastq_ftp_2": "",
        "fastq_md5_1": "",
        "fastq_md5_2": "",
        "sra_download_accession": non_empty(row.get("Run")),
        "_download_path": download_path,
        "_source": "ncbi",
        "_source_row": row,
    }


def merge_rows(source_rows: Sequence[dict], source_order: Sequence[str]) -> dict:
    ordered_rows = []
    for source in source_order:
        ordered_rows.extend(row for row in source_rows if row.get("_source") == source)
    if not ordered_rows:
        return {}
    merged = {}
    for row in ordered_rows:
        for key, value in row.items():
            if key.startswith("_"):
                continue
            if value and not merged.get(key):
                merged[key] = value
    merged["_source_rows"] = ordered_rows
    return merged


def infer_type(row: dict) -> Tuple[str, str, str, bool]:
    strategy = non_empty(row.get("library_strategy")).lower()
    texts = " ".join(
        non_empty(row.get(field))
        for field in ("sample_title", "experiment_title", "study_title", "scientific_name")
    ).lower()
    evidence: List[str] = []

    if strategy == "rna-seq":
        evidence.append("library_strategy=RNA-Seq")
        matches = detect_keyword_types(texts)
        if "riboseq" in matches or "tiseq" in matches:
            evidence.append("conflict=title_keywords")
            return "unknown", "low", "; ".join(evidence), True
        return "rnaseq", "high", "; ".join(evidence), False

    if strategy in {"ribo-seq", "riboseq"}:
        evidence.append(f"library_strategy={row.get('library_strategy')}")
        return "riboseq", "high", "; ".join(evidence), False

    matches = detect_keyword_types(texts)
    if len(matches) == 1:
        inferred = matches.pop()
        evidence.append(f"keyword={inferred}")
        if strategy:
            evidence.append(f"library_strategy={row.get('library_strategy')}")
        return inferred, "medium", "; ".join(evidence), False
    if len(matches) > 1:
        evidence.append("conflict=multiple_keyword_classes")
        if strategy:
            evidence.append(f"library_strategy={row.get('library_strategy')}")
        return "unknown", "low", "; ".join(evidence), True

    if strategy and strategy not in {"other", "other_transcriptomic"}:
        evidence.append(f"library_strategy={row.get('library_strategy')}")
    return "unknown", "low", "; ".join(evidence) or "insufficient_evidence", True


def detect_keyword_types(text: str) -> set:
    matches = set()
    for label, patterns in TYPE_PATTERNS.items():
        for pattern in patterns:
            if re.search(pattern, text, flags=re.IGNORECASE):
                matches.add(label)
                break
    return matches


def infer_group(row: dict) -> Tuple[str, str, str]:
    texts = {
        "sample_title": non_empty(row.get("sample_title")),
        "experiment_title": non_empty(row.get("experiment_title")),
        "study_title": non_empty(row.get("study_title")),
    }
    local_text = " ".join([texts["sample_title"], texts["experiment_title"]]).lower()
    group_candidates: List[str] = []
    evidence: List[str] = []

    for label, patterns in GROUP_HINT_PATTERNS:
        for pattern in patterns:
            if re.search(pattern, local_text, flags=re.IGNORECASE):
                group_candidates.append(label)
                evidence.append(f"keyword={label}")
                break

    replicate_hint = extract_replicate_hint(texts["sample_title"]) or extract_replicate_hint(texts["experiment_title"])
    if replicate_hint:
        evidence.append(f"replicate={replicate_hint}")

    if not group_candidates:
        base_group = derive_group_from_sample_title(texts["sample_title"])
        if base_group:
            group_candidates.append(base_group)
            evidence.append(f"sample_title_base={base_group}")

    if group_candidates:
        unique = []
        seen = set()
        for candidate in group_candidates:
            if candidate not in seen:
                seen.add(candidate)
                unique.append(candidate)
        return unique[0], "; ".join(evidence), replicate_hint

    return "", "; ".join(evidence) if evidence else "no_group_hint_detected", replicate_hint


def extract_replicate_hint(text: str) -> str:
    cleaned = non_empty(text)
    for pattern in REPLICATE_PATTERNS:
        match = pattern.search(cleaned)
        if match:
            return f"rep{match.group(1)}"
    return ""


def derive_group_from_sample_title(text: str) -> str:
    cleaned = non_empty(text)
    if not cleaned:
        return ""
    simplified = re.sub(r"\s+", "_", cleaned)
    simplified = re.sub(r"[^A-Za-z0-9_.-]+", "_", simplified).strip("_.-")
    if not simplified:
        return ""
    if re.search(r"[_\-.]\d$", simplified):
        group = re.sub(r"[_\-.]\d$", "", simplified)
    elif re.search(r"[A-Za-z]0\d$", simplified):
        group = re.sub(r"0\d$", "", simplified)
    elif re.search(r"\d$", simplified):
        group = simplified[:-1]
    else:
        group = simplified
    group = group.rstrip("_.-")
    if not group or group == simplified:
        return ""
    return group


def sanitize_sample_name(row: dict) -> str:
    base = non_empty(row.get("sample_title")) or non_empty(row.get("run_accession"))
    sanitized = re.sub(r"\s+", "_", base)
    sanitized = re.sub(r"[^A-Za-z0-9_.-]+", "_", sanitized).strip("_.-")
    return sanitized or non_empty(row.get("run_accession")) or "sample"


def curated_row(merged: dict) -> dict:
    inferred_type, confidence, type_evidence, manual_review = infer_type(merged)
    suggested_group, group_evidence, replicate_hint = infer_group(merged)
    row = {column: "" for column in CURATED_COLUMNS}
    row.update(
        {
            "input_accession": merged.get("input_accession", ""),
            "run_accession": merged.get("run_accession", ""),
            "study_accession": merged.get("study_accession", ""),
            "project_accession": merged.get("project_accession", ""),
            "experiment_accession": merged.get("experiment_accession", ""),
            "sample_accession": merged.get("sample_accession", ""),
            "biosample_accession": merged.get("biosample_accession", ""),
            "sample_title": merged.get("sample_title", ""),
            "experiment_title": merged.get("experiment_title", ""),
            "study_title": merged.get("study_title", ""),
            "scientific_name": merged.get("scientific_name", ""),
            "library_strategy": merged.get("library_strategy", ""),
            "library_source": merged.get("library_source", ""),
            "library_selection": merged.get("library_selection", ""),
            "library_layout": merged.get("library_layout", ""),
            "platform": merged.get("platform", ""),
            "model": merged.get("model", ""),
            "fastq_ftp_1": merged.get("fastq_ftp_1", ""),
            "fastq_ftp_2": merged.get("fastq_ftp_2", ""),
            "fastq_md5_1": merged.get("fastq_md5_1", ""),
            "fastq_md5_2": merged.get("fastq_md5_2", ""),
            "sra_download_accession": merged.get("sra_download_accession", ""),
            "inferred_type": inferred_type,
            "type_confidence": confidence,
            "type_evidence": type_evidence,
            "suggested_group": suggested_group,
            "group_evidence": group_evidence,
            "replicate_hint": replicate_hint,
            "needs_manual_review": "true" if manual_review else "false",
        }
    )
    return row


def write_tsv(path: Path, rows: Sequence[dict], fieldnames: Sequence[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def write_json(path: Path, rows: Sequence[dict]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(rows, handle, indent=2, ensure_ascii=False)


def build_download_manifest(rows: Sequence[dict]) -> List[dict]:
    manifest: List[dict] = []
    for row in rows:
        ftp_pairs = [
            ("R1", row.get("fastq_ftp_1", ""), row.get("fastq_md5_1", "")),
            ("R2", row.get("fastq_ftp_2", ""), row.get("fastq_md5_2", "")),
        ]
        added = False
        for role, ftp_url, md5 in ftp_pairs:
            if ftp_url:
                manifest.append(
                    {
                        "input_accession": row.get("input_accession", ""),
                        "run_accession": row.get("run_accession", ""),
                        "file_role": role,
                        "file_name": Path(ftp_url).name,
                        "ftp_url": ftp_url,
                        "md5": md5,
                        "recommended_method": "ena-ftp",
                        "fallback_sra_accession": row.get("sra_download_accession", ""),
                    }
                )
                added = True
        if not added:
            manifest.append(
                {
                    "input_accession": row.get("input_accession", ""),
                    "run_accession": row.get("run_accession", ""),
                    "file_role": "SRA",
                    "file_name": f"{row.get('sra_download_accession', '')}.sra",
                    "ftp_url": "",
                    "md5": "",
                    "recommended_method": "ncbi-sra",
                    "fallback_sra_accession": row.get("sra_download_accession", ""),
                }
            )
    return manifest


def build_samplesheet(rows: Sequence[dict], strandedness_default: str) -> List[dict]:
    counts: Dict[str, int] = defaultdict(int)
    sheet: List[dict] = []
    for row in rows:
        base_sample = sanitize_sample_name(row)
        counts[base_sample] += 1
        sample_name = base_sample if counts[base_sample] == 1 else f"{base_sample}_{counts[base_sample]}"
        sheet.append(
            {
                "sample": sample_name,
                "fastq_1": row.get("fastq_ftp_1", "") or row.get("sra_download_accession", ""),
                "fastq_2": row.get("fastq_ftp_2", ""),
                "strandedness": strandedness_default,
                "type": row.get("inferred_type", ""),
                "group": row.get("suggested_group", ""),
                "run_accession": row.get("run_accession", ""),
                "input_accession": row.get("input_accession", ""),
            }
        )
    return sheet


def main() -> int:
    args = parse_args()
    accessions = read_accessions(args)
    source_order = SOURCE_ORDERS[args.source_strategy]
    requested_formats = {item.strip() for item in args.format.split(",") if item.strip()}
    if not requested_formats.issubset({"tsv", "json"}):
        raise SystemExit("Unsupported --format value. Supported values: tsv,json")

    # Expand GEO accessions (GSE/GDS) to SRA study accessions
    expanded_accessions: List[str] = []
    warnings: List[str] = []
    for accession in accessions:
        if is_geo_accession(accession):
            sra_studies = resolve_geo_to_sra(accession, args.timeout)
            if sra_studies:
                for study in sra_studies:
                    if study not in expanded_accessions:
                        expanded_accessions.append(study)
            else:
                warnings.append(f"{accession}\tgeo\tCould not resolve GEO accession to SRA study")
        else:
            if accession not in expanded_accessions:
                expanded_accessions.append(accession)

    if not expanded_accessions:
        if warnings:
            for w in warnings:
                print(w, file=sys.stderr)
        raise SystemExit("No accessions to query after GEO expansion.")

    raw_rows: List[dict] = []
    curated_rows: List[dict] = []

    for accession in expanded_accessions:
        source_rows: List[dict] = []
        for source in source_order:
            try:
                if source == "ena":
                    for row in query_ena(accession, args.timeout):
                        source_rows.append(normalize_ena_row(row, accession))
                elif source == "ncbi":
                    for row in query_ncbi_runinfo(accession, args.timeout, args.email, args.api_key):
                        source_rows.append(normalize_ncbi_row(row, accession))
            except urllib.error.HTTPError as exc:
                warnings.append(f"{accession}\t{source}\tHTTPError {exc.code}: {exc.reason}")
            except urllib.error.URLError as exc:
                warnings.append(f"{accession}\t{source}\tURLError: {exc.reason}")
            except Exception as exc:  # pragma: no cover - defensive I/O guard
                warnings.append(f"{accession}\t{source}\t{type(exc).__name__}: {exc}")

        if not source_rows:
            warnings.append(f"{accession}\tall\tNo records found")
            continue

        grouped: Dict[str, List[dict]] = defaultdict(list)
        for row in source_rows:
            run_accession = row.get("run_accession", "")
            if run_accession:
                grouped[run_accession].append(row)
                raw_row = dict(row)
                raw_row["_source_row"] = json.dumps(row.get("_source_row", {}), ensure_ascii=False)
                raw_rows.append(raw_row)
            else:
                warnings.append(f"{accession}\t{row.get('_source', 'unknown')}\tMissing run accession in source row")

        for run_accession, rows_for_run in grouped.items():
            merged = merge_rows(rows_for_run, source_order)
            if not merged:
                warnings.append(f"{accession}\tmerge\tFailed to merge {run_accession}")
                continue
            curated_rows.append(curated_row(merged))

    curated_rows.sort(key=lambda row: (row["input_accession"], row["run_accession"]))
    raw_rows.sort(key=lambda row: (row["input_accession"], row["run_accession"], row.get("_source", "")))

    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    raw_base = prefix.with_suffix("")

    if "tsv" in requested_formats:
        raw_fieldnames = sorted({key for row in raw_rows for key in row.keys()}) if raw_rows else [
            "input_accession",
            "run_accession",
            "_source",
            "_source_row",
        ]
        write_tsv(raw_base.parent / f"{raw_base.name}.metadata_raw.tsv", raw_rows, raw_fieldnames)
    if "json" in requested_formats:
        json_ready = []
        for row in raw_rows:
            converted = dict(row)
            source_json = converted.get("_source_row", "")
            if source_json:
                converted["_source_row"] = json.loads(source_json)
            json_ready.append(converted)
        write_json(raw_base.parent / f"{raw_base.name}.metadata_raw.json", json_ready)

    write_tsv(raw_base.parent / f"{raw_base.name}.metadata_curated.tsv", curated_rows, CURATED_COLUMNS)

    if args.emit_download_manifest:
        manifest = build_download_manifest(curated_rows)
        write_tsv(raw_base.parent / f"{raw_base.name}.downloads.tsv", manifest, DOWNLOAD_COLUMNS)

    if args.emit_samplesheet_template:
        samplesheet_rows = build_samplesheet(curated_rows, args.strandedness_default)
        with (raw_base.parent / f"{raw_base.name}.samplesheet.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=SAMPLESHEET_COLUMNS)
            writer.writeheader()
            writer.writerows(samplesheet_rows)

    warnings_path = raw_base.parent / f"{raw_base.name}.warnings.tsv"
    with warnings_path.open("w", encoding="utf-8") as handle:
        handle.write("accession\tsource\tmessage\n")
        for line in warnings:
            handle.write(f"{line}\n")

    print(f"Wrote curated metadata: {raw_base.parent / f'{raw_base.name}.metadata_curated.tsv'}")
    if "tsv" in requested_formats:
        print(f"Wrote raw metadata TSV: {raw_base.parent / f'{raw_base.name}.metadata_raw.tsv'}")
    if "json" in requested_formats:
        print(f"Wrote raw metadata JSON: {raw_base.parent / f'{raw_base.name}.metadata_raw.json'}")
    if args.emit_download_manifest:
        print(f"Wrote download manifest: {raw_base.parent / f'{raw_base.name}.downloads.tsv'}")
    if args.emit_samplesheet_template:
        print(f"Wrote samplesheet template: {raw_base.parent / f'{raw_base.name}.samplesheet.csv'}")
    if warnings:
        print(f"Warnings recorded in: {warnings_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
