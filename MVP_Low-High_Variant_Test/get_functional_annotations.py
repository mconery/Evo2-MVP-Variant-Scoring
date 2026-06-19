#!/usr/bin/env python3

"""
Fetch functional annotation scores for genomic variants.

Primary source  : myVariant.info REST API (CADD, dbNSFP, gnomAD)
Secondary source: UCSC REST API (CADD phred bigWig fallback; cCRE bigBed for all variants)
Tertiary source : Local RegulomeDB TSV file (positional regulatory evidence)

Scores retrieved:
  CADD_phred         - CADD Phred-scaled score [0, 99]; higher = more deleterious
  CADD_raw           - CADD raw score (pre-scaling)
  DANN_score         - DANN deep learning score [0, 1]; higher = more deleterious
  gnomAD_AF          - gnomAD genome allele frequency [0, 1]
  FitCons_score      - Fitness consequence score [0, 1]; all variants
  GenoCanyon_score   - GenoCanyon non-coding functional score [0, 1]; all variants
  FATHMM_XF          - FATHMM-XF coding pathogenicity score [0, 1]; all variants (dbNSFP)
  ENCODE_cCRE        - ENCODE cCRE functional category (PLS/pELS/dELS/CTCF-only/etc.) or None
  RegulomeDB_rank    - RegulomeDB rank string (1a = strongest, 7 = none); from local file
  RegulomeDB_prob    - RegulomeDB probability score [0, 1]; from local file
  SIFT_score         - SIFT [0, 1]; lower = more deleterious (missense_variant only)
  PolyPhen2_HVAR     - PolyPhen2 HVAR [0, 1]; higher = more damaging (missense_variant only)
  REVEL_score        - REVEL ensemble [0, 1]; higher = more pathogenic (missense_variant only)
  AlphaMissense      - AlphaMissense pathogenicity [0, 1]; higher = more pathogenic (missense_variant only)
  BayesDel_score     - BayesDel addAF score; higher = more deleterious (missense_variant only)
  PrimateAI_score    - PrimateAI pathogenicity [0, 1]; higher = more pathogenic (missense_variant only)

Usage:
    python get_functional_annotations.py -i MVP_matched_variants.csv -o functional_annotations.csv
    python get_functional_annotations.py --help
"""

import argparse
import os
import sys
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import requests

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MYVARIANT_URL = "https://myvariant.info/v1/variant"
UCSC_API_URL  = "https://api.genome.ucsc.edu/getData/track"
UCSC_GENOME   = "hg38"

MYVARIANT_FIELDS = [
    "cadd.phred",
    "cadd.rawscore",
    "dbnsfp.sift.score",
    "dbnsfp.polyphen2.hvar.score",
    "dbnsfp.revel.score",
    "dbnsfp.dann.score",
    "dbnsfp.alphamissense.score",
    "dbnsfp.fitcons.integrated.score",
    "dbnsfp.genocanyon.score",
    "dbnsfp.fathmm-xf.coding_score",
    "dbnsfp.bayesdel.add_af.score",
    "dbnsfp.primateai.score",
    "gnomad_genome.af.af",
]

UCSC_TRACKS = {
    "cadd_phred": "caddPhred",           # bigWig — CADD phred fallback for missing myVariant hits
    "ccre":       "encodeCcreCombined",  # bigBed — queried for all variants
}

# VEP consequence term that warrants missense-specific scores
MISSENSE_TERM = "missense_variant"

DEFAULT_INPUT      = "C:/Users/mitch/Documents/Argonne/Variant_Scoring/MVP_matched_variants.csv"
DEFAULT_OUTPUT     = "C:/Users/mitch/Documents/Argonne/Variant_Scoring/MVP_functional_annotations.csv"
DEFAULT_REGULOMEDB = "C:/Users/mitch/Documents/Argonne/Variant_Scoring/input_data/regulomedb_scores.tsv"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Fetch functional annotation scores (CADD, SpliceAI, dbNSFP, etc.) for genomic variants.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i MVP_matched_variants.csv -o MVP_functional_annotations.csv
  %(prog)s --input variants.csv --batch-size 500 --delay 2.0
  %(prog)s --no-ucsc-fallback  # skip UCSC queries (no CADD fallback, no cCRE)
        """
    )
    parser.add_argument('-i', '--input',    default=DEFAULT_INPUT,
                        help='Input CSV file path (default: %(default)s)')
    parser.add_argument('-o', '--output',   default=DEFAULT_OUTPUT,
                        help='Output CSV file path (default: %(default)s)')
    parser.add_argument('--batch-size',     type=int,   default=1000,
                        help='Variants per myVariant.info POST request (default: 1000)')
    parser.add_argument('--delay',          type=float, default=1.0,
                        help='Seconds between myVariant.info batches (default: 1.0)')
    parser.add_argument('--max-workers',    type=int,   default=4,
                        help='Threads for UCSC requests (default: 4)')
    parser.add_argument('--ucsc-delay',     type=float, default=0.25,
                        help='Seconds between UCSC requests per thread (default: 0.25)')
    parser.add_argument('--no-ucsc-fallback', action='store_true',
                        help='Disable UCSC REST API queries (no CADD phred fallback, no cCRE, no LINSIGHT)')
    parser.add_argument('--regulomedb',     default=DEFAULT_REGULOMEDB,
                        help='Path to RegulomeDB TSV file; skipped if file not found '
                             '(default: %(default)s)')
    return parser.parse_args()

# ---------------------------------------------------------------------------
# HGVS ID construction  (identical to get_conservation_scores.py)
# ---------------------------------------------------------------------------

def build_hgvs_id(chrom, pos_hg19, ref, alt):
    """
    Build a myVariant.info HGVS-like identifier (hg19 coordinates).

    SNP:       chr4:g.174899324C>T
    Deletion:  chr4:g.174899324_174899326del
    Insertion: chr4:g.174899324_174899325insAT
    """
    chrom_str = f"chr{chrom}"
    pos = int(pos_hg19)
    if len(ref) == 1 and len(alt) == 1:
        return f"{chrom_str}:g.{pos}{ref}>{alt}"
    elif len(ref) > len(alt):
        end = pos + len(ref) - 1
        return f"{chrom_str}:g.{pos}_{end}del"
    else:
        inserted = alt[len(ref):]
        return f"{chrom_str}:g.{pos}_{pos + 1}ins{inserted}"

# ---------------------------------------------------------------------------
# myVariant.info — primary fetcher
# ---------------------------------------------------------------------------

def _coerce_scalar(val, use_min=False, use_max_abs=False):
    """Coerce a dbNSFP value (scalar or list) to a single float."""
    if val is None:
        return None
    if isinstance(val, list):
        valid = []
        for v in val:
            try:
                valid.append(float(v))
            except (TypeError, ValueError):
                pass
        if not valid:
            return None
        if use_max_abs:
            return max(valid, key=abs)
        if use_min:
            return min(valid)
        return max(valid)
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def parse_myvariant_hit(hit):
    """Extract all functional annotation fields from a single myVariant.info hit dict."""
    result = {
        "cadd_phred":       None,
        "cadd_raw":         None,
        "dann_score":       None,
        "sift_score":       None,
        "polyphen2_hvar":   None,
        "revel_score":      None,
        "alphamissense":    None,
        "fitcons_score":    None,
        "genocanyon_score": None,
        "fathmm_xf":       None,
        "bayesdel_score":   None,
        "primateai_score":  None,
        "gnomad_af":        None,
    }
    if not hit or hit.get("notfound"):
        return result

    # CADD fields — cadd.rawscore (no underscore) is the correct field name
    cadd = hit.get("cadd", {})
    if isinstance(cadd, list):
        cadd = cadd[0] if cadd else {}
    if cadd:
        result["cadd_phred"] = _coerce_scalar(cadd.get("phred"))
        result["cadd_raw"]   = _coerce_scalar(cadd.get("rawscore"))

    # dbNSFP fields — all are nested sub-objects
    dbnsfp = hit.get("dbnsfp", {})
    if isinstance(dbnsfp, list):
        dbnsfp = dbnsfp[0] if dbnsfp else {}
    if dbnsfp:
        # DANN: dbnsfp.dann.score
        dann_block = dbnsfp.get("dann", {})
        if isinstance(dann_block, dict):
            result["dann_score"] = _coerce_scalar(dann_block.get("score"))

        # SIFT: list of per-transcript scores; take min (lower = more deleterious)
        sift_block = dbnsfp.get("sift", {})
        if isinstance(sift_block, dict):
            result["sift_score"] = _coerce_scalar(sift_block.get("score"), use_min=True)

        # PolyPhen2: take max over transcripts (higher = more damaging)
        pp2_block = dbnsfp.get("polyphen2", {})
        if isinstance(pp2_block, dict):
            hvar = pp2_block.get("hvar", {})
            if isinstance(hvar, dict):
                result["polyphen2_hvar"] = _coerce_scalar(hvar.get("score"))

        # REVEL: take max over transcripts
        revel_block = dbnsfp.get("revel", {})
        if isinstance(revel_block, list):
            revel_block = revel_block[0] if revel_block else {}
        if isinstance(revel_block, dict):
            result["revel_score"] = _coerce_scalar(revel_block.get("score"))

        # AlphaMissense: take max over transcripts
        am_block = dbnsfp.get("alphamissense", {})
        if isinstance(am_block, dict):
            result["alphamissense"] = _coerce_scalar(am_block.get("score"))

        # FitCons: integrated (cell-line-agnostic) score
        fitcons = dbnsfp.get("fitcons", {})
        if isinstance(fitcons, dict):
            integrated = fitcons.get("integrated", {})
            if isinstance(integrated, dict):
                result["fitcons_score"] = _coerce_scalar(integrated.get("score"))

        # GenoCanyon: non-coding functional score
        genocanyon = dbnsfp.get("genocanyon", {})
        if isinstance(genocanyon, dict):
            result["genocanyon_score"] = _coerce_scalar(genocanyon.get("score"))

        # FATHMM-XF: coding model score (noncoding model not in myVariant.info dbNSFP); note hyphenated key
        fathmm_xf = dbnsfp.get("fathmm-xf", {})
        if isinstance(fathmm_xf, dict):
            result["fathmm_xf"] = _coerce_scalar(fathmm_xf.get("coding_score"))

        # BayesDel (allele-frequency-aware): take scalar score
        bayesdel = dbnsfp.get("bayesdel", {})
        if isinstance(bayesdel, dict):
            add_af = bayesdel.get("add_af", {})
            if isinstance(add_af, dict):
                result["bayesdel_score"] = _coerce_scalar(add_af.get("score"))

        # PrimateAI: take max over transcripts
        primateai = dbnsfp.get("primateai", {})
        if isinstance(primateai, dict):
            result["primateai_score"] = _coerce_scalar(primateai.get("score"))

    # gnomAD genome AF
    gnomad = hit.get("gnomad_genome", {})
    if isinstance(gnomad, list):
        gnomad = gnomad[0] if gnomad else {}
    if isinstance(gnomad, dict):
        af_block = gnomad.get("af", {})
        if isinstance(af_block, dict):
            result["gnomad_af"] = _coerce_scalar(af_block.get("af"))
        else:
            result["gnomad_af"] = _coerce_scalar(af_block)

    return result


def _post_with_retry(ids, fields, retries=3, backoff=2.0):
    """POST to myVariant.info with exponential-backoff retry on transient errors."""
    payload = {"ids": ids, "fields": ",".join(fields)}
    for attempt in range(1, retries + 1):
        try:
            resp = requests.post(MYVARIANT_URL, json=payload, timeout=60)
            if resp.status_code == 200:
                return resp.json()
            if resp.status_code in (429, 500, 502, 503):
                wait = backoff * attempt
                print(f"    HTTP {resp.status_code}; retrying in {wait:.0f}s "
                      f"(attempt {attempt}/{retries})", file=sys.stderr)
                time.sleep(wait)
            else:
                print(f"    HTTP {resp.status_code} from myVariant.info; "
                      f"skipping batch.", file=sys.stderr)
                return []
        except requests.exceptions.RequestException as exc:
            wait = backoff * attempt
            print(f"    Request error: {exc}; retrying in {wait:.0f}s "
                  f"(attempt {attempt}/{retries})", file=sys.stderr)
            time.sleep(wait)
    print("    Max retries exceeded; scores will be NaN for this batch.", file=sys.stderr)
    return []


def fetch_myvariant_scores(hgvs_ids, batch_size=1000, delay=1.0):
    """
    Query myVariant.info for all variants in batches.

    Returns
    -------
    dict : hgvs_id -> dict of score keys (see parse_myvariant_hit)
    """
    scores = {}
    batches = [hgvs_ids[i:i + batch_size] for i in range(0, len(hgvs_ids), batch_size)]

    print(f"Querying myVariant.info: {len(hgvs_ids):,} variants in {len(batches)} batches...")

    for batch_idx, batch in enumerate(batches, 1):
        print(f"  Batch {batch_idx}/{len(batches)} ({len(batch)} variants)...",
              end=" ", flush=True)
        hits = _post_with_retry(batch, MYVARIANT_FIELDS)
        found = 0
        for hit in hits:
            query_id = hit.get("query", "")
            parsed   = parse_myvariant_hit(hit)
            scores[query_id] = parsed
            if any(v is not None for v in parsed.values()):
                found += 1
        print(f"{found}/{len(batch)} with at least one score")
        if batch_idx < len(batches):
            time.sleep(delay)

    return scores

# ---------------------------------------------------------------------------
# UCSC REST API — CADD phred fallback + cCRE for all variants
# ---------------------------------------------------------------------------

def _fetch_ucsc_track(chrom, pos_hg38, ucsc_track_id, delay=0.25):
    """Fetch a single bigWig track value from the UCSC REST API. Returns float or None."""
    params = {
        "genome": UCSC_GENOME,
        "track":  ucsc_track_id,
        "chrom":  f"chr{chrom}",
        "start":  str(int(pos_hg38) - 1),   # 0-based half-open interval
        "end":    str(int(pos_hg38)),
    }
    time.sleep(delay)
    try:
        resp = requests.get(UCSC_API_URL, params=params, timeout=20)
        if resp.status_code != 200:
            return None
        data = resp.json().get("data", [])
        if data:
            raw = data[0].get("value")
            return float(raw) if raw is not None else None
        return None
    except Exception:
        return None


def _fetch_ucsc_ccre(chrom, pos_hg38, delay=0.25):
    """
    Fetch ENCODE cCRE functional category for a variant from UCSC bigBed track.

    bigBed tracks return results under the track name key (not "data"), and the
    element type is in the "encodeLabel" field (e.g. "PLS", "pELS", "dELS", "CTCF-only").
    Returns category string or None if no overlapping cCRE element.
    """
    params = {
        "genome": UCSC_GENOME,
        "track":  UCSC_TRACKS["ccre"],
        "chrom":  f"chr{chrom}",
        "start":  str(int(pos_hg38) - 1),
        "end":    str(int(pos_hg38)),
    }
    time.sleep(delay)
    try:
        resp = requests.get(UCSC_API_URL, params=params, timeout=20)
        if resp.status_code != 200:
            return None
        # bigBed response uses track name as key, not "data"
        data = resp.json().get(UCSC_TRACKS["ccre"], [])
        if not data:
            return None
        label = data[0].get("encodeLabel") or data[0].get("ccreLabel")
        return str(label).strip() if label else None
    except Exception:
        return None


def _fetch_ucsc_all_tracks(row, need_cadd, ucsc_delay=0.25):
    """
    Fetch CADD phred (if needed) and ENCODE cCRE category from UCSC.
    cCRE is queried for all variants; CADD phred only when need_cadd=True.
    """
    chrom  = row["CHR"]
    pos38  = row["BP38"]
    mvp_id = row["MVP ID"]

    result = {"cadd_phred_ucsc": None, "ccre": None}
    if pd.isna(pos38):
        return mvp_id, result

    if need_cadd:
        result["cadd_phred_ucsc"] = _fetch_ucsc_track(
            chrom, pos38, UCSC_TRACKS["cadd_phred"], ucsc_delay
        )
    result["ccre"] = _fetch_ucsc_ccre(chrom, pos38, ucsc_delay)

    return mvp_id, result


def fetch_ucsc_scores(df, need_cadd_mask, max_workers=4, ucsc_delay=0.25):
    """
    Run UCSC queries for all variants (cCRE) and for variants missing CADD phred.

    Parameters
    ----------
    df            : full variant DataFrame
    need_cadd_mask: boolean Series indexed like df; True where CADD phred is missing
    max_workers   : concurrent request threads
    ucsc_delay    : per-request sleep in seconds

    Returns
    -------
    dict : mvp_id -> {"cadd_phred_ucsc": float|None, "ccre": str|None}
    """
    print(f"\nUCSC queries: {len(df):,} variants for cCRE; "
          f"{need_cadd_mask.sum():,} also need CADD phred fallback "
          f"({max_workers} threads)...")

    all_scores = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _fetch_ucsc_all_tracks,
                row,
                bool(need_cadd_mask.loc[idx]),
                ucsc_delay
            ): row["MVP ID"]
            for idx, row in df.iterrows()
        }
        done = 0
        for future in as_completed(futures):
            mvp_id, scores = future.result()
            all_scores[mvp_id] = scores
            done += 1
            if done % 100 == 0 or done == len(df):
                print(f"  UCSC: {done}/{len(df)} complete", flush=True)

    return all_scores

# ---------------------------------------------------------------------------
# RegulomeDB — local file lookup
# ---------------------------------------------------------------------------

def load_regulomedb(path, target_rsids, target_pos):
    """
    Load RegulomeDB scores for target rsIDs/positions from a local TSV file.

    The file uses hg38 BED format (0-based half-open), chr-prefixed chromosomes.
    Primary match: rsID. Fallback: (chrom, 0-based start) positional match.

    Returns
    -------
    dict : rsid or (chrom, start) -> (ranking: str, probability_score: float)
    """
    print(f"Loading RegulomeDB from: {path}")
    rdb = {}
    chunks = pd.read_csv(
        path, sep='\t',
        usecols=['chrom', 'start', 'rsid', 'ranking', 'probability_score'],
        dtype={'chrom': str, 'start': int, 'rsid': str,
               'ranking': str, 'probability_score': float},
        chunksize=500_000
    )
    for chunk in chunks:
        # rsID match (primary)
        for _, row in chunk[chunk['rsid'].isin(target_rsids)].iterrows():
            if row['rsid'] not in rdb:
                rdb[row['rsid']] = (row['ranking'], row['probability_score'])
        # Positional fallback for variants not matched by rsID
        pos_hits = chunk[chunk.apply(
            lambda r: (r['chrom'], r['start']) in target_pos, axis=1
        )]
        for _, row in pos_hits.iterrows():
            pos_key = (row['chrom'], row['start'])
            if pos_key not in rdb:
                rdb[pos_key] = (row['ranking'], row['probability_score'])
    n_rsid = sum(1 for k in rdb if isinstance(k, str))
    n_pos  = sum(1 for k in rdb if isinstance(k, tuple))
    print(f"  {n_rsid:,} rsID matches + {n_pos:,} positional matches = {len(rdb):,} total")
    return rdb

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_arguments()

    # Validate input
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # Load variants
    print(f"Loading variants from: {args.input}")
    df = pd.read_csv(args.input)
    df["REF"] = df["MVP ID"].apply(lambda x: x.split(":")[2])
    df["ALT"] = df["MVP ID"].apply(lambda x: x.split(":")[3])
    set_counts = df["Set"].value_counts().to_dict() if "Set" in df.columns else {}
    print(f"  {len(df):,} variants loaded {set_counts}")

    # Build HGVS IDs from hg19 coordinates (myVariant.info uses hg19)
    df["_hgvs_id"] = df.apply(
        lambda r: build_hgvs_id(r["CHR"], r["BP"], r["REF"], r["ALT"]), axis=1
    )

    # Primary: myVariant.info
    mv_scores = fetch_myvariant_scores(
        df["_hgvs_id"].tolist(),
        batch_size=args.batch_size,
        delay=args.delay,
    )

    # Identify variants with no scores from myVariant.info (need CADD phred fallback)
    def has_any_score(hgvs):
        return any(v is not None for v in mv_scores.get(hgvs, {}).values())

    df["_has_mv"] = df["_hgvs_id"].apply(has_any_score)
    n_found   = df["_has_mv"].sum()
    n_missing = len(df) - n_found
    print(f"\nmyVariant.info: {n_found:,}/{len(df):,} variants with at least one score")
    if n_missing:
        print(f"  {n_missing:,} variants will use UCSC fallback for CADD phred")

    # UCSC: cCRE + LINSIGHT for all variants; CADD phred fallback for missing
    ucsc_scores = {}
    if not args.no_ucsc_fallback:
        need_cadd_mask = ~df["_has_mv"]
        # Use df index as the mask index for safe .loc access inside worker
        need_cadd_mask.index = df.index
        ucsc_scores = fetch_ucsc_scores(
            df,
            need_cadd_mask,
            max_workers=args.max_workers,
            ucsc_delay=args.ucsc_delay,
        )

    # RegulomeDB from local file
    rdb_scores = {}
    if os.path.exists(args.regulomedb):
        target_rsids = set(df["RSID"].dropna().tolist())
        target_pos   = set(
            (f"chr{r['CHR']}", int(r['BP38']) - 1)
            for _, r in df.dropna(subset=["BP38"]).iterrows()
        )
        rdb_scores = load_regulomedb(args.regulomedb, target_rsids, target_pos)
    else:
        print(f"RegulomeDB file not found, skipping: {args.regulomedb}")

    # Merge scores onto DataFrame
    def get_mv(row, key):
        return mv_scores.get(row["_hgvs_id"], {}).get(key)

    def get_ucsc(row, key):
        return ucsc_scores.get(row["MVP ID"], {}).get(key)

    def get_rdb(row, col_idx):
        hit = rdb_scores.get(row["RSID"])
        if hit is None and pd.notna(row.get("BP38")):
            hit = rdb_scores.get((f"chr{row['CHR']}", int(row["BP38"]) - 1))
        return hit[col_idx] if hit else None

    print("\nMerging scores...")
    # CADD: myVariant primary, UCSC fallback
    df["CADD_phred"] = df.apply(
        lambda r: get_mv(r, "cadd_phred") if get_mv(r, "cadd_phred") is not None
        else get_ucsc(r, "cadd_phred_ucsc"),
        axis=1
    )
    df["CADD_raw"]         = df.apply(lambda r: get_mv(r, "cadd_raw"),         axis=1)
    df["DANN_score"]       = df.apply(lambda r: get_mv(r, "dann_score"),        axis=1)
    df["gnomAD_AF"]        = df.apply(lambda r: get_mv(r, "gnomad_af"),         axis=1)
    # Additional scores from myVariant.info dbNSFP (~6% coverage — dbNSFP near-transcript index)
    df["FitCons_score"]    = df.apply(lambda r: get_mv(r, "fitcons_score"),     axis=1)
    df["GenoCanyon_score"] = df.apply(lambda r: get_mv(r, "genocanyon_score"),  axis=1)
    df["FATHMM_XF"]        = df.apply(lambda r: get_mv(r, "fathmm_xf"),        axis=1)
    # ENCODE cCRE: UCSC bigBed
    df["ENCODE_cCRE"]      = df.apply(lambda r: get_ucsc(r, "ccre"),            axis=1)
    # RegulomeDB: local file
    df["RegulomeDB_rank"]  = df.apply(lambda r: get_rdb(r, 0),                 axis=1)
    df["RegulomeDB_prob"]  = df.apply(lambda r: get_rdb(r, 1),                 axis=1)
    # Missense-specific scores (myVariant primary)
    df["SIFT_score"]       = df.apply(lambda r: get_mv(r, "sift_score"),        axis=1)
    df["PolyPhen2_HVAR"]   = df.apply(lambda r: get_mv(r, "polyphen2_hvar"),    axis=1)
    df["REVEL_score"]      = df.apply(lambda r: get_mv(r, "revel_score"),       axis=1)
    df["AlphaMissense"]    = df.apply(lambda r: get_mv(r, "alphamissense"),     axis=1)
    df["BayesDel_score"]   = df.apply(lambda r: get_mv(r, "bayesdel_score"),    axis=1)
    df["PrimateAI_score"]  = df.apply(lambda r: get_mv(r, "primateai_score"),   axis=1)

    # Gate missense-only scores: NaN for any variant that is not missense_variant
    if "VEP Annotation" in df.columns:
        non_missense = df["VEP Annotation"].str.strip() != MISSENSE_TERM
        for col in ["SIFT_score", "PolyPhen2_HVAR", "REVEL_score",
                    "AlphaMissense", "BayesDel_score", "PrimateAI_score"]:
            df.loc[non_missense, col] = np.nan

    # Drop working columns
    df.drop(columns=["REF", "ALT", "_hgvs_id", "_has_mv"], errors="ignore", inplace=True)

    # Save
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"Output saved to: {args.output}")

    # -----------------------------------------------------------------------
    # Summary statistics
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    numeric_cols = [
        "CADD_phred", "CADD_raw", "DANN_score", "gnomAD_AF",
        "FitCons_score", "GenoCanyon_score", "FATHMM_XF",
        "RegulomeDB_prob",
    ]
    missense_cols = [
        "SIFT_score", "PolyPhen2_HVAR", "REVEL_score",
        "AlphaMissense", "BayesDel_score", "PrimateAI_score",
    ]
    categorical_cols = ["ENCODE_cCRE", "RegulomeDB_rank"]

    # Numeric columns — coverage, mean, median, range
    for col in numeric_cols:
        if col not in df.columns:
            continue
        n = df[col].notna().sum()
        pct = 100 * n / len(df)
        print(f"  {col}: {n:,}/{len(df):,} ({pct:.1f}%) variants with scores")
        if n > 0:
            print(f"    mean={df[col].mean():.4f}  "
                  f"median={df[col].median():.4f}  "
                  f"range=[{df[col].min():.4f}, {df[col].max():.4f}]")

    # Missense-specific columns — denominator is missense variants only
    if "VEP Annotation" in df.columns:
        n_missense = (df["VEP Annotation"].str.strip() == MISSENSE_TERM).sum()
    else:
        n_missense = len(df)
    print(f"\n  Missense-specific scores (denominator: {n_missense:,} missense_variant rows):")
    for col in missense_cols:
        if col not in df.columns:
            continue
        n = df[col].notna().sum()
        pct = 100 * n / n_missense if n_missense > 0 else 0.0
        print(f"  {col}: {n:,}/{n_missense:,} ({pct:.1f}%) missense variants with scores")
        if n > 0:
            print(f"    mean={df[col].mean():.4f}  "
                  f"median={df[col].median():.4f}  "
                  f"range=[{df[col].min():.4f}, {df[col].max():.4f}]")

    # Categorical columns — coverage + top value_counts
    print()
    for col in categorical_cols:
        if col not in df.columns:
            continue
        n = df[col].notna().sum()
        pct = 100 * n / len(df)
        print(f"  {col}: {n:,}/{len(df):,} ({pct:.1f}%) variants with annotation")
        if n > 0:
            vc = df[col].value_counts().head(5)
            for val, cnt in vc.items():
                print(f"    {val}: {cnt:,} ({100 * cnt / len(df):.1f}%)")

    # Per-set breakdown
    if "Set" in df.columns:
        print()
        for set_name, group in df.groupby("Set"):
            print(f"  {set_name} ({len(group):,} variants):")
            for col in numeric_cols + missense_cols:
                if col not in group.columns:
                    continue
                valid = group[col].dropna()
                if len(valid):
                    print(f"    {col}: mean={valid.mean():.4f}  median={valid.median():.4f}")
            for col in categorical_cols:
                if col not in group.columns:
                    continue
                vc = group[col].value_counts().head(3)
                if len(vc):
                    top = ", ".join(f"{v}={c}" for v, c in vc.items())
                    print(f"    {col}: {top}")

    print("\nDone!")


if __name__ == "__main__":
    main()
