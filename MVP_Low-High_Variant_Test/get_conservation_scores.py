#!/usr/bin/env python3

"""
Fetch conservation scores (phastCons, phyloP, GERP) for genomic variants.

Primary source : myVariant.info REST API (CADD fields, hg19 coordinates)
Fallback source: UCSC REST API bigwig tracks (position-based, hg38)

Scores retrieved:
  phastCons100way  - conservation probability [0, 1], vertebrate 100-way alignment
  phyloP100way     - conservation rate (positive=conserved, negative=accelerated)
  GERP_RS          - Genomic Evolutionary Rate Profiling RS score [~-12, +6]

Usage:
    python get_conservation_scores.py -i MVP_matched_variants.csv -o conservation_scores.csv
    python get_conservation_scores.py --help
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
    "cadd.phast_cons",
    "cadd.phylop",
    "cadd.gerp",
]

UCSC_TRACKS = {
    "phastCons": "phastCons100way",
    "phyloP":    "phyloP100way",
    "GERP":      "allHg38RS_BW",
}

DEFAULT_INPUT  = "C:/Users/mitch/Documents/Argonne/Variant_Scoring/MVP_matched_variants.csv"
DEFAULT_OUTPUT = "C:/Users/mitch/Documents/Argonne/Variant_Scoring/MVP_conservation_scores.csv"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Fetch conservation scores (phastCons, phyloP, GERP) for genomic variants.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i MVP_matched_variants.csv -o MVP_conservation_scores.csv
  %(prog)s --input variants.csv --batch-size 500 --delay 2.0
  %(prog)s --no-ucsc-fallback  # skip UCSC fallback, faster but lower coverage
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
                        help='Threads for UCSC fallback requests (default: 4)')
    parser.add_argument('--ucsc-delay',     type=float, default=0.25,
                        help='Seconds between UCSC requests per thread (default: 0.25)')
    parser.add_argument('--no-ucsc-fallback', action='store_true',
                        help='Disable UCSC REST API fallback for missing variants')
    return parser.parse_args()

# ---------------------------------------------------------------------------
# HGVS ID construction
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

def _coerce_scalar(val, use_max_abs=False):
    """Coerce a dbNSFP value (scalar or list) to a single float."""
    if val is None:
        return None
    if isinstance(val, list):
        valid = [v for v in val if v is not None]
        if not valid:
            return None
        return max(valid, key=abs) if use_max_abs else max(valid)
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def parse_myvariant_hit(hit):
    """Extract phastCons, phyloP, GERP from a single myVariant.info hit dict (CADD fields)."""
    result = {"phastCons": None, "phyloP": None, "GERP": None}
    if not hit or hit.get("notfound"):
        return result
    cadd = hit.get("cadd", {})
    if not cadd:
        return result
    # phast_cons and phylop are dicts with keys: vertebrate, mammalian, primate
    phast = cadd.get("phast_cons", {})
    phylop = cadd.get("phylop", {})
    gerp = cadd.get("gerp", {})
    result["phastCons"] = _coerce_scalar(phast.get("vertebrate") if isinstance(phast, dict) else phast)
    result["phyloP"]    = _coerce_scalar(phylop.get("vertebrate") if isinstance(phylop, dict) else phylop,
                                         use_max_abs=True)
    # gerp.s is the GERP++ RS score (rejected substitutions); gerp.n is neutral rate
    result["GERP"]      = _coerce_scalar(gerp.get("s") if isinstance(gerp, dict) else gerp)
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
    dict : hgvs_id -> {"phastCons": float|None, "phyloP": float|None, "GERP": float|None}
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
        print(f"{found}/{len(batch)} with scores")
        if batch_idx < len(batches):
            time.sleep(delay)

    return scores

# ---------------------------------------------------------------------------
# UCSC REST API — positional fallback
# ---------------------------------------------------------------------------

def _fetch_ucsc_track(chrom, pos_hg38, ucsc_track_id, delay=0.25):
    """Fetch a single track value from the UCSC REST API. Returns float or None."""
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


def _fetch_ucsc_all_tracks(row, ucsc_delay=0.25):
    """Fetch phastCons, phyloP, and GERP from UCSC for one variant row."""
    chrom  = row["CHR"]
    pos38  = row["BP38"]
    mvp_id = row["MVP ID"]

    if pd.isna(pos38):
        return mvp_id, {"phastCons": None, "phyloP": None, "GERP": None}

    result = {}
    track_map = {
        "phastCons": UCSC_TRACKS["phastCons"],
        "phyloP":    UCSC_TRACKS["phyloP"],
        "GERP":      UCSC_TRACKS["GERP"],
    }
    for score_key, ucsc_id in track_map.items():
        result[score_key] = _fetch_ucsc_track(chrom, pos38, ucsc_id, ucsc_delay)

    return mvp_id, result


def fetch_ucsc_fallback(df_missing, max_workers=4, ucsc_delay=0.25):
    """
    Run UCSC REST fallback for variants missing from myVariant.info.
    Uses ThreadPoolExecutor for concurrent requests.

    Returns
    -------
    dict : mvp_id -> {"phastCons": float|None, "phyloP": float|None, "GERP": float|None}
    """
    if df_missing.empty:
        return {}

    print(f"\nUCSC fallback: {len(df_missing):,} variants "
          f"with {max_workers} threads...")

    all_scores = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(_fetch_ucsc_all_tracks, row, ucsc_delay): row["MVP ID"]
            for _, row in df_missing.iterrows()
        }
        done = 0
        for future in as_completed(futures):
            mvp_id, scores = future.result()
            all_scores[mvp_id] = scores
            done += 1
            if done % 100 == 0 or done == len(df_missing):
                print(f"  UCSC: {done}/{len(df_missing)} complete", flush=True)

    return all_scores

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

    # Identify variants with no scores from myVariant.info
    def has_any_score(hgvs):
        return any(v is not None for v in mv_scores.get(hgvs, {}).values())

    df["_has_mv"] = df["_hgvs_id"].apply(has_any_score)
    n_found = df["_has_mv"].sum()
    n_missing = len(df) - n_found
    print(f"\nmyVariant.info: {n_found:,}/{len(df):,} variants with at least one score")
    if n_missing:
        print(f"  {n_missing:,} variants will use UCSC fallback")

    # Fallback: UCSC REST API
    ucsc_scores = {}
    if not args.no_ucsc_fallback and n_missing > 0:
        ucsc_scores = fetch_ucsc_fallback(
            df[~df["_has_mv"]].copy(),
            max_workers=args.max_workers,
            ucsc_delay=args.ucsc_delay,
        )

    # Merge scores onto DataFrame
    def get_score(row, key):
        val = mv_scores.get(row["_hgvs_id"], {}).get(key)
        if val is None:
            val = ucsc_scores.get(row["MVP ID"], {}).get(key)
        return val

    print("\nMerging scores...")
    df["phastCons100way"] = df.apply(lambda r: get_score(r, "phastCons"), axis=1)
    df["phyloP100way"]    = df.apply(lambda r: get_score(r, "phyloP"),    axis=1)
    df["GERP_RS"]         = df.apply(lambda r: get_score(r, "GERP"),      axis=1)

    # Drop working columns
    df.drop(columns=["REF", "ALT", "_hgvs_id", "_has_mv"], errors="ignore", inplace=True)

    # Save
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"Output saved to: {args.output}")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for col in ["phastCons100way", "phyloP100way", "GERP_RS"]:
        n = df[col].notna().sum()
        pct = 100 * n / len(df)
        print(f"  {col}: {n:,}/{len(df):,} ({pct:.1f}%) variants with scores")
        if n > 0:
            print(f"    mean={df[col].mean():.4f}  "
                  f"median={df[col].median():.4f}  "
                  f"range=[{df[col].min():.4f}, {df[col].max():.4f}]")

    # Per-set breakdown if Set column present
    if "Set" in df.columns:
        print()
        for set_name, group in df.groupby("Set"):
            print(f"  {set_name} ({len(group):,} variants):")
            for col in ["phastCons100way", "phyloP100way", "GERP_RS"]:
                valid = group[col].dropna()
                if len(valid):
                    print(f"    {col}: mean={valid.mean():.4f}  median={valid.median():.4f}")

    print("\nDone!")


if __name__ == "__main__":
    main()
