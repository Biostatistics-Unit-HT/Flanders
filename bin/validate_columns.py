#!/usr/bin/env python

import argparse
import glob
import json
import sys

import pandas as pd


def is_missing(value):
    return pd.isna(value) or str(value).strip().lower() in {"", "na", "nan", "none"}


def fail(message):
    print(message, file=sys.stderr)
    return 1


def validate_one_row(gwas_path, sdY_file, sdY_is_file, record):
    """Validate one row against its staged GWAS and optional sdY file."""
    errors = []

    if is_missing(record.get("freq_lab")) and is_missing(record.get("n_lab")):
        errors.append(
            f"freq_lab and n_lab are both missing/NA for {gwas_path}; at least one is required."
        )

    if sdY_is_file:
        try:
            sdY_columns = pd.read_csv(sdY_file, sep="\t", nrows=1).columns
            required_sdY_columns = {"sdY", "phenotype_id"}
            if not required_sdY_columns.issubset(sdY_columns):
                errors.append(
                    f"sdY file {sdY_file} must contain 'sdY' and 'phenotype_id' columns."
                )
        except Exception as error:
            errors.append(f"Could not read sdY file {sdY_file}: {error}")
    else:
        sdY_value = record.get("sdY")
        if not is_missing(sdY_value) and pd.isna(
            pd.to_numeric(sdY_value, errors="coerce")
        ):
            errors.append(
                f"sdY value {sdY_value} is invalid; provide NA, a number, or an existing file path."
            )

    try:
        gwas_columns = pd.read_csv(
            gwas_path,
            compression="infer",
            engine="python",
            sep=None,
            nrows=1,
        ).columns
    except Exception as error:
        errors.append(f"Could not read GWAS file {gwas_path}: {error}")
        return errors

    labels = [
        "pos_lab",
        "rsid_lab",
        "chr_lab",
        "a0_lab",
        "a1_lab",
        "effect_lab",
        "se_lab",
        "freq_lab",
        "n_lab",
    ]
    if str(record.get("is_molQTL", "")).upper() == "TRUE":
        labels.append("key")

    columns_to_check = [
        record[label]
        for label in labels
        if label in record and not is_missing(record[label])
    ]
    missing_columns = sorted(set(columns_to_check) - set(gwas_columns))

    print(f"Checking columns {columns_to_check} in file {gwas_path}.")
    if missing_columns:
        errors.append(
            f"Columns {missing_columns} from the input table are not present in {gwas_path}."
        )

    return errors


def main():
    parser = argparse.ArgumentParser(description="Validate all summary-statistics input rows.")
    parser.add_argument("--manifest", required=True, help="JSON-lines manifest for all input rows")
    args = parser.parse_args()

    try:
        with open(args.manifest, encoding="utf-8") as manifest_file:
            records = [json.loads(line) for line in manifest_file if line.strip()]
    except (OSError, ValueError, json.JSONDecodeError) as error:
        return fail(f"Could not read validation manifest {args.manifest}: {error}")

    # Each input file is staged into an indexed directory to avoid basename collisions.
    gwas_files = sorted(glob.glob("gwas*/*"))
    sdY_files = sorted(glob.glob("sdy*/*"))
    if len(records) != len(gwas_files) or len(records) != len(sdY_files):
        return fail(
            "The number of manifest rows, staged GWAS files, and staged sdY files must match "
            f"(rows={len(records)}, GWAS={len(gwas_files)}, sdY={len(sdY_files)})."
        )

    errors = []
    for index, (entry, gwas_path, sdY_file) in enumerate(
        zip(records, gwas_files, sdY_files), start=1
    ):
        record = entry["metadata"]
        sdY_is_file = entry["sdY_is_file"]
        print(f"Validating samplesheet row {index}.")
        errors.extend(validate_one_row(gwas_path, sdY_file, sdY_is_file, record))

    if errors:
        for error in errors:
            print(f"Error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
