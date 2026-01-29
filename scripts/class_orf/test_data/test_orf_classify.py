import csv
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_DIR = os.path.dirname(os.path.abspath(__file__))

sys.path.insert(0, SCRIPT_DIR)

from class_ORFtype import load_gene_level_cds, parse_coords, classify_orfs


def read_orfs_tsv(path):
    records = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            records.append({
                "orf_id": row.get("orf_id"),
                "gene_id": row.get("gene_id"),
                "seqnames": row.get("seqnames"),
                "strand": row.get("strand"),
                "exons": parse_coords(row.get("exons") or ""),
            })
    return records


def read_r_results(path):
    results = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            results[row["orf_id"]] = row.get("ORF_type_py")
    return results


def main():
    cds_path = os.path.join(TEST_DIR, "test_cds.txt")
    orf_path = os.path.join(TEST_DIR, "test_orfs.tsv")
    r_path = os.path.join(TEST_DIR, "r_results.csv")

    cds_data = load_gene_level_cds(cds_path)
    records = read_orfs_tsv(orf_path)
    py_results = classify_orfs(records, cds_data)

    if not os.path.exists(r_path):
        print("R results not found. Run test_orf_classify.R first.")
        sys.exit(1)

    r_results = read_r_results(r_path)

    mismatches = []
    for rec, py_cls in zip(records, py_results):
        orf_id = rec["orf_id"]
        r_cls = r_results.get(orf_id)
        if r_cls != py_cls:
            mismatches.append((orf_id, py_cls, r_cls))

    if mismatches:
        print("Mismatches detected:")
        for orf_id, py_cls, r_cls in mismatches:
            print(f"{orf_id}: python={py_cls}, r={r_cls}")
        sys.exit(2)

    print("All classifications match between Python and R.")


if __name__ == "__main__":
    main()
