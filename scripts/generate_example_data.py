"""Generate a small publishable dataset with known phenotype-linked variants."""

from __future__ import annotations

import csv
from pathlib import Path


REFERENCE_NAME = "synthetic_amplicon"
REFERENCE = (
    "ACGTTGCAAGTCGATCGTACGATGCTAGCTAGGCTAACGTTACCGATGATCGTACCGTTAAGC"
    "GATTCGACCTAGGTACCTGATCGTAGCTAGCATGCTACGATTCGGTACAGTCCGATGCTAGCA"
    "TACGGTAGCTAACGTGACCT"
)

SAMPLES = (
    ("black_01", "Black", "AACGTA", 14, "C"),
    ("black_02", "Black", "CCGTAC", 14, "C"),
    ("orange_01", "Orange", "GATTCA", 84, "C"),
    ("orange_02", "Orange", "TGCAGT", 84, "C"),
    ("yellow_01", "Yellow", "AGCTTG", 121, "G"),
    ("yellow_02", "Yellow", "CTAGAC", 121, "G"),
    ("green_01", "Green", "GTACCA", 135, "A"),
    ("green_02", "Green", "TATGCG", 135, "A"),
)


def mutate(sequence: str, position: int, alternate: str) -> str:
    index = position - 1
    if sequence[index] == alternate:
        raise ValueError("Alternate allele must differ from the reference")
    return sequence[:index] + alternate + sequence[index + 1 :]


def generate(output_dir: Path, reads_per_sample: int = 40) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "reference.fa").write_text(
        f">{REFERENCE_NAME}\n{REFERENCE}\n", encoding="utf-8"
    )

    with (output_dir / "samples.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "phenotype", "barcode"])
        for sample, phenotype, barcode, _, _ in SAMPLES:
            writer.writerow([sample, phenotype, barcode])

    with (output_dir / "pooled_reads.fastq").open("w", encoding="utf-8") as handle:
        for sample, _, barcode, position, alternate in SAMPLES:
            variant_sequence = mutate(REFERENCE, position, alternate)
            for read_number in range(reads_per_sample):
                biological_sequence = (
                    variant_sequence if read_number < reads_per_sample // 2 else REFERENCE
                )
                # Three deliberately poor-quality terminal bases verify Q20 trimming.
                sequence = barcode + biological_sequence + "TGA"
                quality = "I" * (len(barcode) + len(biological_sequence)) + "+##"
                handle.write(
                    f"@{sample}_read_{read_number + 1}\n{sequence}\n+\n{quality}\n"
                )

    with (output_dir / "expected_variants.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "phenotype", "position_1_based", "reference", "alternate"])
        for sample, phenotype, _, position, alternate in SAMPLES:
            writer.writerow([sample, phenotype, position, REFERENCE[position - 1], alternate])


if __name__ == "__main__":
    generate(Path(__file__).resolve().parents[1] / "example_data")
