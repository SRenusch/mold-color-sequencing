"""Process barcoded amplicon reads and report reference-based SNVs.

The workflow validates and demultiplexes FASTQ reads, optionally trims their
low-quality tails, aligns each sample with BWA, creates indexed BAM files with
samtools, and calls simple SNVs using configurable evidence thresholds.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
import subprocess
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


DNA_PATTERN = re.compile(r"^[ACGTN]+$", re.IGNORECASE)
DEGRADED_TAIL_PATTERN = re.compile(r"[DF]{2,}$")


@dataclass(frozen=True)
class Sample:
    """Metadata required to assign a barcoded read to a sample."""

    name: str
    group: str
    barcode: str


@dataclass(frozen=True)
class FastqRecord:
    """One validated four-line FASTQ record."""

    header: str
    sequence: str
    separator: str
    quality: str


@dataclass(frozen=True)
class VariantCall:
    """A reference-based single-nucleotide variant observed in one sample."""

    contig: str
    position: int
    reference: str
    alternate: str
    depth: int
    reference_count: int
    alternate_count: int
    alternate_fraction: float


def open_text(path: Path) -> TextIO:
    """Open plain-text or gzip-compressed input in text mode."""

    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_samples(
    metadata_path: Path,
    sample_column: str = "Name",
    group_column: str = "Color",
    barcode_column: str = "Barcode",
) -> dict[str, Sample]:
    """Read a tab-delimited sample sheet and index samples by barcode."""

    samples: dict[str, Sample] = {}
    with open_text(metadata_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {sample_column, group_column, barcode_column}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Metadata is missing required column(s): {', '.join(sorted(missing))}"
            )

        for line_number, row in enumerate(reader, start=2):
            name = row[sample_column].strip()
            group = row[group_column].strip()
            barcode = row[barcode_column].strip().upper()
            if not name or not group or not barcode:
                raise ValueError(f"Metadata line {line_number} contains an empty value")
            if not DNA_PATTERN.fullmatch(barcode):
                raise ValueError(f"Invalid barcode {barcode!r} on metadata line {line_number}")
            if barcode in samples:
                raise ValueError(f"Duplicate barcode {barcode!r} in metadata")
            samples[barcode] = Sample(name=name, group=group, barcode=barcode)

    if not samples:
        raise ValueError("Metadata contains no samples")
    return samples


def iter_fastq(path: Path) -> Iterator[FastqRecord]:
    """Yield validated FASTQ records and report malformed input clearly."""

    with open_text(path) as handle:
        record_number = 0
        while True:
            lines = [handle.readline() for _ in range(4)]
            if not any(lines):
                return
            record_number += 1
            if any(line == "" for line in lines):
                raise ValueError(f"Incomplete FASTQ record {record_number}")

            header, sequence, separator, quality = (line.rstrip("\r\n") for line in lines)
            if not header.startswith("@"):
                raise ValueError(f"FASTQ record {record_number} has an invalid header")
            if not separator.startswith("+"):
                raise ValueError(f"FASTQ record {record_number} has an invalid separator")
            if len(sequence) != len(quality):
                raise ValueError(
                    f"FASTQ record {record_number} has unequal sequence and quality lengths"
                )
            if not DNA_PATTERN.fullmatch(sequence):
                raise ValueError(f"FASTQ record {record_number} contains an invalid nucleotide")
            yield FastqRecord(header, sequence.upper(), separator, quality)


def trim_degraded_tail(sequence: str, quality: str) -> tuple[str, str]:
    """Remove a terminal run of two or more assignment-specific D/F symbols."""

    match = DEGRADED_TAIL_PATTERN.search(quality)
    trim_at = match.start() if match else len(quality)
    return sequence[:trim_at], quality[:trim_at]


def trim_phred_tail(
    sequence: str, quality: str, min_quality: int = 20
) -> tuple[str, str]:
    """Trim consecutive terminal bases whose Phred+33 score is below a threshold."""

    if min_quality < 0:
        raise ValueError("min_quality cannot be negative")
    trim_at = len(quality)
    # FASTQ quality characters use Phred+33 encoding: subtracting 33 from
    # the ASCII value gives the estimated base-quality score.
    while trim_at > 0 and ord(quality[trim_at - 1]) - 33 < min_quality:
        trim_at -= 1
    return sequence[:trim_at], quality[:trim_at]


def safe_sample_name(name: str) -> str:
    """Convert a sample label into a safe, predictable filename component."""

    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip()).strip("._")
    if not safe:
        raise ValueError(f"Sample name {name!r} cannot be used as a filename")
    return safe


def demultiplex_fastq(
    reads_path: Path,
    metadata_path: Path,
    output_dir: Path,
    *,
    sample_column: str = "Name",
    group_column: str = "Color",
    barcode_column: str = "Barcode",
    trimming_mode: str = "phred",
    trim_quality: int = 20,
) -> dict[str, int]:
    """Split pooled reads into validated per-sample FASTQ files."""

    samples = read_samples(metadata_path, sample_column, group_column, barcode_column)
    barcode_lengths = {len(barcode) for barcode in samples}
    if len(barcode_lengths) != 1:
        raise ValueError("All barcodes must have the same length")
    barcode_length = barcode_lengths.pop()

    output_dir.mkdir(parents=True, exist_ok=True)
    counts = {sample.name: 0 for sample in samples.values()}
    unmatched = 0

    with ExitStack() as stack:
        outputs = {
            barcode: stack.enter_context(
                (output_dir / f"{safe_sample_name(sample.name)}.fastq").open(
                    "w", encoding="utf-8"
                )
            )
            for barcode, sample in samples.items()
        }

        for record in iter_fastq(reads_path):
            # This workflow expects the sample barcode inline at the 5' start
            # of each read, rather than in a separate Illumina index-read file.
            barcode = record.sequence[:barcode_length]
            sample = samples.get(barcode)
            if sample is None:
                # An unknown barcode can occur without making the FASTQ invalid.
                # Count it for quality control, but never guess a sample assignment.
                unmatched += 1
                continue

            # Remove the barcode from both sequence and quality strings.
            sequence = record.sequence[barcode_length:]
            quality = record.quality[barcode_length:]
            if trimming_mode == "phred":
                sequence, quality = trim_phred_tail(sequence, quality, trim_quality)
            elif trimming_mode == "classroom-df":
                sequence, quality = trim_degraded_tail(sequence, quality)
            elif trimming_mode != "none":
                raise ValueError(
                    "trimming_mode must be 'phred', 'classroom-df', or 'none'"
                )
            if not sequence:
                continue

            outputs[barcode].write(
                f"{record.header}\n{sequence}\n{record.separator}\n{quality}\n"
            )
            counts[sample.name] += 1

    counts["__unmatched__"] = unmatched
    return counts


def write_summary(counts: dict[str, int], output_path: Path) -> None:
    """Write demultiplexing read counts as a tab-delimited table."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "read_count"])
        for sample, count in sorted(counts.items()):
            writer.writerow([sample, count])


def require_executable(name: str) -> str:
    """Return an executable path or raise an actionable error."""

    executable = shutil.which(name)
    if executable is None:
        raise RuntimeError(
            f"Required program {name!r} was not found. Install it and ensure it is on PATH."
        )
    return executable


def ensure_bwa_index(reference_path: Path, bwa: str) -> None:
    """Create the five BWA index files when any are missing."""

    required_suffixes = (".amb", ".ann", ".bwt", ".pac", ".sa")
    missing = [
        Path(f"{reference_path}{suffix}")
        for suffix in required_suffixes
        if not Path(f"{reference_path}{suffix}").is_file()
    ]
    if not missing:
        return

    completed = subprocess.run(
        [bwa, "index", str(reference_path)],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"BWA indexing failed:\n{completed.stderr.strip()}")


def align_fastq(
    fastq_path: Path,
    reference_path: Path,
    output_bam: Path,
    *,
    bwa: str,
    samtools: str,
    threads: int = 1,
) -> None:
    """Align one FASTQ and stream the result into a sorted, indexed BAM."""

    if threads < 1:
        raise ValueError("threads must be at least 1")
    output_bam.parent.mkdir(parents=True, exist_ok=True)

    # Stream BWA's SAM output directly into samtools sort. This avoids writing
    # a large intermediate SAM file and lets either program's failure be checked.
    bwa_process = subprocess.Popen(
        [bwa, "mem", "-t", str(threads), str(reference_path), str(fastq_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if bwa_process.stdout is None or bwa_process.stderr is None:
        raise RuntimeError("Unable to open BWA output streams")

    samtools_process = subprocess.Popen(
        [samtools, "sort", "-@", str(threads), "-o", str(output_bam), "-"],
        stdin=bwa_process.stdout,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    bwa_process.stdout.close()
    _, samtools_stderr = samtools_process.communicate()
    bwa_stderr = bwa_process.stderr.read()
    bwa_returncode = bwa_process.wait()

    if bwa_returncode != 0 or samtools_process.returncode != 0:
        output_bam.unlink(missing_ok=True)
        details = "\n".join(
            message.decode("utf-8", errors="replace").strip()
            for message in (bwa_stderr, samtools_stderr)
            if message
        )
        raise RuntimeError(f"Alignment failed for {fastq_path.name}:\n{details}")

    completed = subprocess.run(
        [samtools, "index", str(output_bam)],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"BAM indexing failed for {output_bam.name}:\n{completed.stderr.strip()}"
        )


def align_sample_fastqs(
    fastq_dir: Path,
    reference_path: Path,
    bam_dir: Path,
    *,
    threads: int = 1,
) -> list[Path]:
    """Align every non-empty sample FASTQ and return the BAM paths."""

    if not reference_path.is_file():
        raise FileNotFoundError(f"Reference FASTA not found: {reference_path}")
    bwa = require_executable("bwa")
    samtools = require_executable("samtools")
    ensure_bwa_index(reference_path, bwa)

    fastq_paths = sorted(fastq_dir.glob("*.fastq"))
    if not fastq_paths:
        raise ValueError(f"No sample FASTQ files found in {fastq_dir}")

    bam_paths: list[Path] = []
    for fastq_path in fastq_paths:
        # Demultiplexing creates a file for every metadata sample, including
        # samples with zero assigned reads. Empty files have nothing to align.
        if fastq_path.stat().st_size == 0:
            continue
        output_bam = bam_dir / f"{fastq_path.stem}.sorted.bam"
        align_fastq(
            fastq_path,
            reference_path,
            output_bam,
            bwa=bwa,
            samtools=samtools,
            threads=threads,
        )
        bam_paths.append(output_bam)
    return bam_paths


def read_single_reference(reference_path: Path) -> tuple[str, str]:
    """Read one FASTA record and return its name and uppercase sequence."""

    name: str | None = None
    sequence_parts: list[str] = []
    with open_text(reference_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    raise ValueError("Variant calling currently requires a single-record FASTA")
                name = line[1:].split()[0]
                if not name:
                    raise ValueError(f"Reference FASTA line {line_number} has an empty name")
            else:
                if name is None:
                    raise ValueError("Reference FASTA sequence appears before its header")
                sequence_parts.append(line.upper())

    sequence = "".join(sequence_parts)
    if name is None or not sequence:
        raise ValueError("Reference FASTA contains no sequence")
    if not DNA_PATTERN.fullmatch(sequence):
        raise ValueError("Reference FASTA contains an invalid nucleotide")
    return name, sequence


def call_variants(
    bam_path: Path,
    reference_path: Path,
    *,
    min_depth: int = 20,
    min_alt_count: int = 3,
    min_alt_fraction: float = 0.20,
    min_base_quality: int = 20,
    min_mapping_quality: int = 20,
) -> list[VariantCall]:
    """Call simple SNVs by comparing pileup bases directly with the reference.

    This transparent caller is intended for a short teaching/demo amplicon. It
    does not replace a production caller that models uncertainty, indels,
    strand bias, local haplotypes, and other sources of sequencing error.
    """

    import pysam

    if min_depth < 1 or min_alt_count < 1:
        raise ValueError("Depth and alternate-count thresholds must be positive")
    if not 0 < min_alt_fraction <= 1:
        raise ValueError("min_alt_fraction must be greater than 0 and no more than 1")

    reference_name, reference_sequence = read_single_reference(reference_path)
    calls: list[VariantCall] = []
    with pysam.AlignmentFile(bam_path, "rb") as alignment:
        if reference_name not in alignment.references:
            raise ValueError(
                f"Reference {reference_name!r} is absent from {bam_path.name}"
            )
        for column in alignment.pileup(
            reference_name,
            stepper="samtools",
            min_base_quality=min_base_quality,
            min_mapping_quality=min_mapping_quality,
            truncate=True,
        ):
            reference_index = column.reference_pos
            if reference_index >= len(reference_sequence):
                continue
            counts = {base: 0 for base in "ACGT"}
            for pileup_read in column.pileups:
                if pileup_read.is_del or pileup_read.is_refskip:
                    continue
                query_position = pileup_read.query_position
                if query_position is None:
                    continue
                base = pileup_read.alignment.query_sequence[query_position].upper()
                if base in counts:
                    counts[base] += 1

            depth = sum(counts.values())
            if depth < min_depth:
                continue
            reference_base = reference_sequence[reference_index]
            alternate_base, alternate_count = max(
                ((base, count) for base, count in counts.items() if base != reference_base),
                key=lambda item: item[1],
            )
            alternate_fraction = alternate_count / depth
            # A call must have adequate total coverage, a minimum number of
            # alternate-supporting reads, and a sufficient allele fraction.
            # These simple evidence filters are configurable for each experiment.
            if (
                alternate_count >= min_alt_count
                and alternate_fraction >= min_alt_fraction
            ):
                calls.append(
                    VariantCall(
                        contig=reference_name,
                        # pysam pileup coordinates are zero-based; variant tables
                        # conventionally report biological positions as one-based.
                        position=reference_index + 1,
                        reference=reference_base,
                        alternate=alternate_base,
                        depth=depth,
                        reference_count=counts[reference_base],
                        alternate_count=alternate_count,
                        alternate_fraction=alternate_fraction,
                    )
                )
    return calls


def write_variant_table(
    bam_paths: list[Path],
    metadata_path: Path,
    reference_path: Path,
    output_path: Path,
    *,
    sample_column: str = "Name",
    group_column: str = "Color",
    barcode_column: str = "Barcode",
    min_depth: int = 20,
    min_alt_count: int = 3,
    min_alt_fraction: float = 0.20,
) -> int:
    """Call variants for all samples and write a structured TSV report."""

    samples = read_samples(metadata_path, sample_column, group_column, barcode_column)
    samples_by_name = {sample.name: sample for sample in samples.values()}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    call_count = 0
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample",
                "group",
                "contig",
                "position_1_based",
                "reference",
                "alternate",
                "depth",
                "reference_count",
                "alternate_count",
                "alternate_fraction",
            ]
        )
        for bam_path in sorted(bam_paths):
            sample_name = bam_path.name.removesuffix(".sorted.bam")
            sample = samples_by_name.get(sample_name)
            if sample is None:
                raise ValueError(f"No metadata found for BAM sample {sample_name!r}")
            for call in call_variants(
                bam_path,
                reference_path,
                min_depth=min_depth,
                min_alt_count=min_alt_count,
                min_alt_fraction=min_alt_fraction,
            ):
                writer.writerow(
                    [
                        sample.name,
                        sample.group,
                        call.contig,
                        call.position,
                        call.reference,
                        call.alternate,
                        call.depth,
                        call.reference_count,
                        call.alternate_count,
                        f"{call.alternate_fraction:.4f}",
                    ]
                )
                call_count += 1
    return call_count


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Demultiplex barcoded amplicon reads and optionally align them and call SNVs."
        )
    )
    parser.add_argument("--reads", required=True, type=Path, help="Input FASTQ or FASTQ.GZ")
    parser.add_argument(
        "--metadata", required=True, type=Path, help="Tab-delimited sample metadata"
    )
    parser.add_argument("--output", required=True, type=Path, help="Output directory")
    parser.add_argument(
        "--reference",
        type=Path,
        help="Reference FASTA; when supplied, generate sorted and indexed BAM files",
    )
    parser.add_argument(
        "--threads", type=int, default=1, help="Threads used by BWA and samtools (default: 1)"
    )
    parser.add_argument("--min-depth", type=int, default=20)
    parser.add_argument("--min-alt-count", type=int, default=3)
    parser.add_argument("--min-alt-fraction", type=float, default=0.20)
    parser.add_argument("--sample-column", default="Name")
    parser.add_argument("--group-column", default="Color")
    parser.add_argument("--barcode-column", default="Barcode")
    parser.add_argument(
        "--trimming-mode",
        choices=("phred", "classroom-df", "none"),
        default="phred",
        help="Read-tail trimming method (default: phred)",
    )
    parser.add_argument(
        "--trim-quality",
        type=int,
        default=20,
        help="Minimum retained terminal Phred+33 quality (default: 20)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    counts = demultiplex_fastq(
        args.reads,
        args.metadata,
        args.output / "fastqs",
        sample_column=args.sample_column,
        group_column=args.group_column,
        barcode_column=args.barcode_column,
        trimming_mode=args.trimming_mode,
        trim_quality=args.trim_quality,
    )
    summary_path = args.output / "demultiplexing_summary.tsv"
    write_summary(counts, summary_path)
    assigned = sum(count for sample, count in counts.items() if sample != "__unmatched__")
    print(f"Assigned reads: {assigned}")
    print(f"Unmatched reads: {counts['__unmatched__']}")
    print(f"Summary: {summary_path}")
    if args.reference is not None:
        bam_paths = align_sample_fastqs(
            args.output / "fastqs",
            args.reference,
            args.output / "bams",
            threads=args.threads,
        )
        print(f"Aligned samples: {len(bam_paths)}")
        print(f"BAM directory: {args.output / 'bams'}")
        variant_path = args.output / "variants.tsv"
        variant_count = write_variant_table(
            bam_paths,
            args.metadata,
            args.reference,
            variant_path,
            sample_column=args.sample_column,
            group_column=args.group_column,
            barcode_column=args.barcode_column,
            min_depth=args.min_depth,
            min_alt_count=args.min_alt_count,
            min_alt_fraction=args.min_alt_fraction,
        )
        print(f"Variant calls: {variant_count}")
        print(f"Variant table: {variant_path}")


if __name__ == "__main__":
    main()
