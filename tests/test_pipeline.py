import tempfile
import unittest
from pathlib import Path

import pysam

from pipeline import (
    align_sample_fastqs,
    call_variants,
    demultiplex_fastq,
    iter_fastq,
    read_samples,
    trim_degraded_tail,
    trim_phred_tail,
)


class PipelineTests(unittest.TestCase):
    def test_trim_degraded_tail_preserves_matching_lengths(self):
        sequence, quality = trim_degraded_tail("ACGTACGT", "IIIIIDDF")
        self.assertEqual(sequence, "ACGTA")
        self.assertEqual(quality, "IIIII")

    def test_phred_tail_trimming_uses_numeric_quality(self):
        sequence, quality = trim_phred_tail("ACGTACGT", "IIIII5+#", min_quality=20)
        self.assertEqual(sequence, "ACGTAC")
        self.assertEqual(quality, "IIIII5")

    def test_invalid_fastq_lengths_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.fastq"
            path.write_text("@read\nACGT\n+\nIII\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unequal sequence and quality"):
                list(iter_fastq(path))

    def test_incomplete_fastq_record_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "truncated.fastq"
            path.write_text("@read\nACGT\n+\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Incomplete FASTQ record 1"):
                list(iter_fastq(path))

    def test_duplicate_barcodes_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "samples.tsv"
            path.write_text(
                "Name\tColor\tBarcode\nA\tBlack\tAACGT\nB\tGreen\tAACGT\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "Duplicate barcode"):
                read_samples(path)

    def test_demultiplex_removes_barcode_from_sequence_and_quality(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = root / "samples.tsv"
            reads = root / "reads.fastq"
            output = root / "out"
            metadata.write_text(
                "Name\tColor\tBarcode\nSample A\tBlack\tAACGT\n",
                encoding="utf-8",
            )
            reads.write_text(
                "@matched\nAACGTACGTAC\n+\nIIIIIIIIIII\n"
                "@unmatched\nTTTTTACGTAC\n+\nIIIIIIIIIII\n",
                encoding="utf-8",
            )

            counts = demultiplex_fastq(
                reads, metadata, output, trimming_mode="none"
            )

            self.assertEqual(counts["Sample A"], 1)
            self.assertEqual(counts["__unmatched__"], 1)
            self.assertEqual(
                (output / "Sample_A.fastq").read_text(encoding="utf-8"),
                "@matched\nACGTAC\n+\nIIIIII\n",
            )

    def test_unknown_barcode_is_counted_without_misassignment(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = root / "samples.tsv"
            reads = root / "reads.fastq"
            output = root / "out"
            metadata.write_text(
                "Name\tColor\tBarcode\nSample A\tBlack\tAACGT\n",
                encoding="utf-8",
            )
            reads.write_text(
                "@unknown\nTTTTTACGTAC\n+\nIIIIIIIIIII\n",
                encoding="utf-8",
            )

            counts = demultiplex_fastq(
                reads, metadata, output, trimming_mode="none"
            )

            self.assertEqual(counts["Sample A"], 0)
            self.assertEqual(counts["__unmatched__"], 1)
            self.assertEqual((output / "Sample_A.fastq").stat().st_size, 0)

    def test_sample_with_no_reads_gets_empty_output_and_zero_count(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = root / "samples.tsv"
            reads = root / "reads.fastq"
            output = root / "out"
            metadata.write_text(
                "Name\tColor\tBarcode\n"
                "Observed\tBlack\tAACGT\n"
                "Missing\tGreen\tCCGTA\n",
                encoding="utf-8",
            )
            reads.write_text(
                "@observed\nAACGTACGT\n+\nIIIIIIIII\n",
                encoding="utf-8",
            )

            counts = demultiplex_fastq(
                reads, metadata, output, trimming_mode="none"
            )

            self.assertEqual(counts["Observed"], 1)
            self.assertEqual(counts["Missing"], 0)
            self.assertTrue((output / "Missing.fastq").is_file())
            self.assertEqual((output / "Missing.fastq").stat().st_size, 0)

    def test_alignment_rejects_missing_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(FileNotFoundError, "Reference FASTA not found"):
                align_sample_fastqs(
                    root / "fastqs", root / "missing.fa", root / "bams"
                )

    def test_variant_is_compared_with_reference_using_one_based_position(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference.fa"
            bam_path = root / "sample.sorted.bam"
            reference.write_text(">test\nACGTACGT\n", encoding="utf-8")
            header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "test", "LN": 8}]}
            with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
                for index in range(8):
                    read = pysam.AlignedSegment()
                    read.query_name = f"read{index}"
                    read.query_sequence = "ACGCACGT" if index < 4 else "ACGTACGT"
                    read.flag = 0
                    read.reference_id = 0
                    read.reference_start = 0
                    read.mapping_quality = 60
                    read.cigar = ((0, 8),)
                    read.query_qualities = pysam.qualitystring_to_array("IIIIIIII")
                    bam.write(read)
            pysam.index(str(bam_path))

            calls = call_variants(
                bam_path,
                reference,
                min_depth=8,
                min_alt_count=3,
                min_alt_fraction=0.20,
            )

            self.assertEqual(len(calls), 1)
            self.assertEqual(calls[0].position, 4)
            self.assertEqual(calls[0].reference, "T")
            self.assertEqual(calls[0].alternate, "C")
            self.assertEqual(calls[0].alternate_count, 4)
            self.assertEqual(calls[0].alternate_fraction, 0.5)

    def test_low_frequency_variant_respects_configured_threshold(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference.fa"
            bam_path = root / "sample.sorted.bam"
            reference.write_text(">test\nACGTACGT\n", encoding="utf-8")
            header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "test", "LN": 8}]}
            with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
                for index in range(20):
                    read = pysam.AlignedSegment()
                    read.query_name = f"read{index}"
                    read.query_sequence = "ACGCACGT" if index < 2 else "ACGTACGT"
                    read.flag = 0
                    read.reference_id = 0
                    read.reference_start = 0
                    read.mapping_quality = 60
                    read.cigar = ((0, 8),)
                    read.query_qualities = pysam.qualitystring_to_array("IIIIIIII")
                    bam.write(read)
            pysam.index(str(bam_path))

            strict_calls = call_variants(
                bam_path,
                reference,
                min_depth=20,
                min_alt_count=2,
                min_alt_fraction=0.20,
            )
            permissive_calls = call_variants(
                bam_path,
                reference,
                min_depth=20,
                min_alt_count=2,
                min_alt_fraction=0.05,
            )

            self.assertEqual(strict_calls, [])
            self.assertEqual(len(permissive_calls), 1)
            self.assertEqual(permissive_calls[0].position, 4)
            self.assertEqual(permissive_calls[0].alternate_count, 2)
            self.assertEqual(permissive_calls[0].alternate_fraction, 0.10)


if __name__ == "__main__":
    unittest.main()
