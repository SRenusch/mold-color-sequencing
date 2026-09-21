# Mold Color Amplicon Sequencing Pipeline

A small, reproducible Python workflow for processing pooled, inline-barcoded
amplicon reads. It demultiplexes reads by sample, performs quality trimming,
aligns each sample to a reference, and reports supported single-nucleotide
variants (SNVs).

This repository is a portfolio-oriented refactoring of an introductory
bioinformatics assignment. The provenance and generation parameters of the
original classroom dataset were unavailable, and permission to redistribute
those files has not been confirmed. Therefore, this repository uses newly
generated fictional data with documented, known variants. This makes the
workflow reproducible without relying on or publishing the original dataset.

## Workflow

1. Validate the sample metadata and four-line FASTQ records.
2. Match the inline barcode at the beginning of each read to a sample.
3. Remove the barcode and trim consecutive low-quality terminal bases.
4. Align each non-empty sample FASTQ to the reference with BWA-MEM.
5. Sort and index the alignments as BAM files with samtools.
6. Compare aligned bases with the reference and report SNVs passing the chosen
   depth, alternate-read-count, and alternate-allele-fraction thresholds.

## Requirements

- Python 3.12 (the validated version)
- BWA 0.7.17 or a compatible release
- samtools 1.21 or a compatible release
- `pysam==0.24.1`

On macOS, BWA and samtools can be installed with Homebrew:

```bash
brew install bwa samtools
```

Create and activate a project-specific Python environment, then install the
Python dependency:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Reproduce the included example

The checked-in example contains eight fictional samples, six-base inline
barcodes, and 320 reads. Each sample contains a deliberately planted SNV at a
50% alternate-allele fraction. Three poor-quality terminal bases per read test
the default Phred-based trimming behavior.

Regenerate the example data deterministically:

```bash
python scripts/generate_example_data.py
```

Run the complete workflow:

```bash
python pipeline.py \
  --reads example_data/pooled_reads.fastq \
  --metadata example_data/samples.tsv \
  --reference example_data/reference.fa \
  --output results/example \
  --sample-column sample_id \
  --group-column phenotype \
  --barcode-column barcode \
  --threads 2
```

Expected console summary:

```text
Assigned reads: 320
Unmatched reads: 0
Aligned samples: 8
Variant calls: 8
```

The called sample, phenotype, position, reference allele, and alternate allele
should match `example_data/expected_variants.tsv`. Generated FASTQ, BAM, index,
summary, and variant files are written beneath `results/example/`; `results/`
is excluded from version control because these files are reproducible.

## Use different data

Supply a pooled FASTQ file, a single-record reference FASTA, and a tab-delimited
sample sheet. By default the metadata headings are `Name`, `Color`, and
`Barcode`; the three `--*-column` options let you use different headings.

```bash
python pipeline.py \
  --reads path/to/pooled_reads.fastq.gz \
  --metadata path/to/samples.tsv \
  --reference path/to/reference.fa \
  --output results/my_run
```

Use `python pipeline.py --help` to see the quality-trimming and variant-calling
threshold options.

## Tests

```bash
python -m unittest discover -s tests -v
```

The tests cover ordinary processing plus malformed FASTQ records, duplicate and
unknown barcodes, samples with no reads, Phred-tail trimming, missing reference
files, coordinate handling, and configurable low-frequency variant detection.

## Biological and technical limitations

- The included data are synthetic and demonstrate software behavior; they do
  not represent an actual mold experiment or establish biological causality.
- The workflow expects an inline sample barcode at the start of each read. It
  does not process separate Illumina index reads or unique molecular identifiers.
- Exact barcode matching is used. Barcode-error correction is not attempted.
- The reference-based caller reports only the most-supported SNV at a position.
  It is not a replacement for a production variant caller and does not report
  insertions, deletions, haplotypes, or statistical confidence.
- The current FASTA parser supports one reference record.
- The synthetic example is intentionally easy: uniform read length, clean
  mapping, balanced alleles, and no biological contamination or index hopping.

## Provenance and responsible use

The design originated as coursework. The pipeline was subsequently refactored,
documented, tested, and validated with AI-assisted development. The repository's
synthetic reference, sample metadata, reads, and expected results were created
for this demonstration. Original instructor-provided files and the original
submitted script are retained locally but excluded from version control unless
redistribution permission is confirmed.

## License

The repository code is available under the MIT License. See `LICENSE`.
