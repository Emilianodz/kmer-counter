# kmer_counter
An extremely fast and memory-efficient genomic indexing engine for large-scale whole-genome analysis.

Rust library for indexing canonical k-mers from FASTA files (one genome per
file, one record per file). Designed as a reusable engine: the core logic
lives in `lib.rs` and is exposed as a public API to be consumed by other
projects

## Project Status

- [x] FASTA file reading and traversal
- [x] k-mer encoding as `u64` (2 bits per nucleotide)
- [x] Canonicalization (reverse complement)
- [x] Concurrent table (`DashMap`) with occurrence counts, per-file presence,
      and optional physicochemical filters (GC content, Shannon entropy)
- [x] Binary serialization (`bincode`) for saving/loading results
- [x] Parallel processing with `rayon`
- [x] Command-line interface (`clap`) — no code modifications required to use the library

## Requirements

- Rust (Edition 2021) and Cargo.

## Usage

```bash
cargo test                # Run the unit test suite

cargo run -- ./data       # Process a directory using k=13 (default), no filters

cargo run -- ./data -k 21 --min-gc 0.2 --max-gc 0.8

cargo run -- ./data --min-entropy 1.0 --no-txt -o ./results

cargo run -- --help       # Show all available options
```

| Flag | Description |
|------|-------------|
| `path` (positional) | Directory containing the FASTA files to process |
| `-k, --k` | k-mer size (1–32, default: 13) |
| `--min-gc` / `--max-gc` | Acceptable GC content range (0.0–1.0) |
| `--min-entropy` | Minimum Shannon entropy (0.0–2.0) |
| `--no-txt` | Do not generate human-readable `.txt` reports |
| `--no-bin` | Do not generate binary `.bin` files |
| `-o, --output` | Output directory (default: current directory) |

## Output Files

Each execution generates:

- `kmer_table.txt` / `file_registry.txt`: human-readable plain-text reports,
  intended for inspection, validation, or discussion with collaborators.
- `kmer_table.bin` / `file_registry.bin`: compact binary files (`bincode`)
  that can be reloaded later without reprocessing the FASTA files.


## Author

Developed by Emiliano Díaz-López