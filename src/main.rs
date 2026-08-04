use kmer_counter::index::{export_registry_to_txt, export_table_to_txt};
use kmer_counter::{persistence, process_folder, FilterConfig};
// use kmer_counter::{kmers, persistence, process_folder, FilterConfig};
use std::path::PathBuf;

use clap::Parser;

/// Indexador de k-meros canónicos a partir de archivos FASTA.
/// Ejemplos:
///   kmer_indexer ./data
///   kmer_indexer ./data -k 21 --min-gc 0.2 --max-gc 0.8
///   kmer_indexer ./data --min-entropy 1.0 --no-txt -o ./resultados
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Carpeta con archivos FASTA (.fasta, .faa, .fna) o archivo unico
    path: PathBuf,
 
    /// Tamaño de k-mero entre 1 y 32 limite físico (2 bits por base, 32 bases = 64 bits)
    #[arg(short = 'k', long, default_value_t = 13)]
    k: usize,
 
    /// Contenido GC mínimo aceptado (0.0 - 1.0). Si no se pasa, este filtro queda desactivado.
    #[arg(long)]
    min_gc: Option<f32>,
 
    /// Contenido GC máximo aceptado (0.0 - 1.0).
    #[arg(long)]
    max_gc: Option<f32>,
 
    /// Entropía de Shannon mínima aceptada (0.0 - 2.0). Descarta k-meros
    /// de baja diversidad interna
    #[arg(long)]
    min_entropy: Option<f32>,
 
    /// No generar los reportes legibles kmer_table.txt / file_registry.txt.
    #[arg(long)]
    no_txt: bool,
 
    /// No generar los archivos binarios kmer_table.bin / file_registry.bin.
    #[arg(long)]
    no_bin: bool,
 
    /// Carpeta donde escribir los archivos de salida (reportes y binarios).
    #[arg(short, long, default_value = ".")]
    output: PathBuf,
}
 
fn main() -> std::io::Result<()> {
    let cli = Cli::parse();
 
    if cli.k == 0 || cli.k > 32 {
        eprintln!("Error: -k debe estar entre 1 y 32 (recibido: {})", cli.k);
        std::process::exit(1);
    }
 
    let filter_config = FilterConfig {
        min_gc: cli.min_gc,
        max_gc: cli.max_gc,
        min_entropy: cli.min_entropy,
    };
 
    let result = process_folder(&cli.path, cli.k, &filter_config)?;
 
    if result.stats.total_files == 0 {
        println!(
            "No .fasta/.faa/.fna files were found in: {}",
            cli.path.display()
        );
        return Ok(());
    }
 
    println!("=== Summary ===");
    println!("Files analyzed: {}", result.stats.total_files);
    println!("Total sequences: {}", result.stats.total_sequences);
    println!(
        "Sum of distinct k-mers per sequence (not deduplicated across files): {}",
        result.stats.total_distinct_kmers_per_sequence
    );
    println!(
        "Distinct canonical k-mers across the whole dataset: {}",
        result.table.len()
    );
 
    // --- Bloque temporal de depuración
    /* 
    println!("\n=== Debug: primeros 5 k-meros canónicos ===");
    for entry in result.table.iter().take(5) {
        let canonical_kmer = entry.key();
        let data = entry.value();
 
        let decoded = kmers::decode_kmer(*canonical_kmer, cli.k);
        let paths: Vec<String> = data
            .files
            .iter()
            .filter_map(|&id| result.registry.path_for(id))
            .map(|p| p.to_string_lossy().into_owned())
            .collect();
 
        println!(
            "k-mero: {} (count={}) -> presente en: {:?}",
            decoded, data.count, paths
        );
    }*/
 
    std::fs::create_dir_all(&cli.output)?;
 
    if !cli.no_txt {
        export_table_to_txt(&result.table, cli.k, &cli.output.join("kmer_table.txt"))?;
        export_registry_to_txt(&result.registry, &cli.output.join("file_registry.txt"))?;
        println!(
            "Exported kmer_table.txt and file_registry.txt to {}",
            cli.output.display()
        );
    }
 
    if !cli.no_bin {
        persistence::save_table(&result.table, &cli.output.join("kmer_table.bin"))?;
        persistence::save_registry(&result.registry, &cli.output.join("file_registry.bin"))?;
        println!(
            "Saved kmer_table.bin and file_registry.bin to {}",
            cli.output.display()
        );
    }
 
    Ok(())
}