//! Librería núcleo de indexación de k-meros

pub mod readers;
pub mod models;
pub mod kmers;
pub mod index;
pub mod persistence;

// Re-exportar tipos importantes para que los usuarios de la librería puedan acceder a ellos directamente desde `kmer_counter::...`
pub use index::{FileEntry, FileRegistry, FilterConfig, KmerData, KmerTable};

use std::io;
use std::path::{Path, PathBuf};

use dashmap::DashMap;
use rayon::prelude::*;

use index::register_kmer;
use kmers::{count_kmers, get_canonical};
use readers::fasta::{lookfor_fasta, FastaReader};


#[derive(Debug, Clone, Default)]
pub struct ProcessingStats {
    pub total_files: usize,
    pub total_sequences: usize,
    pub total_distinct_kmers_per_sequence: usize,
}

/// Resultado completo de procesar una carpeta de archivos FASTA
pub struct ProcessResult {
    pub table: KmerTable,
    pub registry: FileRegistry,
    pub stats: ProcessingStats,
}

struct FileProcessResult {
    archivo_id: u32,
    header: Option<String>,
    sequences: usize,
    distinct_kmers_this_file: usize,
}
 
/// Procesa un solo archivo FASTA: abre su propio `FastaReader`
fn process_single_file(
    archivo_id: u32,
    path: &Path,
    k: usize,
    filters: &FilterConfig,
    table: &KmerTable,
) -> io::Result<FileProcessResult> {
    let path_str = path.to_str().ok_or_else(|| {
        io::Error::new(io::ErrorKind::InvalidInput, "la ruta no es UTF-8 válido")
    })?;
    let mut reader = FastaReader::new(path_str)?;
 
    let mut header = None;
    let mut sequences = 0;
    let mut distinct_kmers_this_file = 0;
 
    while let Some(record) = reader.next_record()? {
        if sequences == 0 {
            header = Some(record.id.clone());
        }
        sequences += 1;
 
        let kmer_counts = count_kmers(record.sequence.as_bytes(), k);
        distinct_kmers_this_file += kmer_counts.len();
 
        for &raw_kmer in kmer_counts.keys() {
            let canonical = get_canonical(raw_kmer, k);
            register_kmer(table, canonical, archivo_id, k, filters);
        }
    }
 
    Ok(FileProcessResult {
        archivo_id,
        header,
        sequences,
        distinct_kmers_this_file,
    })
}


/// Punto de entrada principal de la librería: recorre todos los archivos
/// FASTA bajo `root`, cuenta y canonicaliza k-meros de tamaño `k`, aplica
/// `filters`, y devuelve la tabla resultante junto con el registro de
/// archivos y estadísticas de la corrida.
/// El trabajo por archivo se hace en paralelo (via `rayon`)
pub fn process_folder(
    root: &Path,
    k: usize,
    filters: &FilterConfig,
) -> io::Result<ProcessResult> {
    let files = lookfor_fasta(root)?;
 
    // --- Fase 1: secuencial ---
    let mut registry = FileRegistry::new();
    let files_with_ids: Vec<(u32, PathBuf)> = files
        .into_iter()
        .map(|path| {
            let id = registry.register(path.clone());
            (id, path)
        })
        .collect();
 
    let table: KmerTable = DashMap::new();
 
    // --- Fase 2: paralela ---
    let results: Vec<io::Result<FileProcessResult>> = files_with_ids
        .par_iter()
        .map(|(archivo_id, path)| process_single_file(*archivo_id, path, k, filters, &table))
        .collect();
 
    // --- Fase 3: secuencial ---
    let mut total_files = 0;
    let mut total_sequences = 0;
    let mut total_distinct_kmers_per_sequence = 0;
 
    for result in results {
        let file_result = result?;
        if let Some(header) = file_result.header {
            registry.set_header(file_result.archivo_id, header);
        }
        total_files += 1;
        total_sequences += file_result.sequences;
        total_distinct_kmers_per_sequence += file_result.distinct_kmers_this_file;
    }
 
    Ok(ProcessResult {
        table,
        registry,
        stats: ProcessingStats {
            total_files,
            total_sequences,
            total_distinct_kmers_per_sequence,
        },
    })
}