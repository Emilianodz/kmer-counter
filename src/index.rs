use dashmap::mapref::entry::Entry;
use dashmap::DashMap;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

use crate::kmers::{decode_kmer, gc_content, shannon_entropy};

/// Configuración de filtros fisicoquímicos, todos opcionales.
///
/// `None` en cualquier campo significa "este filtro está apagado". Con
#[derive(Debug, Clone, Default)]
pub struct FilterConfig {
    pub min_gc: Option<f32>,
    pub max_gc: Option<f32>,
    pub min_entropy: Option<f32>,
}

/// Evalúa si un k-mero pasa la configuración de filtros dada.
///
fn passes_filters(gc: Option<f32>, entropy: Option<f32>, config: &FilterConfig) -> bool {
    if let (Some(min), Some(gc_value)) = (config.min_gc, gc) {
        if gc_value < min {
            return false;
        }
    }
    if let (Some(max), Some(gc_value)) = (config.max_gc, gc) {
        if gc_value > max {
            return false;
        }
    }
    if let (Some(min), Some(entropy_value)) = (config.min_entropy, entropy) {
        if entropy_value < min {
            return false;
        }
    }
    true
}

/// Datos acumulados para un k-mero canónico: cuántas veces se vio en total,
/// en qué archivos, y sus propiedades fisicoquímicas (calculadas una sola
/// vez, la primera vez que se ve el k-mero -- nunca se recalculan)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerData {
    pub count: u32,
    pub files: HashSet<u32>,
    pub gc_content: Option<f32>,
    pub entropy: Option<f32>,
}

/// Tabla global: k-mero canónico -> sus datos acumulados
pub type KmerTable = DashMap<u64, KmerData>;

/// Registra que `kmer` (ya canonicalizado) apareció en el archivo
/// identificado por `archivo_id`, aplicando `config` si el k-mero es nuevo.
/// Devuelve `true` si el k-mero quedó registrado (ya existía, o era nuevo y
/// pasó los filtros), `false` si era nuevo pero fue descartado por no pasar
/// los filtros configurados.

pub fn register_kmer(
    table: &KmerTable,
    kmer: u64,
    archivo_id: u32,
    k: usize,
    config: &FilterConfig,
) -> bool {
    match table.entry(kmer) {
        Entry::Occupied(mut entry) => {
            let data = entry.get_mut();
            data.count += 1;
            data.files.insert(archivo_id);
            true
        }
        Entry::Vacant(entry) => {
            let gc = if config.min_gc.is_some() || config.max_gc.is_some() {
                Some(gc_content(kmer, k))
            } else {
                None
            };
            let entropy = if config.min_entropy.is_some() {
                Some(shannon_entropy(kmer, k))
            } else {
                None
            };

            if !passes_filters(gc, entropy, config) {
                return false;
            }

            let mut files = HashSet::new();
            files.insert(archivo_id);
            entry.insert(KmerData {
                count: 1,
                files,
                gc_content: gc,
                entropy,
            });
            true
        }
    }
}

/// Un archivo registrado: su ruta real, y el header
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileEntry {
    pub path: PathBuf,
    pub header_id: Option<String>,
}

/// Tabla auxiliar para traducir un `archivo_id: u32` de vuelta a su
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct FileRegistry {
    entries: Vec<FileEntry>,
}

impl FileRegistry {
    pub fn new() -> Self {
        FileRegistry { entries: Vec::new() }
    }

    /// Registra un archivo y devuelve el `archivo_id` (u32)
    pub fn register(&mut self, path: PathBuf) -> u32 {
        let id = self.entries.len() as u32;
        self.entries.push(FileEntry { path, header_id: None });
        id
    }

    pub fn set_header(&mut self, archivo_id: u32, header_id: String) {
        if let Some(entry) = self.entries.get_mut(archivo_id as usize) {
            entry.header_id = Some(header_id);
        }
    }

    pub fn path_for(&self, archivo_id: u32) -> Option<&PathBuf> {
        self.entries.get(archivo_id as usize).map(|entry| &entry.path)
    }

    pub fn header_for(&self, archivo_id: u32) -> Option<&str> {
        self.entries
            .get(archivo_id as usize)
            .and_then(|entry| entry.header_id.as_deref())
    }

    pub fn iter(&self) -> impl Iterator<Item = (u32, &FileEntry)> {
        self.entries
            .iter()
            .enumerate()
            .map(|(id, entry)| (id as u32, entry))
    }
}

/// Vuelca la tabla de k-meros a un archivo de texto, una línea por k-mero
/// canónico distinto, ordenado por su valor numérico 
/// orden: k-mer(AAAAAAAACACAG), u64(274), count(1), gc:content(0.23=23%), entrpy(max 2), archivo_id(1) 
pub fn export_table_to_txt(table: &KmerTable, k: usize, path: &Path) -> io::Result<()> {

    let mut entries: Vec<(u64, KmerData)> = table
        .iter()
        .map(|entry| (*entry.key(), entry.value().clone()))
        .collect();
    entries.sort_by_key(|(kmer, _)| *kmer);

    let mut file = File::create(path)?;
    for (kmer, data) in entries {
        let mut ids: Vec<u32> = data.files.into_iter().collect();
        ids.sort();
        let ids_str = ids
            .iter()
            .map(|id| id.to_string())
            .collect::<Vec<_>>()
            .join(",");

        let gc_str = data
            .gc_content
            .map(|v| format!("{:.3}", v))
            .unwrap_or_else(|| "NA".to_string());
        let entropy_str = data
            .entropy
            .map(|v| format!("{:.3}", v))
            .unwrap_or_else(|| "NA".to_string());

        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}",
            decode_kmer(kmer, k),
            kmer,
            data.count,
            gc_str,
            entropy_str,
            ids_str
        )?;
    }

    Ok(())
}

pub fn export_registry_to_txt(registry: &FileRegistry, path: &Path) -> io::Result<()> {
    let mut file = File::create(path)?;
    for (id, entry) in registry.iter() {
        let filename = entry
            .path
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_else(|| entry.path.display().to_string());
        let header = entry.header_id.as_deref().unwrap_or("NA");
 
        writeln!(file, "{}\t{}\t{}", id, filename, header)?;
    }
    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;
 
    #[test]
    fn register_kmer_accumulates_count_and_files_without_filters() {
        let table: KmerTable = DashMap::new();
        let config = FilterConfig::default(); // todo en None -> nada se filtra
 
        register_kmer(&table, 42, 0, 3, &config);
        register_kmer(&table, 42, 1, 3, &config);
 
        let data = table.get(&42).unwrap();
        assert_eq!(data.count, 2);
        assert_eq!(data.files.len(), 2);
        assert!(data.files.contains(&0));
        assert!(data.files.contains(&1));
        // Sin filtros activos, gc/entropy no se calculan -- deben quedar None.
        assert_eq!(data.gc_content, None);
        assert_eq!(data.entropy, None);
    }
 
    #[test]
    fn register_kmer_does_not_duplicate_same_file() {
        let table: KmerTable = DashMap::new();
        let config = FilterConfig::default();
 
        register_kmer(&table, 42, 0, 3, &config);
        register_kmer(&table, 42, 0, 3, &config); // mismo archivo_id otra vez
 
        let data = table.get(&42).unwrap();
        assert_eq!(data.files.len(), 1); // no se duplicó
        assert_eq!(data.count, 2); // el conteo SÍ sube -- son cosas distintas
    }
 
    #[test]
    fn register_kmer_rejects_kmer_below_min_gc() {
        let table: KmerTable = DashMap::new();
        let config = FilterConfig {
            min_gc: Some(0.5),
            ..Default::default()
        };
 
        // "AAA" (code 0) tiene 0% GC -> no debe pasar min_gc = 0.5
        let accepted = register_kmer(&table, 0, 0, 3, &config);
 
        assert!(!accepted);
        assert!(table.get(&0).is_none()); // nunca se insertó
    }
 
    #[test]
    fn register_kmer_accepts_kmer_that_passes_filters_and_computes_metrics() {
        let table: KmerTable = DashMap::new();
        let config = FilterConfig {
            min_gc: Some(0.0),
            ..Default::default()
        };
 
        let accepted = register_kmer(&table, 0, 0, 3, &config); // "AAA", 0% GC
 
        assert!(accepted);
        let data = table.get(&0).unwrap();
        assert_eq!(data.gc_content, Some(0.0)); // sí se calculó, porque min_gc estaba activo
        assert_eq!(data.entropy, None); // min_entropy no estaba activo, no se calculó
    }
 
#[cfg(test)]
mod tests {
    use super::*;
 
    #[test]
    fn register_kmer_accumulates_count_and_files_without_filters() {
        let table: KmerTable = DashMap::new();
        let config = FilterConfig::default(); // todo en None -> nada se filtra
 
        register_kmer(&table, 42, 0, 3, &config);
        register_kmer(&table, 42, 1, 3, &config);
 
        let data = table.get(&42).unwrap();
        assert_eq!(data.count, 2);
        assert_eq!(data.files.len(), 2);
        assert!(data.files.contains(&0));
        assert!(data.files.contains(&1));
        // Sin filtros activos, gc/entropy no se calculan -- deben quedar None.
        assert_eq!(data.gc_content, None);
        assert_eq!(data.entropy, None);
    }
 
    #[test]
    fn register_kmer_does_not_duplicate_same_file() {
        let table: KmerTable = DashMap::new();
        let config = FilterConfig::default();
 
        register_kmer(&table, 42, 0, 3, &config);
        register_kmer(&table, 42, 0, 3, &config); // mismo archivo_id otra vez
 
        let data = table.get(&42).unwrap();
        assert_eq!(data.files.len(), 1); // no se duplicó
        assert_eq!(data.count, 2); // el conteo SÍ sube -- son cosas distintas
    }
 
    #[test]
    fn register_kmer_rejects_kmer_below_min_gc() {
        let table: KmerTable = DashMap::new();
        let config = FilterConfig {
            min_gc: Some(0.5),
            ..Default::default()
        };
 
        // "AAA" (code 0) tiene 0% GC -> no debe pasar min_gc = 0.5
        let accepted = register_kmer(&table, 0, 0, 3, &config);
 
        assert!(!accepted);
        assert!(table.get(&0).is_none()); // nunca se insertó
    }
 
    #[test]
    fn register_kmer_accepts_kmer_that_passes_filters_and_computes_metrics() {
        let table: KmerTable = DashMap::new();
        let config = FilterConfig {
            min_gc: Some(0.0),
            ..Default::default()
        };
 
        let accepted = register_kmer(&table, 0, 0, 3, &config); // "AAA", 0% GC
 
        assert!(accepted);
        let data = table.get(&0).unwrap();
        assert_eq!(data.gc_content, Some(0.0)); // sí se calculó, porque min_gc estaba activo
        assert_eq!(data.entropy, None); // min_entropy no estaba activo, no se calculó
    }
 
    #[test]
    fn export_table_to_txt_writes_expected_lines() {
        let table: KmerTable = DashMap::new();
        let config = FilterConfig::default();
        // "AAA" -> code 0 (para k=3); ya canonicalizado a mano para no
        // depender de kmers::get_canonical en este test.
        register_kmer(&table, 0, 0, 3, &config);
        register_kmer(&table, 0, 1, 3, &config);
 
        let tmp_path = std::env::temp_dir().join("test_export_table.txt");
        export_table_to_txt(&table, 3, &tmp_path).unwrap();
 
        let content = std::fs::read_to_string(&tmp_path).unwrap();
        // count=2, gc=NA, entropy=NA (filtros apagados), archivos=0,1
        assert_eq!(content.trim(), "AAA\t0\t2\tNA\tNA\t0,1");
 
        std::fs::remove_file(&tmp_path).ok(); // limpieza, ignora error si ya no existe
    }
 
    #[test]
    fn export_registry_to_txt_writes_filename_and_header() {
        let mut registry = FileRegistry::new();
        let id_a = registry.register(PathBuf::from("/ruta/completa/a.fasta"));
        registry.set_header(id_a, "chr1".to_string());
        registry.register(PathBuf::from("/otra/ruta/b.fasta")); // sin header todavía
 
        let tmp_path = std::env::temp_dir().join("test_export_registry.txt");
        export_registry_to_txt(&registry, &tmp_path).unwrap();
 
        let content = std::fs::read_to_string(&tmp_path).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        // Nota: solo el NOMBRE del archivo, no la ruta completa.
        assert_eq!(lines[0], "0\ta.fasta\tchr1");
        assert_eq!(lines[1], "1\tb.fasta\tNA"); // header desconocido -> "NA"
 
        std::fs::remove_file(&tmp_path).ok();
    }
 
    #[test]
    fn file_registry_roundtrip() {
        let mut registry = FileRegistry::new();
        let id_a = registry.register(PathBuf::from("archivo_a.fasta"));
        let id_b = registry.register(PathBuf::from("archivo_b.fasta"));
 
        assert_eq!(id_a, 0);
        assert_eq!(id_b, 1);
        assert_eq!(registry.path_for(0), Some(&PathBuf::from("archivo_a.fasta")));
        assert_eq!(registry.path_for(1), Some(&PathBuf::from("archivo_b.fasta")));
        assert_eq!(registry.path_for(99), None);
    }
 
    #[test]
    fn file_registry_header_starts_none_and_can_be_set() {
        let mut registry = FileRegistry::new();
        let id = registry.register(PathBuf::from("genoma.fasta"));
 
        assert_eq!(registry.header_for(id), None); // todavía no se leyó el archivo
 
        registry.set_header(id, "NC_045512.2".to_string());
        assert_eq!(registry.header_for(id), Some("NC_045512.2"));
    }
}}