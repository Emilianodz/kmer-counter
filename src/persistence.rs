use std::fs;
use std::io;
use std:: path::Path;

use crate::index::{FileRegistry, KmerTable};

fn to_io_error<E: std::fmt::Display>(e: E) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, e.to_string())
}

pub fn save_table(table: &KmerTable, path: &Path) -> io::Result<()> {
    let bytes = bincode::serde::encode_to_vec(table, bincode::config::standard()).map_err(to_io_error)?;
    fs::write(path, bytes)
}

pub fn load_table(path: &Path) -> io::Result<KmerTable> {
    let bytes = fs::read(path)?;
    let (table, _bytes_read): (KmerTable, usize) =
        bincode::serde::decode_from_slice(&bytes, bincode::config::standard()).map_err(to_io_error)?;
    Ok(table)
}

pub fn save_registry(registry: &FileRegistry, path: &Path) -> io::Result<()> {
    let bytes = bincode::serde::encode_to_vec(registry, bincode::config::standard())
        .map_err(to_io_error)?;
    fs::write(path, bytes)
}

pub fn load_registry(path: &Path) -> io::Result<FileRegistry> {
    let bytes = fs::read(path)?;
    let (registry, _bytes_read): (FileRegistry, usize) =
        bincode::serde::decode_from_slice(&bytes, bincode::config::standard())
            .map_err(to_io_error)?;
    Ok(registry)
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::{register_kmer, FilterConfig};
    use dashmap::DashMap;
    use std::path::PathBuf;
 
    #[test]
    fn table_survives_roundtrip_through_disk() {
        let table: KmerTable = DashMap::new();
        let config = FilterConfig::default();
        register_kmer(&table, 42, 0, 3, &config);
        register_kmer(&table, 42, 1, 3, &config);
 
        let tmp_path = std::env::temp_dir().join("test_table_roundtrip.bin");
        save_table(&table, &tmp_path).unwrap();
        let loaded = load_table(&tmp_path).unwrap();
 
        let original = table.get(&42).unwrap();
        let restored = loaded.get(&42).unwrap();
        assert_eq!(original.count, restored.count);
        assert_eq!(original.files, restored.files);
 
        std::fs::remove_file(&tmp_path).ok();
    }
 
    #[test]
    fn registry_survives_roundtrip_through_disk() {
        let mut registry = FileRegistry::new();
        let id = registry.register(PathBuf::from("genoma_a.fasta"));
        registry.set_header(id, "chr1_ejemplo".to_string());
 
        let tmp_path = std::env::temp_dir().join("test_registry_roundtrip.bin");
        save_registry(&registry, &tmp_path).unwrap();
        let loaded = load_registry(&tmp_path).unwrap();
 
        assert_eq!(loaded.path_for(id), registry.path_for(id));
        assert_eq!(loaded.header_for(id), Some("chr1_ejemplo"));
 
        std::fs::remove_file(&tmp_path).ok();
    }
}