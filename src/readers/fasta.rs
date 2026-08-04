use std::fs::File;
use std::io::{BufReader};
use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::fs;
use crate::models::FastaRecord;

pub struct FastaReader{
    reader: BufReader<File>,
    pending_header: Option<String> // hader ya leido
}

impl FastaReader {
    pub fn new(path: &str) -> Result<Self, std::io::Error> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        Ok(Self{
            reader,
            pending_header: None,
        })
        
    }

    pub fn next_line(&mut self) -> Result<Option<String>, std::io::Error> {
        //todo!()
        let mut line = String::new();
        let bytes_read = self.reader.read_line(&mut line)?;

        if bytes_read == 0 {
            return Ok(None);
        }

        if line.ends_with('\n') {
            line.pop();
            if line.ends_with('\r') {
                line.pop();
            }
        }

        Ok(Some(line))
        //self.reader.read_line();

    }

    pub fn next_record(&mut self) -> Result<Option<FastaRecord>, std::io::Error> {
        // conseguir necabezado actual
        let header_line = match self.pending_header.take() {
            Some(h) => h, //guardada
            None => {
                loop {
                    match self.next_line()? {
                        None => return Ok(None),
                        Some(line) => {
                            if line.starts_with('>') {
                                break line;
                            }
                        }
                        
                    }
                }
            }
        };

        // parsear header
        let no_greather = &header_line[1..];
        let mut part = no_greather.splitn(2, ' ');
        let id = part.next().unwrap_or("").to_string();
        let description = part.next().unwrap_or("").to_string();

        // leer secuencia hasta siguiente elemnto
        let mut sequence = String::new();
        loop {
            match self.next_line()? {
                None => {
                    self.pending_header = None;
                    break;
                }
                Some(line) => {
                    if line.starts_with('>') {
                        self.pending_header = Some(line); //proximo record
                        break;
                    }
                    sequence.push_str(&line);
                }
            }
        }
        Ok(Some(FastaRecord {header: header_line, sequence: sequence, id, description: description }))
    }

}


    pub fn lookfor_fasta(base: &Path) -> std::io::Result<Vec<PathBuf>> {
        let mut fasta_founds = Vec::new();
        recursive_looking(base, &mut fasta_founds)?;
        Ok(fasta_founds)
    }

    fn recursive_looking(path: &Path, fasta_founds: &mut Vec<PathBuf>) -> std::io::Result<()> {
        if path.is_file() {
            if has_valid_extension(path) {
                fasta_founds.push(path.to_path_buf());
            }
            return Ok(());
        }
        if path.is_dir() {
            for entry in fs::read_dir(path)?{
                let entry = entry?;
                let road = entry.path();
                recursive_looking(&road, fasta_founds)?;
            }
        }
        Ok(())
    }

    fn has_valid_extension(path: &Path) -> bool {
        match path.extension().and_then(|ext| ext.to_str()) {
            Some("fasta") |Some("faa") |Some("fna") => true,
            _ => false,
        }
    }
