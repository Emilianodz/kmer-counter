use std::collections::HashMap;


// Representación interna: cada k-mero se codifica como `u64`, usando 2 bits
// por base (A=00, C=01, G=10, T=11). Esto limita k a un máximo de 32
/// Codifica una ventana de bases (`&[u8]`) en un único `u64`.
/// 
/// 
pub fn encode_kmer(bases: &[u8]) -> Option<u64> {
    let mut code: u64 = 0;

    for &base in bases {
        let bits: u64 = match base {
            b'A' => 0b00,
            b'C' => 0b01,
            b'G' => 0b10,
            b'T' => 0b11,
            _ => return None, // early return: cualquier base inválida invalida toda la ventana
        };
        code = (code << 2) | bits;
    }

    Some(code)
}

/// Decodifica un `u64` de vuelta a su representación de bases (String)
/// Solo para depuración
/// No forma parte del pipeline de producción
pub fn decode_kmer(code: u64, k: usize) -> String {
    let mut bases = Vec::with_capacity(k);

    for i in 0..k {
        let shift = (k - 1 - i) * 2;
        let bits = (code >> shift) & 0b11;
        let base = match bits {
            0b00 => b'A',
            0b01 => b'C',
            0b10 => b'G',
            0b11 => b'T',
            _ => unreachable!("una máscara de 2 bits solo puede dar 0..=3"),
        };
        bases.push(base);
    }

    String::from_utf8(bases).unwrap()
}

/// Cuenta la frecuencia de cada k-mero válido en una secuencia.
pub fn count_kmers(sequence: &[u8], k: usize) -> HashMap<u64, usize> {
    let mut counts = HashMap::new();

    for window in sequence.windows(k) {
        let Some(code) = encode_kmer(window) else {
            continue;
        };
        // entry(...).or_insert(0) evita el patrón "check-then-insert":
        // busca la entrada, y si no existe, la crea con valor 0
        *counts.entry(code).or_insert(0) += 1;
    }

    counts
}

/// Calcula el complemento reverso de un k-mero codificado, operando
pub fn reverse_complement(code: u64, k: usize) -> u64 {
    let mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let complemented = !code & mask;

    let mut reversed: u64 = 0;
    let mut remaining = complemented;
    for _ in 0..k {
        let base_bits = remaining & 0b11;
        reversed = (reversed << 2) | base_bits;
        remaining >>= 2;
    }

    reversed
}

pub fn get_canonical(code: u64, k: usize) -> u64 {
    let rc = reverse_complement(code, k);
    code.min(rc)
}

pub fn gc_content(code: u64, k: usize) -> f32 {
    let xor = code ^ (code >> 1);
    const EVEN_BITS_MASK: u64 = 0x5555_5555_5555_5555;
    let gc_bits = xor & EVEN_BITS_MASK;
    let gc_count = gc_bits.count_ones();
    gc_count as f32 / k as f32
}

pub fn shannon_entropy(code: u64, k: usize) -> f32 {
    let mut counts = [0u32; 4]; // índices 0,1,2,3 = A,C,G,T (mismo mapeo que encode_kmer)
    let mut remaining = code;
    for _ in 0..k {
        let base = (remaining & 0b11) as usize;
        counts[base] += 1;
        remaining >>= 2;
    }

    let k_f = k as f32;
    let mut entropy = 0.0f32;
    for &count in counts.iter() {
        if count > 0 {
            let p = count as f32 / k_f;
            entropy -= p * p.log2();
        }
    }
    entropy
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_decode_roundtrip() {
        let bases = b"ATGC";
        let code = encode_kmer(bases).expect("bases válidas no deberían dar None");
        assert_eq!(decode_kmer(code, 4), "ATGC");
    }

    #[test]
    fn encode_rejects_invalid_base() {
        assert_eq!(encode_kmer(b"ATNC"), None);
        assert_eq!(encode_kmer(b"atgc"), None); // minúsculas también se rechazan
    }

    #[test]
    fn count_kmers_known_sequence() {
        // "ATGCATGC", k=3 -> ventanas: ATG, TGC, GCA, CAT, ATG, TGC
        // ATG aparece 2 veces, TGC aparece 2 veces, GCA y CAT 1 vez cada uno.
        let seq = b"ATGCATGC";
        let counts = count_kmers(seq, 3);

        let atg = encode_kmer(b"ATG").unwrap();
        let tgc = encode_kmer(b"TGC").unwrap();
        let gca = encode_kmer(b"GCA").unwrap();
        let cat = encode_kmer(b"CAT").unwrap();

        assert_eq!(counts.get(&atg), Some(&2));
        assert_eq!(counts.get(&tgc), Some(&2));
        assert_eq!(counts.get(&gca), Some(&1));
        assert_eq!(counts.get(&cat), Some(&1));
        assert_eq!(counts.len(), 4);
    }

    #[test]
    fn gc_content_known_values() {
        // ATG: A(00) y T(11) no cuentan, G(10) sí -> 1 de 3 -> ~0.333
        let atg = encode_kmer(b"ATG").unwrap();
        assert!((gc_content(atg, 3) - (1.0 / 3.0)).abs() < 1e-6);

        // GCGCGCG: las 7 bases son G o C -> 100% GC
        let gcgcgcg = encode_kmer(b"GCGCGCG").unwrap();
        assert!((gc_content(gcgcgcg, 7) - 1.0).abs() < 1e-6);

        // AAAA: ninguna es G/C -> 0% GC
        let aaaa = encode_kmer(b"AAAA").unwrap();
        assert_eq!(gc_content(aaaa, 4), 0.0);
    }

    #[test]
    fn shannon_entropy_extremes() {
        // Homopolímero puro: una sola categoría con toda la probabilidad
        // -> entropía 0 (mínima diversidad posible).
        let aaaa = encode_kmer(b"AAAA").unwrap();
        assert_eq!(shannon_entropy(aaaa, 4), 0.0);

        // Las 4 bases en proporciones exactamente iguales -> entropía
        // máxima = log2(4) = 2.0.
        let acgt = encode_kmer(b"ACGT").unwrap();
        assert!((shannon_entropy(acgt, 4) - 2.0).abs() < 1e-6);
    }

    #[test]
    fn shannon_entropy_between_extremes() {
        // "AATT": dos categorías (A y T) en proporciones iguales (50/50)
        // -> entropía = 1.0 (log2(2), no log2(4), porque solo 2 categorías
        // tienen probabilidad distinta de cero).
        let aatt = encode_kmer(b"AATT").unwrap();
        assert!((shannon_entropy(aatt, 4) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn reverse_complement_known_value() {
        // ATG -> complemento base a base: T,A,C -> invertir orden: C,A,T
        // Reverso complementario de "ATG" es "CAT".
        let atg = encode_kmer(b"ATG").unwrap();
        let cat = encode_kmer(b"CAT").unwrap();
        assert_eq!(reverse_complement(atg, 3), cat);
    }

    #[test]
    fn reverse_complement_is_involution() {
        // Aplicar reverse_complement dos veces regresa al valor original,
        // porque tanto "complementar" como "invertir orden" son operaciones
        // que se deshacen a sí mismas al repetirlas.
        let kmer = encode_kmer(b"GATTACA").unwrap();
        let rc = reverse_complement(kmer, 7);
        let rc_rc = reverse_complement(rc, 7);
        assert_eq!(kmer, rc_rc);
    }

    #[test]
    fn get_canonical_is_idempotent() {
        let kmer = encode_kmer(b"GATTACA").unwrap();
        let canon = get_canonical(kmer, 7);
        assert_eq!(get_canonical(canon, 7), canon);
    }

    #[test]
    fn get_canonical_matches_for_kmer_and_its_reverse_complement() {
        // Un k-mero y su reverso complementario deben mapear a LA MISMA
        // forma canónica -- ese es justamente el propósito de canonicalizar.
        let atg = encode_kmer(b"ATG").unwrap();
        let cat = encode_kmer(b"CAT").unwrap();
        assert_eq!(get_canonical(atg, 3), get_canonical(cat, 3));
    }

    #[test]
    fn count_kmers_excludes_windows_with_n() {
        // "ATGNATGC", k=3
        // Ventanas: ATG, TGN, GNA, NAT, ATG, TGC
        // Las que contienen la N (TGN, GNA, NAT) se descartan.
        // Sobreviven: ATG (x2), TGC (x1)
        let seq = b"ATGNATGC";
        let counts = count_kmers(seq, 3);

        let atg = encode_kmer(b"ATG").unwrap();
        let tgc = encode_kmer(b"TGC").unwrap();

        assert_eq!(counts.get(&atg), Some(&2));
        assert_eq!(counts.get(&tgc), Some(&1));
        assert_eq!(counts.len(), 2);
    }
}