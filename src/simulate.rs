//! Utility functions for simulating random sequences.

use rand::prelude::*;

/// All 20 amino acids.
pub static AMINO_ACIDS: [u8; 20] = [
    b'A', b'C', b'D', b'E', b'F', b'G', b'H', b'I', b'K', b'L',
    b'M', b'N', b'P', b'Q', b'R', b'S', b'T', b'V', b'W', b'Y'
];

/// All 4 nucleotides and the undetermined/placeholder nucleotide.
pub static NUC: [u8; 5] = [
    b'A', b'C', b'G', b'N', b'T'
];

/// Given an input byte string, create a randomly mutated copy and
/// add random suffixes to both strings.
pub fn rand_mutate_suffix<R: Rng>(a: &mut Vec<u8>, k: usize, alpha: &[u8], suffix_len: usize, rng: &mut R) -> Vec<u8> {
    let mut b = rand_mutate(a, k, alpha, rng);
    let a_suffix = rand_str(suffix_len, alpha, rng);
    let b_suffix = rand_str(suffix_len, alpha, rng);
    a.extend_from_slice(&a_suffix);
    b.extend_from_slice(&b_suffix);
    b
}

/// Given an input byte string, create a randomly mutated copy with
/// a single long random insert.
pub fn rand_mutate_insert<R: Rng>(a: &[u8], k: usize, alpha: &[u8], insert_len: usize, rng: &mut R) -> Vec<u8> {
    let b = rand_mutate(a, k, alpha, rng);
    let insert = rand_str(insert_len, alpha, rng);
    let idx = rng.gen_range(1..b.len());
    let mut res = Vec::with_capacity(b.len() + insert_len);
    // insert the long insert string
    res.extend_from_slice(&b[..idx]);
    res.extend_from_slice(&insert);
    res.extend_from_slice(&b[idx..]);
    res
}

/// Given an input byte string, craete a randomly mutated copy.
pub fn rand_mutate<R: Rng>(a: &[u8], k: usize, alpha: &[u8], rng: &mut R) -> Vec<u8> {
    let mut b = vec![];
    let mut i = 0;
    let mut edits = 0;

    while i < a.len() {
        if edits < k && rng.gen_range(0..a.len()) < k {
            let edit = rng.gen_range(1..4);

            match edit {
                1 => { // mismatch
                    let mut iter = alpha.choose_multiple(rng, 2);
                    let first = *iter.next().unwrap();
                    let second = *iter.next().unwrap();
                    b.push(if first == a[i] { second } else { first });
                    i += 1;
                },
                2 => { // delete
                    i += 1;
                },
                3 => { // insert
                    b.push(*alpha.choose(rng).unwrap());
                },
                _ => panic!("Not possible!")
            }

            edits += 1;
        } else {
            b.push(a[i]);
            i += 1;
        }
    }

    b
}

/// Generate a random string of a certain length, with a certain
/// alphabet.
pub fn rand_str<R: Rng>(length: usize, alpha: &[u8], rng: &mut R) -> Vec<u8> {
    let mut res = vec![0u8; length];

    for i in 0..length {
        res[i] = *alpha.choose(rng).unwrap();
    }

    res
}
