//! SIMD-accelerated library for computing global and X-drop affine
//! gap penalty sequence-to-sequence or sequence-to-profile alignments
//! using an adaptive block-based algorithm.
//!
//! Currently, AVX2, Neon, and WASM SIMD are supported.
//!
//! ## Example
//! ```
//! use block_aligner::scan_block::*;
//! use block_aligner::scores::*;
//! use block_aligner::cigar::*;
//!
//! let block_size = 16;
//! let gaps = Gaps { open: -2, extend: -1 };
//! let r = PaddedBytes::from_bytes::<NucMatrix>(b"TTAAAAAAATTTTTTTTTTTT", block_size);
//! let q = PaddedBytes::from_bytes::<NucMatrix>(b"TTTTTTTTAAAAAAATTTTTTTTT", block_size);
//!
//! // Align with traceback, but no x drop threshold.
//! let mut a = Block::<true, false>::new(q.len(), r.len(), block_size);
//! a.align(&q, &r, &NW1, gaps, block_size..=block_size, 0);
//! let res = a.res();
//!
//! assert_eq!(res, AlignResult { score: 7, query_idx: 24, reference_idx: 21 });
//!
//! let mut cigar = Cigar::new(res.query_idx, res.reference_idx);
//! a.trace().cigar_eq(&q, &r, res.query_idx, res.reference_idx, &mut cigar);
//!
//! assert_eq!(cigar.to_string(), "2=6I16=3D");
//! ```
//!
//! ## Tuning block sizes
//!
//! For long, noisy Nanopore reads, a min block size of ~1% sequence length and a max block size
//! of ~10% sequence length performs well (tested with reads up to ~50kbps).
//! For proteins, a min block size of 32 and a max block size of 256 performs well.
//! Using a minimum block size of 32 is recommended for most applications.
//! Let me know how block aligner performs on your data!
//!
//! When building your code that uses this library, it is important to specify the
//! correct feature flags: `simd_avx2`, `simd_neon`, or `simd_wasm`.
//! More information on specifying different features for different platforms
//! with the same dependency [here](https://doc.rust-lang.org/cargo/reference/specifying-dependencies.html#platform-specific-dependencies).

// special SIMD instruction set modules adapted for this library
// their types and lengths are abstracted out

#[cfg(feature = "simd_avx2")]
#[macro_use]
#[doc(hidden)]
/// cbindgen:ignore
pub mod avx2;

#[cfg(feature = "simd_avx2")]
pub use avx2::L;

#[cfg(feature = "simd_wasm")]
#[macro_use]
#[doc(hidden)]
/// cbindgen:ignore
pub mod simd128;

#[cfg(feature = "simd_wasm")]
pub use simd128::L;

#[cfg(feature = "simd_neon")]
#[macro_use]
#[doc(hidden)]
/// cbindgen:ignore
pub mod neon;

#[cfg(feature = "simd_neon")]
pub use neon::L;

#[cfg(any(feature = "simd_avx2", feature = "simd_wasm", feature = "simd_neon"))]
pub mod scan_block;
#[cfg(any(feature = "simd_avx2", feature = "simd_wasm", feature = "simd_neon"))]
pub mod scores;
#[cfg(any(feature = "simd_avx2", feature = "simd_wasm", feature = "simd_neon"))]
pub mod cigar;

#[cfg(any(feature = "simd_avx2", feature = "simd_neon"))]
#[doc(hidden)]
pub mod ffi;

/// Calculate the percentage of a length, rounded to the next power of two.
///
/// This is useful for computing the min and max block sizes for sequences of a certain
/// length by using percentages. The returned value is at least 32.
pub fn percent_len(len: usize, p: f32) -> usize {
    ((p * (len as f32)).round() as usize).max(32).next_power_of_two()
}
