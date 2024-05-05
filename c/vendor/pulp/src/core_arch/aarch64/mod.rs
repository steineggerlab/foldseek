macro_rules! __impl {
    ($name: ident, $feature: tt) => {
        #[derive(Clone, Copy)]
        #[repr(transparent)]
        pub struct $name {
            __private: (),
        }

        impl $name {
            #[inline(always)]
            pub fn new_unchecked() -> Self {
                Self { __private: () }
            }
            #[inline(always)]
            pub fn try_new() -> Option<Self> {
                if feature_detected!($feature) {
                    Some(Self { __private: () })
                } else {
                    None
                }
            }

            #[inline(always)]
            pub fn is_available() -> bool {
                feature_detected!($feature)
            }
        }

        impl ::core::fmt::Debug for $name {
            #[inline]
            fn fmt(&self, f: &mut ::core::fmt::Formatter<'_>) -> core::fmt::Result {
                f.write_str(stringify!($name))
            }
        }
    };
}

__impl!(Neon, "neon");
__impl!(Pmull, "pmull");
__impl!(Fp, "fp");
__impl!(Fp16, "fp16");
__impl!(Sve, "sve");
__impl!(Crc, "crc");
__impl!(Lse, "lse");
__impl!(Lse2, "lse2");
__impl!(Rdm, "rdm");
__impl!(Rcpc, "rcpc");
__impl!(Rcpc2, "rcpc2");
__impl!(Dotprod, "dotprod");
__impl!(Tme, "tme");
__impl!(Fhm, "fhm");
__impl!(Dit, "dit");
__impl!(Flagm, "flagm");
__impl!(Ssbs, "ssbs");
__impl!(Sb, "sb");
__impl!(Paca, "paca");
__impl!(Pacg, "pacg");
__impl!(Dpb, "dpb");
__impl!(Dpb2, "dpb2");
__impl!(Sve2, "sve2");
__impl!(Sve2Aes, "sve2-aes");
__impl!(Sve2Sm4, "sve2-sm4");
__impl!(Sve2Sha3, "sve2-sha3");
__impl!(Sve2Bitperm, "sve2-bitperm");
__impl!(Frintts, "frintts");
__impl!(I8mm, "i8mm");
__impl!(F32mm, "f32mm");
__impl!(F64mm, "f64mm");
__impl!(Bf16, "bf16");
__impl!(Rand, "rand");
__impl!(Bti, "bti");
__impl!(Mte, "mte");
__impl!(Jsconv, "jsconv");
__impl!(Fcma, "fcma");
__impl!(Aes, "aes");
__impl!(Sha2, "sha2");
__impl!(Sha3, "sha3");
__impl!(Sm4, "sm4");
__impl!(Asimd, "asimd");
