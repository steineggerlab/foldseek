#[derive(Clone)]
pub struct Hasher {
}

impl Hasher {
    pub fn new() -> Self {
        return Hasher{}
    }

    pub fn new_with_initial(_init: u32) -> Self {
        Self::new()
    }

    pub fn new_with_initial_len(_init: u32, _amount: u64) -> Self {
        Self::new()
    }

    pub fn update(&mut self, _buf: &[u8]) {
    }

    pub fn finalize(self) -> u32 {
        0
    }

    pub fn reset(&mut self) {
    }

    pub fn combine(&mut self, _other: &Self) {
    }
}

impl Default for Hasher {
    fn default() -> Self {
        Self::new()
    }
}