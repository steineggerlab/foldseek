pub struct Regex;
use std::iter::{self, Once};
use std::io::{self, ErrorKind};

impl Regex {
    pub fn new(_: &str) -> Result<Self, std::io::Error> {
        Ok(Regex)
    }
    pub fn find_iter<'a>(&self, _text: &'a str) -> Once<Result<String, io::Error>> {
        iter::once(Err(io::Error::new(ErrorKind::Other, "no matches found")))
    }
}