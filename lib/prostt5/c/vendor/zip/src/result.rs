//! Error types that can be emitted from this library

use std::error::Error;
use std::fmt;
use std::io;
use std::io::IntoInnerError;
use std::num::TryFromIntError;

/// Generic result type with ZipError as its error variant
pub type ZipResult<T> = Result<T, ZipError>;

/// Error type for Zip
#[derive(Debug)]
#[non_exhaustive]
pub enum ZipError {
    /// An Error caused by I/O
    Io(io::Error),

    /// This file is probably not a zip archive
    InvalidArchive(&'static str),

    /// This archive is not supported
    UnsupportedArchive(&'static str),

    /// The requested file could not be found in the archive
    FileNotFound,

    /// The password provided is incorrect
    InvalidPassword,
}

impl From<io::Error> for ZipError {
    fn from(err: io::Error) -> ZipError {
        ZipError::Io(err)
    }
}

impl<W> From<IntoInnerError<W>> for ZipError {
    fn from(value: IntoInnerError<W>) -> Self {
        ZipError::Io(value.into_error())
    }
}

impl fmt::Display for ZipError {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ZipError::Io(err) => write!(fmt, "{err}"),
            ZipError::InvalidArchive(err) => write!(fmt, "invalid Zip archive: {err}"),
            ZipError::UnsupportedArchive(err) => write!(fmt, "unsupported Zip archive: {err}"),
            ZipError::FileNotFound => write!(fmt, "specified file not found in archive"),
            ZipError::InvalidPassword => write!(fmt, "incorrect password for encrypted file"),
        }
    }
}

impl Error for ZipError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            ZipError::Io(err) => Some(err),
            _ => None,
        }
    }
}

impl ZipError {
    /// The text used as an error when a password is required and not supplied
    ///
    /// ```rust,no_run
    /// # use zip::result::ZipError;
    /// # let mut archive = zip::ZipArchive::new(std::io::Cursor::new(&[])).unwrap();
    /// match archive.by_index(1) {
    ///     Err(ZipError::UnsupportedArchive(ZipError::PASSWORD_REQUIRED)) => eprintln!("a password is needed to unzip this file"),
    ///     _ => (),
    /// }
    /// # ()
    /// ```
    pub const PASSWORD_REQUIRED: &'static str = "Password required to decrypt file";
}

impl From<ZipError> for io::Error {
    fn from(err: ZipError) -> io::Error {
        let kind = match &err {
            ZipError::Io(err) => err.kind(),
            ZipError::InvalidArchive(_) => io::ErrorKind::InvalidData,
            ZipError::UnsupportedArchive(_) => io::ErrorKind::Unsupported,
            ZipError::FileNotFound => io::ErrorKind::NotFound,
            ZipError::InvalidPassword => io::ErrorKind::InvalidInput,
        };

        io::Error::new(kind, err)
    }
}

/// Error type for time parsing
#[derive(Debug)]
pub struct DateTimeRangeError;

// TryFromIntError is also an out-of-range error.
impl From<TryFromIntError> for DateTimeRangeError {
    fn from(_value: TryFromIntError) -> Self {
        DateTimeRangeError
    }
}

impl fmt::Display for DateTimeRangeError {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        write!(
            fmt,
            "a date could not be represented within the bounds the MS-DOS date range (1980-2107)"
        )
    }
}

impl Error for DateTimeRangeError {}
