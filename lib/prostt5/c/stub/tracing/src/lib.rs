pub mod tracing {
    // Mock representation of the Level enum used for the span!
    #[derive(Debug, Clone)]
    pub enum Level {
        TRACE,
        DEBUG,
        INFO,
        WARN,
        ERROR,
    }

    // Mock representation of a Span, which does nothing
    #[derive(Debug, Clone)]
    pub struct Span;

    impl Span {
        pub fn new() -> Self {
            Span
        }

        pub fn enter(&self) {

        }
    }

    // Define a macro that mimics the tracing::span! macro
    #[macro_export]
    macro_rules! span {
        ($level:expr, $name:expr) => {
            {
                // println!("Span created at level {:?} with name {}", $level, $name);
                $crate::tracing::Span::new()
            }
        };
    }

    // Make the macro accessible as part of the `tracing` module
    pub use crate::span;
}

pub use tracing::{Level, Span};