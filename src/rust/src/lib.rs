#![recursion_limit = "256"]
use extendr_api::prelude::*;

mod bas;
mod grts_h3;
mod grts_sobol;

// ---------------------------------------------------------------------------
// extendr module registration
// ---------------------------------------------------------------------------

extendr_module! {
    mod h3sampling;
    use bas;
    use grts_h3;
    use grts_sobol;
}
