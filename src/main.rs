use clap::Parser;
use rsomics_common::run;
use rsomics_pdb_contacts::cli::{Cli, META};

fn main() -> std::process::ExitCode {
    let cli = Cli::parse();
    let common = cli.common.clone();
    run(&common, META, || cli.compute())
}
