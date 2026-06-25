use std::fs;
use std::io::{self, BufWriter, Read};
use std::path::PathBuf;

use clap::{Parser, ValueEnum};
use rsomics_common::{CommonFlags, Result, ToolMeta};

use crate::contacts::{ContactSet, Level, contacts};
use crate::{output, pdb};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum LevelArg {
    Residue,
    Atom,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum MetricArg {
    Euclidean,
}

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-pdb-contacts",
    version,
    about = "Residue/atom contact pairs within a distance cutoff from a PDB (Bio.PDB.NeighborSearch-compatible)"
)]
pub struct Cli {
    /// Input PDB file, or - for stdin.
    pub pdb: PathBuf,

    /// Contact distance cutoff in Angstrom (inclusive: distance <= cutoff).
    #[arg(long, default_value_t = 8.0)]
    pub cutoff: f64,

    /// Emit residue pairs or atom pairs.
    #[arg(long, value_enum, default_value = "residue")]
    pub level: LevelArg,

    /// Distance metric (only euclidean is defined).
    #[arg(long, value_enum, default_value = "euclidean")]
    pub metric: MetricArg,

    /// Restrict the atom universe to Cα atoms (Cα-Cα contact convention).
    #[arg(long = "ca-only")]
    pub ca_only: bool,

    /// Skip intra-chain pairs with |resseq_i - resseq_j| below this value
    /// (0 = no sequence-separation filter).
    #[arg(long = "min-seq-sep", default_value_t = 0)]
    pub min_seq_sep: i32,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    fn level(&self) -> Level {
        match self.level {
            LevelArg::Residue => Level::Residue,
            LevelArg::Atom => Level::Atom,
        }
    }

    pub fn compute(self) -> Result<ContactSet> {
        let text = if self.pdb.as_os_str() == "-" {
            let mut s = String::new();
            io::stdin().read_to_string(&mut s)?;
            s
        } else {
            fs::read_to_string(&self.pdb)?
        };
        let atoms = pdb::parse(&text)?;
        let set = contacts(
            &atoms,
            self.cutoff,
            self.level(),
            self.ca_only,
            self.min_seq_sep,
        );

        if !self.common.json {
            let stdout = io::stdout();
            let mut w = BufWriter::new(stdout.lock());
            output::write_table(&mut w, &set)?;
        }
        Ok(set)
    }
}

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
