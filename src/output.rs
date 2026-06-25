use std::io::Write;

use rsomics_common::Result;
use rsomics_common::error::RsomicsError;

use crate::contacts::ContactSet;

fn restag(resseq: i32, icode: char) -> String {
    if icode == ' ' {
        resseq.to_string()
    } else {
        format!("{resseq}{icode}")
    }
}

pub fn write_table(out: &mut dyn Write, set: &ContactSet) -> Result<()> {
    let mut buf = String::new();
    match set {
        ContactSet::Residue(pairs) => {
            buf.push_str("chainA\tresA\tresnameA\tchainB\tresB\tresnameB\n");
            for p in pairs {
                buf.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\n",
                    p.chain_a,
                    restag(p.resseq_a, p.icode_a),
                    p.resname_a,
                    p.chain_b,
                    restag(p.resseq_b, p.icode_b),
                    p.resname_b,
                ));
            }
        }
        ContactSet::Atom(pairs) => {
            buf.push_str("chainA\tresA\tresnameA\tatomA\tchainB\tresB\tresnameB\tatomB\n");
            for p in pairs {
                buf.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    p.chain_a,
                    restag(p.resseq_a, p.icode_a),
                    p.resname_a,
                    p.atom_a,
                    p.chain_b,
                    restag(p.resseq_b, p.icode_b),
                    p.resname_b,
                    p.atom_b,
                ));
            }
        }
    }
    out.write_all(buf.as_bytes()).map_err(RsomicsError::Io)?;
    Ok(())
}
