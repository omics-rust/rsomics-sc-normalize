use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_sc_normalize::{NormalizeParams, open_output, parse_target_sum, run};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-sc-normalize", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// 10x MTX directory (matrix.mtx[.gz], genes×cells).
    pub input: PathBuf,

    #[arg(short = 'o', long, default_value = "-")]
    output: String,

    /// Per-cell target count, or `median` for scanpy's default.
    #[arg(long = "target-sum", default_value = "median")]
    target_sum: String,

    /// Skip the log1p step, emitting only the size-normalized matrix.
    #[arg(long = "no-log", default_value_t = false)]
    no_log: bool,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        self.common.install_rayon_pool()?;
        let params = NormalizeParams {
            target_sum: parse_target_sum(&self.target_sum)?,
            log1p: !self.no_log,
        };
        let out = open_output(&self.output)?;
        let (genes, cells) = run(&self.input, &params, out)?;
        if !self.common.quiet {
            eprintln!("normalized {cells} cells × {genes} genes");
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Library-size normalization + log1p of a single-cell count matrix.",
    origin: Some(Origin {
        upstream: "scanpy sc.pp.normalize_total + sc.pp.log1p",
        upstream_license: "BSD-3-Clause",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1186/s13059-017-1382-0"),
    }),
    usage_lines: &["<10x-mtx-dir> [--target-sum median|<float>] [--no-log] [-o out.mtx]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("String"),
                required: false,
                default: Some("-"),
                description: "Output MTX path (genes×cells, real values); '-' for stdout.",
                why_default: Some("Streams to stdout for pipeline composition."),
            },
            FlagSpec {
                short: None,
                long: "target-sum",
                aliases: &[],
                value: Some("<median|float>"),
                type_hint: Some("String"),
                required: false,
                default: Some("median"),
                description: "Per-cell target count; 'median' uses the median of cell totals.",
                why_default: Some("Matches scanpy's target_sum=None default."),
            },
            FlagSpec {
                short: None,
                long: "no-log",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Emit the size-normalized matrix without applying log1p.",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "scanpy-default normalization (median target, log1p)",
            command: "rsomics-sc-normalize filtered_feature_bc_matrix/ -o norm.mtx",
        },
        Example {
            description: "CPM normalization (target 1e6), no log",
            command: "rsomics-sc-normalize mtx_dir/ --target-sum 1e6 --no-log -o cpm.mtx",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
