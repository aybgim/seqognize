use seqognize::nt_aligner::{GlobalNtAligner, NtAlignmentConfig};
use seqognize::aligner::Aligner;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Reference sequence
    #[arg(short = 'r', long = "ref")]
    reference: String,

    /// Subject sequence
    #[arg(short = 's', long = "sub")]
    subject: String,

    /// Match score
    #[arg(short = 'm', long = "match", default_value_t = 1)]
    match_score: i32,

    /// Mismatch penalty
    #[arg(short = 'x', long = "mismatch", default_value_t = -1)]
    mismatch_penalty: i32,

    /// Subject gap opening
    #[arg(long = "sg", default_value_t = -1)]
    subject_gap: i32,

    /// Reference gap opening
    #[arg(long = "rg", default_value_t = -1)]
    reference_gap: i32,

    /// Vertical output
    #[arg(long)]
    vertical: bool,
}

fn main() {
    let args = Args::parse();

    let aligner: GlobalNtAligner = GlobalNtAligner {
        config: NtAlignmentConfig {
            match_score: args.match_score,
            mismatch_penalty: args.mismatch_penalty,
            subject_gap_penalty: args.subject_gap,
            reference_gap_penalty: args.reference_gap,
        }
    };

    let alignment = aligner.align(args.subject.as_bytes(), args.reference.as_bytes());
    println!("Score: {:?}", alignment.score);
    if args.vertical {
        alignment.print_vertical();
    } else {
        alignment.print_horizontal();
    }
}