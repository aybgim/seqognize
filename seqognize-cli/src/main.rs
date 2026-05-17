#![allow(dead_code)]

use seqognize::nt_aligner::{GlobalNtAligner, NtAlignmentConfig};
use seqognize::aligner::Aligner;
use clap::{App, Arg, ArgMatches};
use std::str::FromStr;
use std::fmt::Debug;

fn main() {
    let matches = App::new("Seqognize")
        .version("1.0")
        .author("Albert Gevorgyan. <ablertus@yahoo.com>")
        .about("Sequence analysis tool.")
        .arg(Arg::with_name("reference")
            .short("r")
            .long("ref")
            .help("Reference sequence")
            .required(true)
            .takes_value(true))
        .arg(Arg::with_name("subject")
            .short("s")
            .long("sub")
            .help("Subject sequence")
            .required(true)
            .takes_value(true))
        .arg(Arg::with_name("match")
            .short("m")
            .long("match")
            .help("Match score")
            .takes_value(true))
        .arg(Arg::with_name("mismatch")
            .short("x")
            .long("mismatch")
            .help("Mismatch penalty")
            .takes_value(true))
        .arg(Arg::with_name("subject_gap")
            // .short("sg")
            .long("sg")
            .help("Subject gap opening")
            .takes_value(true))
        .arg(Arg::with_name("reference_gap")
            // .short("rg")
            .long("rg")
            .help("Reference gap opening")
            .takes_value(true))
        .arg(Arg::with_name("vertical")
            .long("vertical")
            .help("Vertical output")
            .takes_value(false))
        .get_matches();

    let reference = matches.value_of("reference").unwrap().as_bytes();
    let subject = matches.value_of("subject").unwrap().as_bytes();

    let aligner: GlobalNtAligner = GlobalNtAligner {
        config: NtAlignmentConfig {
            match_score: arg(&matches, "match", 1.0),
            mismatch_penalty: arg(&matches, "mismatch", -1.0),
            subject_gap_penalty: arg(&matches, "subject_gap", -1.0),
            reference_gap_penalty: arg(&matches, "reference_gap", -1.0),
        }
    };

    let alignment = aligner.align(&subject, &reference);
    println!("Score: {:?}", alignment.score);
    if matches.is_present("vertical") {
        alignment.print_vertical();
    } else {
        alignment.print_horizontal();
    }
}

fn arg<T: FromStr + Debug>(matches: &ArgMatches, argname: &str, default: T) -> T
    where <T as std::str::FromStr>::Err: std::fmt::Debug {
    matches.value_of(argname).map(|s| s.parse().unwrap()).unwrap_or(default)
}