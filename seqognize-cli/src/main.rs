#![allow(dead_code)]

use clap::{App, Arg, ArgMatches};
use seqognize::aligner::Aligner;
use seqognize::nt_aligner::{GlobalNtAligner, NtAlignmentConfig};
use std::fmt::Debug;
use std::str::FromStr;

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

    let reference = matches.value_of("reference").expect("reference is required").as_bytes();
    let subject = matches.value_of("subject").expect("subject is required").as_bytes();

    let mut aligner = GlobalNtAligner::<_>::new(
        NtAlignmentConfig::new(
            arg(&matches, "match", 1i16),
            arg(&matches, "mismatch", -1i16),
            arg(&matches, "subject_gap", -1i16),
            arg(&matches, "reference_gap", -1i16),
        ),
        reference.to_vec()
    ).expect("Failed to create aligner");

    let alignment = aligner.align(&subject).expect("Alignment failed");
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