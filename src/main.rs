
extern crate chrono;
extern crate clap;
extern crate num_cpus;
extern crate slog;
extern crate slog_term;

use clap::{crate_authors, crate_version, App, Arg};
use slog::{o, Drain};
use arms;





fn main() {

    let max_num_threads: String = (num_cpus::get() as u32).to_string();
    let crate_authors = crate_authors!("\n");
    let version = crate_version!();
    // [] add command for just counting barcode frequency
    // [] add other algorithms for determining barcode cutoff

    let convert_app = App::new("convert")
        .about("Convert a BAM file to a RAD file")
        .version(version)
        .author(crate_authors)
        .arg(Arg::from("-b, --bam=<bam-file> 'input SAM/BAM file'"))
        .arg(
            Arg::from("-t, --threads 'number of threads to use for processing'")
                .default_value(&max_num_threads),
        )
        .arg(Arg::from("-o, --output=<rad-file> 'output RAD file'"));


    let filter_app = App::new("filter")
    .about("remove alignments outside terminal kilobase")
    .version(version)
    .author(crate_authors)
    .arg(Arg::from("-b, --bam=<bam-file> 'input SAM/BAM file'"))
    .arg(Arg::from("-l, --txplen=<txplen-file> 'input txplen tsv file'"))
    .arg(
        Arg::from("-t, --threads 'number of threads to use for processing'")
            .default_value(&max_num_threads),
    );

    let opts = App::new("fishgill")
    .version(version)
    .author(crate_authors)
    .about("Cook BAM for alevin-fry from the command line")
    .subcommand(filter_app)
    .subcommand(convert_app)
    .get_matches();

    
    let decorator = slog_term::TermDecorator::new().build();
    let drain = slog_term::CompactFormat::new(decorator)
        .use_custom_timestamp(|out: &mut dyn std::io::Write| {
            write!(out, "{}", chrono::Local::now().format("%Y-%m-%d %H:%M:%S")).unwrap();
            Ok(())
        })
        .build()
        .fuse();
    let drain = slog_async::Async::new(drain).build().fuse();

    let log = slog::Logger::root(drain, o!());

    if let Some(ref t) = opts.subcommand_matches("convert") {
        let input_file: String = t.value_of_t("bam").unwrap();
        let rad_file: String = t.value_of_t("txplen").unwrap();
        let num_threads: u32 = t.value_of_t("threads").unwrap();
        arms::convert::bam2rad(input_file, rad_file, num_threads, &log)
    }

    if let Some(ref t) = opts.subcommand_matches("filter") {
        let bam_file: String = t.value_of_t("bam").unwrap();
        let txplen_file: String = t.value_of_t("txplen").unwrap();
        let num_threads: usize = t.value_of_t("threads").unwrap();
        arms::filter::filter_bam(&bam_file, &txplen_file, num_threads, &log)
    }

}
