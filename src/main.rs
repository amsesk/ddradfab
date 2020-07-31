use bio::io::fasta::Reader;
use clap::{App, Arg};
use std::io::Error;

fn main() -> Result<(), Error> {
    let args = App::new("ddradfab")
        .version("0.1")
        .author("Kevin Amses")
        .about("Generates fabricated EcoRI-MseI ddRAD tags from a nucleotide FASTA.")
        .arg(
            Arg::with_name("fasta")
                .short('f')
                .long("fasta")
                .value_name("FASTA")
                .about("Path to nucelteotide fasta.")
                .takes_value(true)
                .required(true),
        )
        .get_matches();

    let ecori = b"GAATTC";
    let msei = b"TTAA";
    let MAXSCANDIST = 600;

    let fasta = Reader::from_file(args.value_of("fasta").unwrap())?;
    let mut fa_head_idx = 1;
    for record in fasta.records() {
        let sequence = record.unwrap().seq().to_owned();
        for p in 0..sequence.len() - 6 {
            if &sequence[p..p + 6] == ecori {
                let scandist = |p| match p + MAXSCANDIST <= sequence.len() {
                    true => return p + MAXSCANDIST,
                    false => return sequence.len(),
                };
                let potential_tag = &sequence[p..scandist(p)];
                for i in 0..potential_tag.len() - 4 {
                    //Print from 1..i instead of 0..i+1 in order to cut off leading G (EcoRI cut) and trailing T (MseI cut)
                    //Actually print from 1.151 because of 150bp read length
                    if &potential_tag[i..i + 4] == msei {
                        if i <= 450 && i >= 350 {
                            println!(
                                "@fabricated_rad_tag_{}\n{}\n+\n{}",
                                fa_head_idx,
                                String::from_utf8_lossy(&potential_tag[1..141]).to_uppercase(),
                                (0..140).map(|_| "?").collect::<String>()
                            );
                            fa_head_idx += 1;
                        }
                    }
                }
            }
        }
    }
    Ok(())
}
