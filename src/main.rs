use anyhow::Error;
use std::collections::HashSet;

/// Main Function
fn main() {
    //Create a reader for a reference genome
    let reader_genome = match rust_htslib::faidx::Reader::from_path("inst/genome1.fasta") {
        Ok(reader) => reader,
        Err(e) => {
            eprintln!("Failed to read the genome file {}", e);
            std::process::exit(1)
        }
    };

    let valid_chroms = fetch_all_chromosome_names_from_genome(&reader_genome);
    let sequence_context =
        fetch_sequence_from_genome(&reader_genome, "chr1", 100000, 100000, &valid_chroms, true);

    println!("Sequence Context: {:?}", sequence_context);
}

/// A simple function to retrieve a sequence based on 0-based start and end (both-end inclusive)
fn fetch_sequence_from_genome(
    reader: &rust_htslib::faidx::Reader,
    chrom: &str,
    begin: usize,
    end: usize,
    valid_chroms: &HashSet<String>,
    uppercase: bool,
) -> anyhow::Result<String> {
    if !valid_chroms.contains(chrom) {
        return Err(Error::msg("Chrom is not valid!"));
    }

    reader
        .fetch_seq_string(chrom, begin, end)
        .map(|seq| if uppercase { seq.to_uppercase() } else { seq })
        .map_err(anyhow::Error::new)
}

fn fetch_all_chromosome_names_from_genome(reader: &rust_htslib::faidx::Reader) -> HashSet<String> {
    (0..reader.n_seqs())
        .map(|i| reader.seq_name(i as i32))
        .collect::<Result<HashSet<String>, _>>()
        .unwrap_or_else(|e| {
            eprintln!("Failed to create our list of valid chromosome names from the reference genome fasta: {:?}", e);
            std::process::exit(1)
        })
}
