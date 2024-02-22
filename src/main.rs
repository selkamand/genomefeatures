use anyhow::Error;
use std::collections::HashSet;
use std::collections::HashMap;
use std::fmt;

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

    let chrom_lengths = fetch_all_chromosome_names_and_lengths(&reader_genome);
    println!("{:?}", chrom_lengths);

    let sequence =
        fetch_sequence_from_genome(&reader_genome, "chr1", 11, 11, &chrom_lengths, true);

    let sequence_context = fetch_flanking_sequence(&reader_genome, "chr1", 1, 2, 2, &chrom_lengths, true);
    println!("Sequence: {:?}", sequence);
    println!("Sequence Context: {:?}", sequence_context)
}

/// Retrieve a sequence based on 0-based start and end (both-end inclusive)
fn fetch_sequence_from_genome(
    reader: &rust_htslib::faidx::Reader,
    chrom: &str,
    begin: usize,
    end: usize,
    chrom_lengths: &HashMap<String, u64>,
    uppercase: bool,
) -> anyhow::Result<String> {
    // Check if Seqnames are in the genome
    if !chrom_lengths.contains_key(chrom) {
        return Err(Error::msg("Chrom is not valid!"));
    }

    // Cheeck if full range of sequences are in the genome
    let chrom_length = match chrom_lengths.get(chrom){
        Some(len) => len,
        None => {
            eprintln!("Failed to find a valid chromosome length for chrom {}", chrom);
            std::process::exit(1);
        }
    };
    if *chrom_length <= end as u64 {
        return Err(Error::msg("End position is longer than chromosome"))
    }
    if *chrom_length <= begin as u64 {
        return Err(Error::msg("Begin position is longer than chromosome"))
    }
    // Fetch the sequence
    reader
        .fetch_seq_string(chrom, begin, end)
        .map(|seq| if uppercase { seq.to_uppercase() } else { seq })
        .map_err(anyhow::Error::new)
}

/// Fetch All Chromosome Names From Genome
fn fetch_all_chromosome_names_from_genome(reader: &rust_htslib::faidx::Reader) -> HashSet<String> {
    (0..reader.n_seqs())
        .map(|i| reader.seq_name(i as i32))
        .collect::<Result<HashSet<String>, _>>()
        .unwrap_or_else(|e| {
            eprintln!("Failed to create our list of valid chromosome names from the reference genome fasta: {:?}", e);
            std::process::exit(1)
        })
}



// Assuming `Error` is from rust_htslib::errors::Error, but adjust as necessary




fn fetch_all_chromosome_names_and_lengths(reader: &rust_htslib::faidx::Reader) -> HashMap<String, u64> {
    let chromosomes = fetch_all_chromosome_names_from_genome(reader);
    let lengths: Vec<u64> = chromosomes.clone().into_iter().map(|name| reader.fetch_seq_len(name) ).collect();
    let map: HashMap<String, u64> = chromosomes.into_iter().zip(lengths).collect();


    // Return our map
    map
}


/// Fetch sequence region from position-upstream to postion + downstream
/// Position should be 0-based
fn fetch_flanking_sequence(
    reader: &rust_htslib::faidx::Reader,
    chrom: &str,
    position: usize,
    upstream: usize,
    downstream: usize,
    chrom_lengths: &HashMap<String, u64>,
    uppercase: bool
) -> anyhow::Result<String>{
    // Check upstream is NOT to large
    if upstream >= position{
        return Err(Error::msg("Cannot fetch enough bases upstream since position is too close to the start of sequence"))
    }
    let begin: usize = position - upstream;
    let end: usize = position + downstream;

    fetch_sequence_from_genome(reader, chrom, begin, end, chrom_lengths, uppercase)
}

/// Produce Count Matrix 
/// For each nucleotide type - count the number of occurences at each position in a String
fn count(sequence: Vec<String>, expected_length: usize) -> anyhow::Result<> {
    if sequence.len() != expected_length {
        return Err(Error::msg(format!("Sequence was not the expected length {}", expected_length)))
    }
}

