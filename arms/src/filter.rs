use std::collections::{HashMap, HashSet};
use std::str::from_utf8;
use std::io::stdout;
use std::io::Write;

use itertools::Itertools;

use crate::config::MIL;

use rust_htslib::bam;
use rust_htslib::bam::{Read, Record};
use rust_htslib::bam::record::Aux;
// use rust_htslib::bam::record::{Cigar};
use slog::{info, warn};

#[derive(PartialEq, Eq, Hash)]
enum KeyT {
    HI(i64),
    FEATURE((i32, i64, i64, bool)),
}


pub fn filter_bam(in_bam_file: &str, 
                  out_bam_file: &str, 
                  txplen_file: &str, 
                  num_threads: usize,
                  fltr_unsplcd:bool,
                  log: &slog::Logger,
                ) {
    info!(log, "Using Input BAM file: {}", in_bam_file);

    // opening input BAM
    let mut input_bam = bam::Reader::from_path(in_bam_file)
        .expect("Can't open BAM file");
    let header = bam::Header::from_template(input_bam.header());

    // creating the output BAM
    // let mut out_bam_file = in_bam_file.to_string();
    // let bam_name_offset = out_bam_file.find(".bam")
    //     .unwrap_or(out_bam_file.len());
    // out_bam_file.replace_range(bam_name_offset..,
    //                        "_filtered.bam");
    let mut output_bam = bam::Writer::from_path(out_bam_file, &header, bam::Format::BAM)
        .expect("can't open BAM file to dump output");

    input_bam.set_threads(2)
        .unwrap();
    output_bam.set_threads(num_threads-3)
        .unwrap();

    // Read in ref names and connect it to tid
    let hdrv = input_bam.header().to_owned();

    // I need rname_to_id to process the txp names in txp_len file
    let mut rname_to_id:HashMap<String, u32> = HashMap::with_capacity(hdrv.inner().n_targets as usize);
    // tid to type is used to write some statistics
    let mut tid_to_type = vec![0u8;hdrv.inner().n_targets as usize];

    // construct the two vectors 
    for (tid, tname) in hdrv.target_names().iter().enumerate() {
        // get ref txp name
        let rname = std::str::from_utf8(tname).unwrap().to_string();
        // check the type of txp and write to the tid_to_type
        if rname.find("-U").is_none() {
            tid_to_type[tid] = 1;
        } else {
            tid_to_type[tid] = 2;
        }
        // insert to rname_to_id
        rname_to_id.insert(rname, tid as u32);
    }
    
    let txplen = std::fs::File::open(txplen_file).expect("couldn't open file");

    // we read in txp lenth file 
    // it has txp length in the third field
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(txplen);

    // we build txp length vector
    let mut txplen_vec:Vec<i64> = vec![i64::MAX;rname_to_id.len() as usize];

    // we build txp end exon length vectror
    let mut tx_end_exon_len_vec:Vec<i64> = vec![i64::MAX;rname_to_id.len() as usize];

    // There are three fields in each row: 
    // 1 txp name 
    // 2 end exon length
    // 3 txp length

    type TSVRec = (String, i64,i64);

    // read in records in tsv file
    for result in rdr.deserialize() {
        let record: TSVRec = result.unwrap();
        // if the txp name in the record is in the ref name list, we write down the two lengths
        if let Some(tid) = rname_to_id.get(&record.0){
            tx_end_exon_len_vec[*tid as usize] = record.1;
            // As the pos is 0-based in bam file
            txplen_vec[*tid as usize] = record.2 - 1; 
        }
    }

    // now we start to go through the records grouped by read name
    info!(log, "Starting to filter");
    let mut read_count = 0;
    for (qname, read_group) in input_bam
        .records()
        .map(|res| res.unwrap())
        .group_by(|rec| from_utf8(rec.qname()).unwrap().to_string())
        .into_iter()
    {
        read_count += 1;
        if read_count % MIL == 0 {
            println!("\rDone processing {}M reads", read_count / MIL);
            stdout().flush().expect("Can't flush output");
        }

        let mut alignments: Vec<Record> = read_group.collect();
        let num_alignments = alignments.len();
        let num_reads: i64 = num_alignments as i64 / 2;
        let is_multi_aligned = num_reads > 1;

        let mut is_discordant: bool = false;
        for i in (0..num_alignments).step(2){
            if alignments[i].pos() == alignments[i+1].mpos() &&
                alignments[i+1].pos() == alignments[i].mpos() &&
                alignments[i].pos() != alignments[i].mpos() &&
                alignments[i+1].pos() != alignments[i+1].mpos()
            { continue; }

            is_discordant = true;
            break;
        }//end-for

        let has_hi_tag = match alignments.first()
            .unwrap()
            .aux("HI".as_bytes()) {
                Some(_) => true,
                None => false,
            };

        let has_nh_tag = match alignments.first()
            .unwrap()
            .aux("NH".as_bytes()) {
                Some(_) => true,
                None => false,
            };

        if is_discordant {
            let mut bucket = HashMap::<KeyT, Vec<usize>>::new();
            for (index,alignment) in alignments.iter().enumerate() {
                let key: KeyT;
                if has_hi_tag {
                    key = KeyT::HI(alignment.aux("HI".as_bytes())
                                   .expect("Some alignment doesn't have HI tag")
                    .integer());
                } else {
                    let m1 = std::cmp::min(alignment.pos(), alignment.mpos());
                    let m2 = std::cmp::max(alignment.pos(), alignment.mpos());
                    key = KeyT::FEATURE( (alignment.tid(), m1, m2, alignment.is_secondary()));
                }

                // insert the alignments into bucket
                bucket.entry(key)
                    .or_insert(Vec::new())
                    .push(index);
            }//end-bucket for
            
            let mut skip_read = false;
            let mut skip_alignments = HashSet::new();
            for (key, vals) in &bucket {
                if vals.len() !=2 {
                    warn!(log, "Wrong PE mapping for {}; Skipping it", qname);
                    skip_read = true;
                }

                //iterate over  alignments
                for val in vals {
                    // if the alignment is outside of the terminal kilobase, we skip it
                    // first get the record from alignments
                    let record = &alignments[*val];
                    // compute the distance of pos to 3 prime end
                    let read_start_pos = txplen_vec[record.tid() as usize] - record.pos();
                    // if the distance is over 1000, we skip it
                    if read_start_pos > 1000 && is_multi_aligned {
                        if  fltr_unsplcd{
                            skip_alignments.insert(key);
                        } else if tid_to_type[record.tid() as usize] == 1 {
                            skip_alignments.insert(key);
                        }
                    }

                    // let mut read_length = 0;
                    // for (_segment_index, cigar) in alignments[*val].cigar()
                    //     .iter()
                    //     .enumerate()
                    // {
                    //     match cigar {
                    //         &Cigar::Ins(_) => {
                    //             if read_length < 23 {
                    //                 skip_alignments.insert(key);
                    //                 break;
                    //             }
                    //         },

                    //         &Cigar::Match(l) | &Cigar::Equal(l) |
                    //         &Cigar::RefSkip(l) | &Cigar::Del(l) |
                    //         &Cigar::Pad(l) | &Cigar::SoftClip(l) |
                    //         &Cigar::HardClip(l) | &Cigar::Diff(l) => {
                    //             read_length += l;
                    //             continue
                    //         },
                    //     }//end-match
                    // }//end for
                }//end-for
            }//end-for

            // Skip a read if not enough information to filter
            if skip_read { continue; }
            for (key, vals) in &bucket {
                if skip_alignments.contains(&key) { continue; }

                for id in vals {
                    let alignment = &mut alignments[*id];

                    if has_nh_tag {
                        alignment.remove_aux("NH".as_bytes());
                    }

                    alignment.push_aux("NH".as_bytes(),
                                       &Aux::Integer(num_reads-skip_alignments.len() as i64));
                        // .expect("can't add NH tag");

                    output_bam.write(&alignment)
                        .expect("can't write bam record");
                }
            }
        } else {
            let mut buffer = Vec::new();
            let mut filtered_algns = Vec::new();

            for (index, mut alignment) in alignments.into_iter().enumerate() {
                if has_nh_tag {
                    alignment.remove_aux("NH".as_bytes());
                }

                let mut skip_alignment = false;
                // let mut read_length = 0;

                let read_start_pos = txplen_vec[alignment.tid() as usize] - alignment.pos();
                // if the distance is over 1000, we skip it
                if read_start_pos > 1000 && is_multi_aligned{
                    if  fltr_unsplcd{
                        skip_alignment = true;
                    } else if tid_to_type[alignment.tid() as usize] == 1 {
                        skip_alignment = true;
                    }
                }

                // for (_segment_index, cigar) in alignment.cigar()
                //     .iter()
                //     .enumerate()
                // {
                //     match cigar {
                //         &Cigar::Ins(_) => {
                //             if read_length < 23 {
                //                 skip_alignment = true;
                //                 break;
                //             }
                //         },

                //         &Cigar::Match(l) | &Cigar::Equal(l) |
                //         &Cigar::RefSkip(l) | &Cigar::Del(l) |
                //         &Cigar::Pad(l) | &Cigar::SoftClip(l) |
                //         &Cigar::HardClip(l) | &Cigar::Diff(l) => {
                //             read_length += l;
                //             continue
                //         },
                //     }//end-match
                // }//end for

                if index%2 == 0 {
                    if !skip_alignment { buffer.push(alignment); }
                } else {
                    if !skip_alignment && buffer.len() > 0 {
                        filtered_algns.push(buffer.iter().nth(0).unwrap().clone());
                        filtered_algns.push(alignment);

                        buffer.pop();
                        assert!(buffer.len() == 0, "parsing error");
                    }
                }
            }

            let nh_tag: i64 = filtered_algns.len() as i64 / 2;
            for mut alignment in filtered_algns {
                alignment.push_aux("NH".as_bytes(),
                                   &Aux::Integer(nh_tag));
                    // .expect("can't add NH tag");

                output_bam.write(&alignment)
                    .expect("can't write bam record");
            }
        }//end-else
    }//iterate input bam
}