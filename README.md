# fishgill
filter BAM file and convert it to RAD

## filter
filter function filter the alignments that are outside the corresponding txp terminal kilobase. 
According to the biological fact of drop seq, all reads are generated near $$3'$$ end, so we filter out the unrealistic alignments to increase our confidence.

## convert
codes are from [alevin-fry](https://github.com/COMBINE-lab/alevin-fry#:~:text=Fork%200-,%F0%9F%90%9F%20%F0%9F%94%AC%20%F0%9F%A6%80%20alevin%2Dfry%20is%20an%20efficient%20and%20flexible,cell%20transcriptomics%20and%20feature%20barcoding.)

This function converts a BAM file to a RAD file.

To use fishgill, just run

```
git clone https://github.com/DongzeHE/fishgill.git
cd fishgill
cargo build --release
```

To run filter function, an input BAM file `-b input_bam.bam` and a `-l txplen.tsv` file, in which each row is the transcript name and its txp length are needed. `-t` specifies how many threads will be used in the filtering process. This function will return a filtered bam file in the input bam file's folter, named `input_bam_filtered.bam`

```
target/release/fishgill -b <bam file> -l <txplen.tsv> -t 12
```

To run convert code, input bam file `-b input_bam_filtered.bam`, output RAD file `-o map.rad` are needed, `-t` specifies how many threads will be used in the filtering process.

```
target/release/fishgill -b <bam file> -o <map.rad> -t 12
```
