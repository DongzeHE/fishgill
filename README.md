# fishgill
filter BAM file and convert it to RAD

## filter
filter function filter the alignments that are outside the corresponding txp terminal kilobase. 
According to the biological fact of drop seq, all reads are generated near $$3'$$ end, so we filter out the unrealistic alignments to increase our confidence.

## convert
codes are from [alevin-fry](https://github.com/COMBINE-lab/alevin-fry#:~:text=Fork%200-,%F0%9F%90%9F%20%F0%9F%94%AC%20%F0%9F%A6%80%20alevin%2Dfry%20is%20an%20efficient%20and%20flexible,cell%20transcriptomics%20and%20feature%20barcoding.)

This function converts a BAM file to a RAD file.
