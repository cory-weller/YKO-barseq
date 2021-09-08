# 1011-barcodes

## Building singularity image
On a computer with sudo access, run:
```bash
# while in directory containing Singularity file
sudo singularity build bartender.sif Singularity
```

## Running bartender extractor
```bash
singularity exec \
    -B /data/SBGE/cory/1011-barcodes:/home/wellerca/ \
    bartender/bartender.sif bartender_extractor_com \
    -f seq/MS3059950-600V3_19109498_S7_L001_R1_001.fastq \
    -o pre  \
    -p CGAGC[34]C -m 1
```

Output:
```
Running bartender extractor
bartender_extractor seq/MS3059950-600V3_19109498_S7_L001_R1_001.fastq pre 1 "(CGAG.|CGA.C|CG.GC|C.AGC|.GAGC)([ATCGN]{34})(C)" CGAGC C 3 1
Totally there are 1187764 reads in seq/MS3059950-600V3_19109498_S7_L001_R1_001.fastq file!
Totally there are 1118562 valid barcodes from seq/MS3059950-600V3_19109498_S7_L001_R1_001.fastq file
Totally there are 1118562 valid barcodes whose quality pass the quality condition
The estimated sequence error from the prefix and suffix parts is 0.0311966
```


## Formatting barcodes
The `extracted_barcode.txt` file contains a 34-mer nucleotide sequence, but we only
want the 20 nucleotide barcode sequence contained within.
```python3
python3 format_barcodes.py pre_barcode.txt > barcodes.txt
```

## Running bartender cluster
```bash
singularity exec \
    -B /data/SBGE/cory/1011-barcodes:/home/wellerca/ \
    bartender/bartender.sif bartender_single_com  \
    -f barcodes.txt \
    -o barcode_clusters  \
    -d 2 \
    -s 5
```

output:
```
Running bartender
Loading barcodes from the file
It takes 00:00:01 to load the barcodes from barcodes.txt
Shortest barcode length: 20
Longest barcode length: 20
Start to group barcode with length 20
Using two sample unpooled test
Transforming the barcodes into seed clusters
Initial number of unique reads:  64431
The distance threshold is 2
Clustering iteration 1
Clustering iteration 2
Clustering iteration 3
Clustering iteration 4
Identified 18272 barcodes with length 20
The clustering process takes 00:00:01
Start to dump clusters to file with prefix barcode_clusters
Start to remove pcr effects
***(Overall error rate estimated from the clustering result)***
Total number of clusters after removing PCR effects: 18272
The estimated error rate is 0.00340786
The overall running time 00:00:05 seconds.
```


## Take most abundant seq (consensus) per cluster and plot
```R
library(data.table)
library(ggplot2)
library(ggrepel)

dat <- fread('barcode_clusters_barcodes.csv')

consensus <- dat[, .SD[which.max(Frequency)], by=Cluster.ID]

setnames(consensus, "Unique.reads", "consensus")
consensus[, Frequency := NULL]
setkey(consensus, Cluster.ID)
setkey(dat, Cluster.ID)

dat.merge <- merge(dat, consensus)
consensus_counts <- dat.merge[, list("N" = sum(Frequency)), by=consensus][order(-N)]
consensus_counts[, abundance_rank := 1:.N]

# ggplot(consensus_counts, aes(x=abundance_rank, y=N)) + geom_point() +
# scale_y_continuous(trans='log10', 
#                     breaks=c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6),
#                     labels=c("1", "10", "100", "1000", "1000", "10000", "100000")) +
# labs(x="Barcodes ranked by abundance",
#         y="Abundance") +
# theme_few(12)

consensus_counts[N > 3000, text_label := consensus]


ggplot(consensus_counts[N>=2], aes(x=abundance_rank, y=N)) + geom_point() +
scale_y_continuous(trans='log10', 
                    breaks=c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6),
                    labels=c("1", "10", "100", "1000", "10000", "100000", "1000000")) +
labs(x="Barcodes ranked by abundance",
        y="Abundance",
        title="Barcodes with cluster counts >= 2") +
theme_few(12) +
geom_text_repel(aes(label=text_label))

ggplot(consensus_counts[N>=2 & N < 10000], aes(x=abundance_rank, y=N)) + geom_point() +
labs(x="Barcodes ranked by abundance",
        y="Abundance",
        title="Barcodes with cluster counts >= 2 and <= 100000") +
theme_few(12) +
geom_text_repel(aes(label=text_label))

```
