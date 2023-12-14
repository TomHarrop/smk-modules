References were downloaded from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001215.4 and subset to Chromosome 2L

Reads were downloaded from zenodo

```bash
wget -O - https://zenodo.org/record/6457007/files/GSM461177_1_subsampled.fastqsanger | pigz -9 > GSM461177_subsampled.r1.fq.gz &
wget -O - https://zenodo.org/record/6457007/files/GSM461177_2_subsampled.fastqsanger | pigz -9 > GSM461177_subsampled.r2.fq.gz &
wget -O - https://zenodo.org/record/6457007/files/GSM461180_1_subsampled.fastqsanger | pigz -9 > GSM461180_subsampled.r1.fq.gz & 
wget -O - https://zenodo.org/record/6457007/files/GSM461180_2_subsampled.fastqsanger | pigz -9 > GSM461180_subsampled.r2.fq.gz &
wait
```