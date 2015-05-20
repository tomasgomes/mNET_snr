# mNET_snr

Python script that, given an aligned bam file and sample name, returns a bam file of the last base from read 2 (of a read pair), with strand information of read 1. Designed to be used with mNET-seq data, although it can also work with other types of data (if the protocol produces a pair-end secondstranded library).

Described as part of an analysis pipeline in *insert reference*

For usage, see get_SNR_bam.py -h
