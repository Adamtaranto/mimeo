# LASTZ2Repeat

Detection of high copy number elements within OR between species by whole-genome
alignment with LASTZ.

## Prepare input

### Split genomes into individual fasta by scaffold

Use [subset-fasta-by-name](https://github.com/Adamtaranto/subset-fasta-by-name)
script to split multifasta.

```bash
./SubsetFastaByName.py -i genome.fa --splitMode -d genome_split
```

### Get scaffold lengths

Requires samtools.

```bash
samtools faidx genome.fa
cut -f1,2 genome.fa.fai > chromlens.txt
```

## Find repeats **WITHIN** a genome

**Arguments:**
  > -z Path to LASTZ executable
  > -g Directory containing scaffold fasta files
  > -i Min identity threshold
  > -p Output file prefix
  > -l Min repeat length
  > -d Output directory
  > -C Min hits from all non-self scaffolds to call repeat
  > -c Min hits from within same scaffold to call repeat
  > -L Scaffold length file
  > -s Skip comparison of scaffolds with same name

### Within genome example usage

Search scaffolds for features with > 3X coverage across all other genomic scaffolds,
or 4X coverage of hits within same scaffold.
Hits must have > 80% identity, and be > 100 bp in length.

```bash
./Self_genome_align_repeat_finder.sh \
-z /usr/local/bin/lastz \
-g genome_split \
-i 80 \
-p HighCopy_genomeA_id70_cv4 \
-l 100 \
-d outdir \
-C 3 \
-c 4 \
-L chromlens.txt
```

## Find features that are high-copy in a **non-self** genome

**Arguments:**
  > -z Path to LASTZ executable
  > -t Target scaffold directory
  > -q Query scaffold directory
  > -i Min identity threshold
  > -p Output file prefix
  > -l Min repeat length
  > -d Output directory
  > -c Hit coverage threshold
  > -L Scaffold length file
  > -s Skip comparison of scaffolds with same name

### Non-self genome example usage

Search target scaffolds for features with at least 5X coverage by matches present
in query scaffold set.
Matches between genomes must have > 60% identity, and be > 100 bp in length.

```bash
./Interspecies_Repeat_Scan.sh \
-z /usr/local/bin/lastz \
-t target_split \
-q query_split \
-i 60 \
-p query2target_id60_cov5 \
-l 100 \
-d outdir \
-c 5 \
-L chromlens.txt
```
