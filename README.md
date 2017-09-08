# Mimeo

Scan genomes for internally repeated sequences, elements which are repetitive in another species, or high-identity HGT candidate regions between species.  

# Table of contents

* [Modules](#modules)
* [Getting started](#getting-started)
    * [Installing Mimeo](#installing-mimeo)
    * [Example usage](#example-usage)
* [Standard options](#standard-options)
  *  [mimeo-self](#mimeo-self)
  *  [mimeo-x](#mimeo-x)
  *  [mimeo-map](#mimeo-map)
* [Alternative alignment engines](#importing-alignments)
* [License](#license)

# Modules

Mimeo comprises three tools for processing whole-genome alignments:

## mimeo-self  

**Internal repeat finder.** Mimeo-self aligns a genome to itself and extracts high-identity segments above
an coverage threshold. This method is less sensitive to disruption by indels and repeat-directed point mutations than
kmer-based methods such as RepeatScout. Reported annotations indicate overlaping segments above the coverage threshold,
mimeo-self does not attempt to separate nested repeats. Use this tool to identify candidate repeat regions for curated 
annotation.

## mimeo-x  

**Cross-species repeat finder.** A newly acquired or low-copy transposon may slip past copy-numner based 
annotation tools. Mimeo-x searches for features which are abundant in an external reference genome, allowing for
annotation of complete elements as they occur in a horizontal-transfer donor species, or of conserved coding segements
of related transposon families.

## mimeo-map  

**Find all high-identity segments shared between genomes.** Mimeo-map identifies candidate horizontally
transferred segments between sufficiently diverged species. When comparing isolates of a single species, aligned 
segments correspond to directly homologous sequences and internally repetative features.  


Intra/Inter-genomic alignments from Mimeo-self or Mimeo-x can be reprocessed with Mimeo-map to generate annotations of
unfiltered/uncollapsed alignments. These raw alignment annotations can be used to interrogate repetitive-segments for 
coverage breakpoints corresponding to nested transposons with differing abundances across the genome.


# Getting started

## Installing Mimeo

Requirements: 
  * [LASTZ](http://www.bx.psu.edu/~rsharris/lastz/) genome alignment tool from the Miller Lab, Penn State.
  * [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html)

Install from PyPi:
```bash
pip install mimeo
```

Clone and install from this repository:
```bash
git clone https://github.com/Adamtaranto/mimeo.git && cd mimeo && pip install -e .
```

## Example usage 

### mimeo-self

Annotate features in genome A which are > 100bp and occur with >=
80% identity at least 3 times on other scaffolds OR at least 4 times
on the same scaffold.

```bash
mimeo-self --adir data/A_genome_Split --afasta data/A_genome.fasta \
-d MS_outdir --gffout A_genome_Inter3_Intra4_id80_len_100.gff3 \
--outfile A_genome_Self_Align.tab --label A_Rep3 --prefix A_Self --minIdt 80 \
--minLen 100 --minCov 3 --intraCov 4 --strictSelf
```

Output: 
  - MS_outdir/A_genome_Inter3_Intra4_id80_len_100.gff3
  - MS_outdir/A_genome_Self_Align.tab
  - data/A_genome_Split/*.fa

### mimeo-x

Annotate features in genome A which are > 100bp and occur with >=
80% identity at least 5 times in genome B.

```bash
mimeo-x --afasta data/A_genome.fasta --bfasta data/B_genome.fasta \
-d MX_outdir --gffout B_Rep5_in_A.gff3 --outfile B_Reps_in_A_id80_len100.tab \
--label B_Rep5 --prefix B_Rep5 --minIdt 80 --minLen 100 --minCov 5
```

Output: 
  - MX_outdir/B_Rep5_in_A.gff3
  - MX_outdir/B_Reps_in_A_id80_len100.tab

### mimeo-map

Annotate features in genome A which are > 100bp and occur with >=
98% identity in genome B. No coverage filter, all alignments are reported.

```bash
mimeo-map --afasta data/A_genome.fasta --bfasta data/B_genome.fasta \
-d MM_outdir --gffout B_in_A_id98.gff3 --outfile B_in_A_id98.tab \
--label B_98 --prefix B_98 --minIdt 98 --minLen 100
```

Output: 
  - MM_outdir/B_in_A_id98.gff3
  - MM_outdir/B_in_A_id98.tab


## Standard options

### mimeo-self

```
usage: mimeo-self [-h] [--adir ADIR] [--afasta AFASTA] [-r] [-d OUTDIR]
                  [--gffout GFFOUT] [--outfile OUTFILE] [--verbose]
                  [--label LABEL] [--prefix PREFIX] [--lzpath LZPATH]
                  [--bedtools BEDTOOLS] [--minIdt MINIDT] [--minLen MINLEN]
                  [--minCov MINCOV] [--hspthresh HSPTHRESH]
                  [--intraCov INTRACOV] [--strictSelf]

Internal repeat finder. Mimeo-self aligns a genome to itself and extracts
high-identity segments above an coverage threshold.

optional arguments:
  -h, --help            show this help message and exit
  --adir ADIR           Name of directory containing sequences from genome.
                        Write split files here if providing genome as
                        multifasta.
  --afasta AFASTA       Genome as multifasta.
  -r, --recycle         Use existing alignment "--outfile" if found.
  -d OUTDIR, --outdir OUTDIR
                        Write output files to this directory. (Default: cwd)
  --gffout GFFOUT       Name of GFF3 annotation file.
  --outfile OUTFILE     Name of alignment result file.
  --verbose             If set report LASTZ progress.
  --label LABEL         Set annotation TYPE field in gff.
  --prefix PREFIX       ID prefix for internal repeats.
  --lzpath LZPATH       Custom path to LASTZ executable if not in $PATH.
  --bedtools BEDTOOLS   Custom path to bedtools executable if not in $PATH.
  --minIdt MINIDT       Minimum alignment identity to report.
  --minLen MINLEN       Minimum alignment length to report.
  --minCov MINCOV       Minimum depth of aligned segments to report repeat
                        feature.
  --hspthresh HSPTHRESH
                        Set HSP min score threshold for LASTZ.
  --intraCov INTRACOV   Minimum depth of aligned segments from same scaffold
                        to report feature. Used if "--strictSelf" mode is
                        selected.
  --strictSelf          If set process same-scaffold alignments separately
                        with option to use higher "--intraCov" threshold.
                        Sometime useful to avoid false repeat calls from
                        staggered alignments over SSRs or short tandem
                        duplication.
```

### mimeo-x

```
usage: mimeo-x [-h] [--adir ADIR] [--bdir BDIR] [--afasta AFASTA]
               [--bfasta BFASTA] [-r] [-d OUTDIR] [--gffout GFFOUT]
               [--outfile OUTFILE] [--verbose] [--label LABEL]
               [--prefix PREFIX] [--lzpath LZPATH] [--bedtools BEDTOOLS]
               [--minIdt MINIDT] [--minLen MINLEN] [--minCov MINCOV]
               [--hspthresh HSPTHRESH]

Cross-species repeat finder. Mimeo-x searches for features which are abundant
in an external reference genome.

optional arguments:
  -h, --help            show this help message and exit
  --adir ADIR           Name of directory containing sequences from A genome.
  --bdir BDIR           Name of directory containing sequences from B genome.
  --afasta AFASTA       A genome as multifasta.
  --bfasta BFASTA       B genome as multifasta.
  -r, --recycle         Use existing alignment "--outfile" if found.
  -d OUTDIR, --outdir OUTDIR
                        Write output files to this directory. (Default: cwd)
  --gffout GFFOUT       Name of GFF3 annotation file.
  --outfile OUTFILE     Name of alignment result file.
  --verbose             If set report LASTZ progress.
  --label LABEL         Set annotation TYPE field in gff.
  --prefix PREFIX       ID prefix for B-genome repeats annotated in A-genome.
  --lzpath LZPATH       Custom path to LASTZ executable if not in $PATH.
  --bedtools BEDTOOLS   Custom path to bedtools executable if not in $PATH.
  --minIdt MINIDT       Minimum alignment identity to report.
  --minLen MINLEN       Minimum alignment length to report.
  --minCov MINCOV       Minimum depth of B-genome hits to report feature in
                        A-genome.
  --hspthresh HSPTHRESH
                        Set HSP min score threshold for LASTZ.
```

### mimeo-map

```
usage: mimeo-map [-h] [--adir ADIR] [--bdir BDIR] [--afasta AFASTA]
                 [--bfasta BFASTA] [-r] [-d OUTDIR] [--gffout GFFOUT]
                 [--outfile OUTFILE] [--verbose] [--label LABEL]
                 [--prefix PREFIX] [--lzpath LZPATH] [--minIdt MINIDT]
                 [--minLen MINLEN] [--hspthresh HSPTHRESH]

Find all high-identity segments shared between genomes.

optional arguments:
  -h, --help            show this help message and exit
  --adir ADIR           Name of directory containing sequences from A genome.
  --bdir BDIR           Name of directory containing sequences from B genome.
  --afasta AFASTA       A genome as multifasta.
  --bfasta BFASTA       B genome as multifasta.
  -r, --recycle         Use existing alignment "--outfile" if found.
  -d OUTDIR, --outdir OUTDIR
                        Write output files to this directory. (Default: cwd)
  --gffout GFFOUT       Name of GFF3 annotation file.
  --outfile OUTFILE     Name of alignment result file.
  --verbose             If set report LASTZ progress.
  --label LABEL         Set annotation TYPE field in gff.
  --prefix PREFIX       ID prefix for B-genome hits annotated in A-genome.
  --lzpath LZPATH       Custom path to LASTZ executable if not in $PATH.
  --minIdt MINIDT       Minimum alignment identity to report.
  --minLen MINLEN       Minimum alignment length to report.
  --hspthresh HSPTHRESH
                        Set HSP min score threshold for LASTZ.
```

# Importing alignments

Whole genome alignments generated by alternative tools (i.e. BLAT) can be provided to any of the Mimeo modules
as a tab-delimited file with the columns:

```
[1]   name1     = Name of target sequence in genome A
[2]   strand1   = Strand of alignment in target sequence
[3]   start1    = 5-prime position of alignment in target (lower value irrespective of strand)
[4]   end1      = 3-prime position of alignment in target (higher value irrespective of strand)
[5]   name2     = Name of source sequence in genome B
[6]   strand2   = Strand of alignment in source
[7]   start2+   = 5-prime position of alignment in source (lower value irrespective of strand)
[8]   end2+     = 3-prime position of alignment in source (higher value irrespective of strand)
[9]   score     = Alignment score as int
[10]  identity  = Identity of alignment as float
```

File should be sorted by columns 1,3,4

# License

Software provided under MIT license.



