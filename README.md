[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI version](https://badge.fury.io/py/mimeo.svg)](https://badge.fury.io/py/mimeo)

# Mimeo

**A tool for finding and annotating repeats in whole-genome alignments.**

## Table of contents

* [Modules](#modules)
* [Installing Mimeo](#installing-mimeo)
* [Example usage](#example-usage)
* [Standard options](#standard-options)
  * [mimeo-self](#mimeo-self)
  * [mimeo-x](#mimeo-x)
  * [mimeo-map](#mimeo-map)
  * [mimeo-filter](#mimeo-filter)
* [Alternative alignment engines](#importing-alignments)
* [License](#license)

## Modules

Mimeo comprises three tools for parsing repeats from whole-genome alignments:

### mimeo-self

**Internal repeat finder.** Mimeo-self aligns a genome to itself and extracts high-identity segments above
a coverage threshold. This method is less sensitive to disruption by indels and repeat-directed point mutations than
kmer-based methods such as RepeatScout. Reported annotations indicate overlapping segments above the coverage threshold,
mimeo-self does not attempt to separate nested repeats. Use this tool to identify candidate repeat regions for curated annotation.

### mimeo-x

**Cross-species repeat finder.** A newly acquired or low-copy transposon may slip past copy-number based annotation tools. Mimeo-x searches for features which are abundant in an external reference genome, allowing for
annotation of complete elements as they occur in a horizontal-transfer donor species, or of conserved coding segments
of related transposon families.

### mimeo-map

**Find all high-identity segments shared between genomes.** Mimeo-map identifies candidate horizontally
transferred segments between sufficiently diverged species. When comparing isolates of a single species, aligned segments correspond to directly homologous sequences and internally repetitive features.

Intra/Inter-genomic alignments from Mimeo-self or Mimeo-x can be reprocessed with Mimeo-map to generate annotations of
unfiltered/uncollapsed alignments. These raw alignment annotations can be used to interrogate repetitive-segments for coverage breakpoints corresponding to nested transposons with differing abundances across the genome.

### mimeo-filter

An additional tool **mimeo-filter** is now included to allow post-filtering of SSR-rich sequences from FASTA formatted
candidate-repeat libraries.

## Installing Mimeo

Requirements:

* [LASTZ](http://www.bx.psu.edu/~rsharris/lastz/) genome alignment tool from the Miller Lab, Penn State.
* [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html)
* [trf](https://tandem.bu.edu/trf/trf.html)

Install from Bioconda:

```bash
conda install mimeo
```

Install from PyPi:

```bash
pip install mimeo
```

Clone and install from this repository:

```bash
git clone https://github.com/Adamtaranto/mimeo.git && cd mimeo

pip install -e '.[dev]'
```

## Example usage

### Demo: mimeo-self

Annotate features in genome A which are > 100bp and occur with >=
80% identity at least 3 times on other scaffolds OR at least 4 times
on the same scaffold.

```bash
mimeo self --adir data/A_genome_Split --afasta data/A_genome.fasta \
-d MS_outdir --gffout A_genome_Inter3_Intra4_id80_len_100.gff3 \
--outfile A_genome_Self_Align.tab --label A_Rep3 --prefix A_Self --minIdt 80 \
--minLen 100 --minCov 3 --intraCov 4 --strictSelf
```

Output:

* MS_outdir/A_genome_Inter3_Intra4_id80_len_100.gff3
* MS_outdir/A_genome_Self_Align.tab
* data/A_genome_Split/*.fa

### Demo: mimeo-x

Annotate features in genome A which are > 100bp and occur with >=
80% identity at least 5 times in genome B.

```bash
mimeo x --afasta data/A_genome.fasta --bfasta data/B_genome.fasta \
-d MX_outdir --gffout B_Rep5_in_A.gff3 --outfile B_Reps_in_A_id80_len100.tab \
--label B_Rep5 --prefix B_Rep5 --minIdt 80 --minLen 100 --minCov 5
```

Output:

* MX_outdir/B_Rep5_in_A.gff3
* MX_outdir/B_Reps_in_A_id80_len100.tab

### Demo: mimeo-map

Annotate features in genome A which are > 100bp and occur with >=
90% identity in genome B. No coverage filter, all alignments are reported.

```bash
mimeo map --afasta data/A_genome.fasta --bfasta data/B_genome.fasta \
-d MM_outdir --gffout B_in_A_id90.gff3 --outfile B_in_A_id90.tab \
--label B_90 --prefix B_90 --minIdt 90 --minLen 100
```

Output:

* MM_outdir/B_in_A_id90.gff3
* MM_outdir/B_in_A_id90.tab

### mimeo-map + SSR filter

Annotate features in genome A which are > 100bp and occur with >=
98% identity in genome B. Reuse B to A-genome alignment from the previous run.

Filter out hits which are >= 40% tandem repeats. Write filtered hits
as tab file and GFF3 annotation.

```bash
mimeo map --afasta data/A_genome.fasta --bfasta data/B_genome.fasta \
-d MM_outdir --gffout B_in_A_id98_maxSSR40.gff3 --outfile B_in_A_id98.tab \
--label B_98 --prefix B_98 --minIdt 98 --minLen 100 \
--recycle --maxtandem 40 --writeTRF
```

Output:

* MM_outdir/B_in_A_id98_maxSSR40.gff3
* MM_outdir/B_in_A_id98.tab.trf

### Demo: mimeo-filter

Filter sequences comprised of >= 40% short tandem repeats from a multifasta
library of candidate transposons.

```bash
mimeo filter --infile data/candidate_TEs.fa
```

Output:

* candidate_TEs_filtered.fa

## Standard options

### mimeo-self

```code
Usage: mimeo self [-h] [--adir ADIR] [--afasta AFASTA] [-r] [-d OUTDIR]
                  [--gffout GFFOUT] [--outfile OUTFILE] [--verbose]
                  [--label LABEL] [--prefix PREFIX] [--lzpath LZPATH]
                  [--bedtools BEDTOOLS] [--minIdt MINIDT] [--minLen MINLEN]
                  [--minCov MINCOV] [--hspthresh HSPTHRESH]
                  [--intraCov INTRACOV] [--strictSelf]

Internal repeat finder. Mimeo-self aligns a genome to itself and extracts
high-identity segments above a coverage threshold.

Optional arguments:
  -h, --help      Show this help message and exit.
  --adir          Name of the directory containing sequences from the genome.
                  Write split files here if providing genome as
                  multifasta.
  --afasta        Genome as multifasta.
  -r, --recycle   Use existing alignment "--outfile" if found.
  -d , --outdir   Write output files to this directory. (Default: cwd)
  --gffout        Name of GFF3 annotation file.
  --outfile       Name of alignment result file.
  --verbose       If set report LASTZ progress.
  --label         Set annotation TYPE field in gff.
  --prefix        ID prefix for internal repeats.
  --lzpath        Custom path to LASTZ executable if not in $PATH.
  --bedtools      Custom path to bedtools executable if not in $PATH.
  --minIdt        Minimum alignment identity to report.
  --minLen        Minimum alignment length to report.
  --minCov        Minimum depth of aligned segments to report repeat
                  feature.
  --hspthresh     Set HSP min score threshold for LASTZ.
  --intraCov      Minimum depth of aligned segments from the same scaffold
                  to report feature. Used if "--strictSelf" mode is
                  selected.
  --strictSelf    If set process same-scaffold alignments separately
                  with the option to use higher "--intraCov" threshold.
                  Sometimes useful to avoid false repeat calls from
                  staggered alignments over SSRs or short tandem
                  duplication.
```

### mimeo-x

```code
Usage: mimeo x [-h] [--adir ADIR] [--bdir BDIR] [--afasta AFASTA]
               [--bfasta BFASTA] [-r] [-d OUTDIR] [--gffout GFFOUT]
               [--outfile OUTFILE] [--verbose] [--label LABEL]
               [--prefix PREFIX] [--lzpath LZPATH] [--bedtools BEDTOOLS]
               [--minIdt MINIDT] [--minLen MINLEN] [--minCov MINCOV]
               [--hspthresh HSPTHRESH]

Cross-species repeat finder. Mimeo-x searches for features which are abundant
in an external reference genome.

Optional arguments:
  -h, --help      Show this help message and exit.
  --adir          Name of the directory containing sequences from A genome.
  --bdir          Name of the directory containing sequences from B genome.
  --afasta        A genome as multifasta.
  --bfasta        B genome as multifasta.
  -r, --recycle   Use existing alignment "--outfile" if found.
  -d , --outdir   Write output files to this directory. (Default: cwd)
  --gffout        Name of GFF3 annotation file.
  --outfile       Name of alignment result file.
  --verbose       If set report LASTZ progress.
  --label         Set annotation TYPE field in GFF.
  --prefix        ID prefix for B-genome repeats annotated in A-genome.
  --lzpath        Custom path to LASTZ executable if not in $PATH.
  --bedtools      Custom path to bedtools executable if not in $PATH.
  --minIdt        Minimum alignment identity to report.
  --minLen        Minimum alignment length to report.
  --minCov        Minimum depth of B-genome hits to report feature in
                  A-genome.
  --hspthresh     Set HSP min score threshold for LASTZ.
```

### mimeo-map

```code
Usage: mimeo map [-h] [--adir ADIR] [--bdir BDIR] [--afasta AFASTA]
                 [--bfasta BFASTA] [-r] [-d OUTDIR] [--gffout GFFOUT]
                 [--outfile OUTFILE] [--verbose] [--label LABEL]
                 [--prefix PREFIX] [--keeptemp] [--lzpath LZPATH]
                 [--minIdt MINIDT] [--minLen MINLEN] [--hspthresh HSPTHRESH]
                 [--TRFpath TRFPATH] [--tmatch TMATCH] [--tmismatch TMISMATCH]
                 [--tdelta TDELTA] [--tPM TPM] [--tPI TPI]
                 [--tminscore TMINSCORE] [--tmaxperiod TMAXPERIOD]
                 [--maxtandem MAXTANDEM] [--writeTRF]

Find all high-identity segments shared between genomes.

Optional arguments:
  -h, --help      Show this help message and exit.
  --adir          Name of the directory containing sequences from A genome.
  --bdir          Name of the directory containing sequences from B genome.
  --afasta        A genome as multifasta.
  --bfasta        B genome as multifasta.
  -r, --recycle   Use existing alignment "--outfile" if found.
  -d, --outdir    Write output files to this directory. (Default: cwd)
  --gffout        Name of GFF3 annotation file. If not set, suppress
                  output.
  --outfile       Name of alignment result file.
  --verbose       If set report LASTZ progress.
  --label         Set annotation TYPE field in GFF.
  --prefix        ID prefix for B-genome hits annotated in A-genome.
  --keeptemp      If set does not remove temp files.
  --lzpath        Custom path to LASTZ executable if not in $PATH.
  --minIdt        Minimum alignment identity to report.
  --minLen        Minimum alignment length to report.
  --hspthresh     Set HSP min score threshold for LASTZ.
  --TRFpath       Custom path to TRF executable if not in $PATH.
  --tmatch        TRF matching weight.
  --tmismatch     TRF mismatching penalty.
  --tdelta        TRF indel penalty.
  --tPM           TRF match probability.
  --tPI           TRF indel probability.
  --tminscore     TRF minimum alignment score to report.
  --tmaxperiod    TRF maximum period size to report.
  --maxtandem     Max percentage of an A-genome alignment which may be masked by TRF.
                  If exceeded, the alignment will be discarded.
  --writeTRF      If set write TRF filtered alignment file for use with
                  other mimeo modules.
```

### mimeo-filter

```code
Usage: mimeo filter [-h] --infile INFILE [-d OUTDIR] [--outfile OUTFILE]
                    [--keeptemp] [--verbose] [--TRFpath TRFPATH]
                    [--tmatch TMATCH] [--tmismatch TMISMATCH]
                    [--tdelta TDELTA] [--tPM TPM] [--tPI TPI]
                    [--tminscore TMINSCORE] [--tmaxperiod TMAXPERIOD]
                    [--maxtandem MAXTANDEM]

Filter SSR containing sequences from FASTA library of repeats.

Optional arguments:
  -h, --help            Show this help message and exit.
  --infile              Name of the directory containing sequences from A genome.
  -d, --outdir          Write output files to this directory. (Default: cwd)
  --outfile             Name of alignment result file.
  --keeptemp            If set does not remove temp files.
  --verbose             If set report LASTZ progress.
  --TRFpath             Custom path to TRF executable if not in $PATH.
  --tmatch              TRF matching weight
  --tmismatch           TRF mismatching penalty.
  --tdelta              TRF indel penalty.
  --tPM                 TRF match probability.
  --tPI                 TRF indel probability.
  --tminscore           TRF minimum alignment score to report.
  --tmaxperiod          TRF maximum period size to report. Note: Setting this
                        score too high may exclude some LTR retrotransposons.
                        Optimal len to exclude only SSRs is 10-50bp.
  --maxtandem           Max percentage of a sequence which may be masked by
                        TRF. If exceeded, the element will be discarded.

```

## Importing alignments

Whole genome alignments generated by alternative tools (i.e. BLAT) can be provided to any of the Mimeo modules
as a tab-delimited file with the columns:

```code
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

## License

Software provided under MIT license.
