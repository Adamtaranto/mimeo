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

**mimeo-self**: Internal repeat finder. Mimeo-self aligns a genome to itself and extracts high-identity segments above
an coverage threshold. This method is less sensitive to disruption by indels and repeat-directed point mutations than
kmer-based methods such as RepeatScout. Reported annotations indicate overlaping segments above the coverage threshold,
mimeo-self does not attempt to separate nested repeats. Use this tool to identify candidate repeat regions for curated 
annotation.

**mimeo-x**: Cross-species repeat finder. A newly acquired or low-copy transposon may slip past copy-numner based 
annotation tools. Mimeo-x searches for features which are abundant in an external reference genome, allowing for
annotation of complete elements as they occur in a horizontal-transfer donor species, or of conserved coding segements
of related transposon families.

**mimeo-map**: Find all high-identity segments shared between genomes. Mimeo-map identifies candidate horizontally
transferred segments between sufficiently diverged species. When comparing isolates of a single species, aligned 
segments correspond to directly homologous sequences and internally repetative features.  


Intra/Inter-genomic alignments from Mimeo-self or Mimeo-x can be reprocessed with Mimeo-map to generate annotations of
unfiltered/uncollapsed alignments. These raw alignment annotations can be used to interrogate repetitive-segments for 
coverage breakpoints corresponding to nested transposons with differing abundances across the genome.


# Getting started

### Installing Mimeo

Requirements: 
  * [LASTZ](http://www.bx.psu.edu/~rsharris/lastz/) genome alignment tool from the Miller Lab, Penn State.
  * [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html)

Install from PyPi:
```
pip install mimeo
```

Clone and install from this repository:
```
git clone https://github.com/Adamtaranto/mimeo.git && cd mimeo && pip install -e .
```

### Example usage 

**Find repetitive features within a genome**  

```
mimeo-self 
```
Output: 
xxx

****

**Annotate features which are repetative in a non-self genome**

```
mimeo-x 
```
Output: 
xxx

**Find high-identity segments shared between genomes**

```
mimeo-map
```
Output: 
xxx




### Standard options

#### mimeo-self

```
```

#### mimeo-x

```
```

#### mimeo-map

```
```

# Importing alignments

Whole genome alignments generated by alternative tool (i.e. BLAT) can be provided to any of the Mimeo modules
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



