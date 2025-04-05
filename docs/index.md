# Mimeo

## Table of contents

* [Modules](#modules)
* [Installing Mimeo](#installing-mimeo)
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

## License

Software provided under MIT license.
