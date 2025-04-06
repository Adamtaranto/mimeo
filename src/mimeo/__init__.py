"""
Mimeo: Toolkit for horizontal gene transfer and repetitive element analysis in genomic data.

Mimeo provides a suite of tools for identifying and analyzing horizontal gene transfer events,
repetitive elements, and sequence homology in genomic data. The package offers:

- High-identity segment detection between genomes
- Self-alignment analysis for repetitive element discovery
- Integration with common genomic alignment tools like LASTZ
- Tandem repeat filtering using TRF (Tandem Repeats Finder)
- Customizable alignment parameters and filtering thresholds
- Output in standard formats (GFF3, BED, tabular)

The package is designed for comparative genomics workflows and supports analysis
of both closely and distantly related genomes. Mimeo can identify potential horizontal
gene transfer events, characterize repetitive elements, and quantify sequence homology
across genomes.

Main components:
- utils: Core utilities for sequence handling and system operations
- wrappers: Interface to external tools like LASTZ and TRF
- run_self: Self-alignment analysis for repetitive element discovery
- run_xspecies: Cross-species sequence comparison

For detailed usage information, see the documentation or run the command-line tools
with the --help flag.
"""
