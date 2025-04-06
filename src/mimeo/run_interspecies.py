"""
Cross-species repeat identification through comparative genomics.

This module provides command-line functionality for identifying repetitive elements
that are abundant in an external reference genome (B) and present in the target
genome (A). It detects segments where regions from genome B align multiple times
to genome A, indicating potentially mobile or repetitive elements.

The tool can:
1. Split multiFASTA files into individual sequence files if needed
2. Perform cross-species alignments using LASTZ
3. Identify regions in genome A covered by multiple genome B segments
4. Filter results based on identity, length, and coverage thresholds
5. Output results in both tabular and GFF3 formats

This is part of the mimeo package for comparative genomics analysis.
"""

import argparse
import logging
import os
import shutil
import sys
from typing import Dict, List

from ._version import __version__
from .logs import init_logging
from .utils import (
    chromlens,
    get_all_pairs,
    missing_tool,
    run_cmd,
    set_paths,
)
from .wrappers import xspecies_LZ_cmds


def mainArgs() -> argparse.Namespace:
    """
    Parse command-line arguments for the mimeo-x tool.

    Sets up the argument parser with options for input files/directories,
    output formats, alignment parameters, and thresholds for identifying
    regions in genome A that are covered by multiple segments from genome B.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Cross-species repeat finder. Mimeo-x searches for features which are abundant in an external reference genome.',
        prog='mimeo-x',
    )
    # Add version information
    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {__version__}',
        help='Show program version and exit.',
    )
    # Input options
    parser.add_argument(
        '--adir',
        type=str,
        default=None,
        help='Name of directory containing sequences from A genome.',
    )
    parser.add_argument(
        '--bdir',
        type=str,
        default=None,
        help='Name of directory containing sequences from B genome.',
    )
    parser.add_argument(
        '--afasta', type=str, default=None, help='A genome as multifasta.'
    )
    parser.add_argument(
        '--bfasta', type=str, default=None, help='B genome as multifasta.'
    )
    parser.add_argument(
        '-r',
        '--recycle',
        action='store_true',
        help='Use existing alignment "--outfile" if found.',
    )
    # Output options
    parser.add_argument(
        '-d',
        '--outdir',
        type=str,
        default=None,
        help='Write output files to this directory. (Default: cwd)',
    )
    parser.add_argument(
        '--gffout',
        type=str,
        default='mimeo_B_in_A.gff3',
        help='Name of GFF3 annotation file.',
    )
    parser.add_argument(
        '--outfile',
        type=str,
        default='mimeo_alignment.tab',
        help='Name of alignment result file.',
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        default=False,
        help='If set report LASTZ progress.',
    )
    parser.add_argument(
        '--label',
        type=str,
        default='B_Repeat',
        help='Set annotation TYPE field in gff.',
    )
    parser.add_argument(
        '--prefix',
        type=str,
        default='B_Repeat',
        help='ID prefix for B-genome repeats annotated in A-genome.',
    )
    parser.add_argument(
        '--keeptemp',
        action='store_true',
        default=False,
        help='If set do not remove temp files.',
    )
    # Alignment options
    parser.add_argument(
        '--lzpath',
        type=str,
        default='lastz',
        help='Custom path to LASTZ executable if not in $PATH.',
    )
    parser.add_argument(
        '--bedtools',
        type=str,
        default='bedtools',
        help='Custom path to bedtools executable if not in $PATH.',
    )
    parser.add_argument(
        '--minIdt', type=int, default=60, help='Minimum alignment identity to report.'
    )
    parser.add_argument(
        '--minLen', type=int, default=100, help='Minimum alignment length to report.'
    )
    parser.add_argument(
        '--minCov',
        type=int,
        default=5,
        help='Minimum depth of B-genome hits to report feature in A-genome.',
    )
    parser.add_argument(
        '--hspthresh',
        type=int,
        default=3000,
        help='Set HSP min score threshold for LASTZ.',
    )
    parser.add_argument(
        '--loglevel',
        type=str,
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        help='Set the logging level.',
    )
    args = parser.parse_args()
    return args


def main() -> None:
    """
    Execute the main cross-species repeat identification workflow.

    This function:
    1. Validates required external tools
    2. Sets up output paths
    3. Identifies sequence files to align between genomes A and B
    4. Calculates chromosome lengths for coverage calculations
    5. Runs LASTZ alignments (or reuses existing alignments)
    6. Processes alignments to identify regions in genome A covered by multiple
       genome B segments, which may represent repetitive elements
    7. Outputs results in tabular and GFF3 formats

    Returns
    -------
    None
        This function does not return a value. It performs file I/O and logging
        to report progress and results.
    """
    # Get command line arguments
    args = mainArgs()

    # Check for required external programs
    tools: List[str] = [args.lzpath, args.bedtools]
    missing_tools: List[str] = []
    for tool in tools:
        missing_tools += missing_tool(tool)
    if missing_tools:
        print(
            'WARNING: Some tools required by mimeo could not be found: '
            + ', '.join(missing_tools),
            file=sys.stderr,
        )
        print('You may need to install them to use all features.', file=sys.stderr)

    # Initialize logging
    init_logging(loglevel=args.loglevel)
    logging.info('Starting cross-species repeat identification...')
    logging.debug('Command line arguments: %s', args)

    # Set output paths for files and directories
    adir_path, bdir_path, outdir, outtab, gffout, tempdir = set_paths(
        adir=args.adir,
        bdir=args.bdir,
        afasta=args.afasta,
        bfasta=args.bfasta,
        outdir=args.outdir,
        outtab=args.outfile,
        gffout=args.gffout,
    )

    # Get all possible alignment pairs between genomes A and B
    pairs: List[tuple] = get_all_pairs(Adir=adir_path, Bdir=bdir_path)

    # Generate chromosome length file for genome A (needed for coverage calculations)
    lenPathA: str = os.path.join(outdir, 'A_gen_lens.txt')
    _chrLensA: Dict[str, int] = chromlens(seqDir=adir_path, outfile=lenPathA)

    # Compose commands for alignment and repeat identification
    cmds: List[str] = xspecies_LZ_cmds(
        lzpath=args.lzpath,
        bdtlsPath=args.bedtools,
        pairs=pairs,
        Adir=adir_path,
        Bdir=bdir_path,
        outtab=outtab,
        outgff=gffout,
        minIdt=args.minIdt,
        minLen=args.minLen,
        minCov=args.minCov,
        AchrmLens=lenPathA,
        reuseTab=args.recycle,
        label=args.label,
        prefix=args.prefix,
    )

    # Run the alignment and processing commands
    logging.info('Running alignments...')
    run_cmd(cmds, verbose=args.verbose, keeptemp=args.keeptemp)

    # Clean up temporary files unless keeptemp flag is set
    if tempdir and os.path.isdir(tempdir) and not args.keeptemp:
        shutil.rmtree(tempdir)

    logging.info('Finished!')
