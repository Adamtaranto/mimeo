"""
Internal repeat detection through self-alignment.

This module provides command-line functionality for identifying repetitive regions
within a single genome by aligning it to itself. It detects high-identity segments
that exceed specified coverage thresholds, which may represent duplicated regions,
transposons, or other repetitive elements.

The tool can:
1. Split multiFASTA files into individual sequence files if needed
2. Perform self-alignment using LASTZ
3. Identify regions covered by multiple alignments
4. Filter based on identity, length, and coverage thresholds
5. Output results in both tabular and GFF3 formats

This is part of the mimeo package for genome structure analysis.
"""

import argparse
import logging
import os
import shutil
import sys
from typing import List

from ._version import __version__
from .logs import init_logging
from .utils import chromlens, get_all_pairs, missing_tool, run_cmd, set_paths
from .wrappers import self_LZ_cmds


def mainArgs() -> argparse.Namespace:
    """
    Parse command-line arguments for the mimeo-self tool.

    Sets up the argument parser with options for input files/directories,
    output formats, alignment parameters, and filtering thresholds for
    identifying internal repeats within a genome.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Internal repeat finder. Mimeo-self aligns a genome to itself and extracts high-identity segments above an coverage threshold.',
        prog='mimeo-self',
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
        help='Name of directory containing sequences from genome. Write split files here if providing genome as multifasta.',
    )
    parser.add_argument(
        '--afasta', type=str, default=None, help='Genome as multifasta.'
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
        default='mimeo-self_repeats.gff3',
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
        default='Self_Repeat',
        help='Set annotation TYPE field in gff.',
    )
    parser.add_argument(
        '--prefix',
        type=str,
        default='Self_Repeat',
        help='ID prefix for internal repeats.',
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
        default=3,
        help='Minimum depth of aligned segments to report repeat feature.',
    )
    parser.add_argument(
        '--hspthresh',
        type=int,
        default=3000,
        help='Set HSP min score threshold for LASTZ.',
    )
    parser.add_argument(
        '--intraCov',
        type=int,
        default=5,
        help='Minimum depth of aligned segments from same scaffold to report feature. Used if "--strictSelf" mode is selected.',
    )
    parser.add_argument(
        '--strictSelf',
        action='store_true',
        help='If set process same-scaffold alignments separately with option to use higher "--intraCov" threshold. Sometime useful to avoid false repeat calls from staggered alignments over SSRs or short tandem duplication.',
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
    Execute the main self-alignment workflow for repeat detection.

    This function:
    1. Validates required external tools
    2. Sets up output paths
    3. Identifies genome sequence files to self-align
    4. Runs LASTZ alignments (or reuses existing alignments)
    5. Processes alignments to identify repeat regions
    6. Filters regions based on coverage thresholds
    7. Outputs results in tabular and GFF3 formats

    Returns
    -------
    None
        This function does not return any value but performs file I/O and logging.
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
    logging.info('Starting self-alignment workflow.')
    # Log the command line arguments
    logging.debug('Command line arguments: %s', args)

    # Set output paths for files and directories
    adir_path, bdir_path, outdir, outtab, gffout, tempdir = set_paths(
        adir=args.adir,
        afasta=args.afasta,
        outdir=args.outdir,
        outtab=args.outfile,
        gffout=args.gffout,
        suppresBdir=True,  # No B directory needed for self-alignment
    )

    # Get file names to align (self-alignment uses the same files)
    pairs = get_all_pairs(Adir=adir_path, Bdir=bdir_path)

    # Get chromosome lengths for coverage calculations and GFF header
    lenPathA = os.path.join(outdir, 'A_gen_lens.txt')
    _chrLensA = chromlens(seqDir=adir_path, outfile=lenPathA)

    # Generate commands for self-alignment and repeat identification
    cmds: List[str] = self_LZ_cmds(
        lzpath=args.lzpath,
        bdtlsPath=args.bedtools,
        pairs=pairs,
        Adir=adir_path,
        Bdir=bdir_path,
        outtab=outtab,
        outgff=gffout,
        minIdt=args.minIdt,
        minLen=args.minLen,
        hspthresh=args.hspthresh,
        minCov=args.minCov,
        intraCov=args.intraCov,
        splitSelf=args.strictSelf,
        AchrmLens=lenPathA,
        reuseTab=args.recycle,
        label=args.label,
        prefix=args.prefix,
    )

    # Execute the commands for alignment and processing
    logging.info('Running alignments...')
    run_cmd(cmds, verbose=args.verbose, keeptemp=args.keeptemp)

    # Clean up temporary files unless keeptemp flag is set
    if tempdir and os.path.isdir(tempdir) and not args.keeptemp:
        shutil.rmtree(tempdir)

    logging.info('Finished!')
