"""
Simple Sequence Repeat (SSR) filtering for sequence libraries.

This module provides command-line functionality for identifying and removing
sequences that contain an excessive proportion of tandem repeats from FASTA files.
It is particularly useful for cleaning repeat libraries of SSR-contaminated entries
before downstream analyses.

The tool can:
1. Process FASTA files containing genomic sequences or repeat libraries
2. Identify tandem repeats using the Tandem Repeats Finder (TRF) tool
3. Filter out sequences where the repeat content exceeds a specified threshold
4. Output a cleaned FASTA file containing only non-SSR-dominated sequences

This is part of the mimeo package for genomic sequence analysis.
"""

import argparse
import logging
import os
import sys
from typing import List

from ._version import __version__
from .logs import init_logging
from .utils import missing_tool
from .wrappers import trfFasta


def mainArgs() -> argparse.Namespace:
    """
    Parse command-line arguments for the mimeo-filter tool.

    Sets up the argument parser with options for input files, output locations,
    and Tandem Repeats Finder (TRF) parameters for identifying and filtering
    sequences with excessive simple sequence repeats.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Filter SSR containing sequences from fasta library of repeats.',
        prog='mimeo-filter',
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
        '--infile',
        type=str,
        required=True,
        help='Name of directory containing sequences from A genome.',
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
        '--outfile', type=str, default=None, help='Name of alignment result file.'
    )
    parser.add_argument(
        '--keeptemp',
        action='store_true',
        default=False,
        help='If set do not remove temp files.',
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        default=False,
        help='If set report LASTZ progress.',
    )
    # TRF filtering
    parser.add_argument(
        '--TRFpath',
        type=str,
        default='trf',
        help='Custom path to TRF executable if not in $PATH.',
    )
    parser.add_argument('--tmatch', type=int, default=2, help='TRF matching weight')
    parser.add_argument(
        '--tmismatch', type=int, default=7, help='TRF mismatching penalty'
    )
    parser.add_argument('--tdelta', type=int, default=7, help='TRF indel penalty')
    parser.add_argument('--tPM', type=int, default=80, help='TRF match probability')
    parser.add_argument('--tPI', type=int, default=10, help='TRF indel probability')
    parser.add_argument(
        '--tminscore',
        type=int,
        default=50,
        help='TRF minimum alignment score to report',
    )
    parser.add_argument(
        '--tmaxperiod',
        type=int,
        default=50,
        help='TRF maximum period size to report. Note: Setting this score too high may exclude some LTR retrotransposons. Optimal len to exclude only SSRs is 10-50bp.',
    )
    parser.add_argument(
        '--maxtandem',
        type=float,
        default=40,
        help='Max percentage of a sequence which may be masked by TRF. If exceeded, element will be discarded.',
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
    Execute the main sequence filtering workflow.

    This function:
    1. Validates the required external TRF tool
    2. Sets up output paths
    3. Processes the input FASTA file using TRF to identify tandem repeats
    4. Filters out sequences where tandem repeat content exceeds the threshold
    5. Writes the filtered sequences to the output file

    Returns
    -------
    None
        Function does not return any value. It writes the filtered sequences
        to the specified output file.
    """
    # Get command line arguments
    args = mainArgs()

    # Check for required external programs
    tools: List[str] = [args.TRFpath]
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
    logging.info('Starting SSR filtering process.')
    logging.debug('Command line arguments: %s', args)

    # Determine output file path
    if not args.outfile:
        # If no output filename specified, create one based on input filename
        outname = os.path.splitext(os.path.basename(args.infile))[0] + '_filtered.fa'
    else:
        outname = args.outfile

    # Set full output path, either in specified directory or current working directory
    if args.outdir:
        outfile = os.path.join(os.path.abspath(args.outdir), outname)
    else:
        outfile = os.path.join(os.getcwd(), outname)

    # Convert input file path to absolute path
    infile = os.path.abspath(args.infile)

    # Run the TRF-based filtering process
    trfFasta(
        fasta=infile,
        outfile=outfile,
        TRFpath=args.TRFpath,
        tmatch=args.tmatch,
        tmismatch=args.tmismatch,
        tdelta=args.tdelta,
        tPM=args.tPM,
        tPI=args.tPI,
        tminscore=args.tminscore,
        tmaxperiod=args.tmaxperiod,
        maxtandem=args.maxtandem,
        verbose=args.verbose,
        keeptemp=args.keeptemp,
    )

    logging.info('Finished!')
