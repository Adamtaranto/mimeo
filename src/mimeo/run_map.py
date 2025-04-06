"""
Genome alignment and high-identity segment detection.

This module provides command-line functionality for identifying high-identity
segments shared between two genomes using LASTZ alignments. It supports filtering
alignments based on identity percentage, length, and tandem repeat content.

The tool can:
1. Align genome sequences from directories or multiFASTA files
2. Filter alignments by identity percentage and length
3. Filter out regions with excessive tandem repeats using TRF
4. Output results in both tabular and GFF3 formats

This is part of the mimeo package for horizontal gene transfer analysis.
"""

import argparse
import logging
import os
import shutil
import sys
from typing import List

from ._version import __version__
from .logs import init_logging
from .utils import (
    chromlens,
    get_all_pairs,
    missing_tool,
    run_cmd,
    set_paths,
)
from .wrappers import import_Align, map_LZ_cmds, trfFilter, writeGFFlines, writetrf


def mainArgs() -> argparse.Namespace:
    """
    Parse command-line arguments for the mimeo-map tool.

    Sets up the argument parser with options for input files/directories,
    output formats, alignment parameters, and tandem repeat filtering.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Find all high-identity segments shared between genomes.',
        prog='mimeo-map',
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
        default=None,
        help='Name of GFF3 annotation file. If not set, suppress output.',
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
        '--label', type=str, default='BHit', help='Set annotation TYPE field in gff.'
    )
    parser.add_argument(
        '--prefix',
        type=str,
        default='BHit',
        help='ID prefix for B-genome hits annotated in A-genome.',
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
        '--minIdt', type=int, default=60, help='Minimum alignment identity to report.'
    )
    parser.add_argument(
        '--minLen', type=int, default=100, help='Minimum alignment length to report.'
    )
    parser.add_argument(
        '--hspthresh',
        type=int,
        default=3000,
        help='Set HSP min score threshold for LASTZ.',
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
        '--tmaxperiod', type=int, default=50, help='TRF maximum period size to report'
    )
    parser.add_argument(
        '--maxtandem',
        type=float,
        default=None,
        help='Max percentage of an A-genome alignment which may be masked by TRF. If exceeded, alignment will be discarded.',
    )
    parser.add_argument(
        '--writeTRF',
        action='store_true',
        default=False,
        help='If set write TRF filtered alignment file for use with other mimeo modules.',
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
    Execute the main genome mapping workflow.

    This function:
    1. Validates required external tools
    2. Sets up output paths
    3. Identifies genome sequence pairs to align
    4. Runs LASTZ alignments (or reuses existing alignments)
    5. Filters alignments by identity and length
    6. Optionally filters out regions with excessive tandem repeats
    7. Outputs results in tabular and GFF3 formats

    Returns
    -------
    None
        This function does not return a value but performs file I/O and logging.
    """
    # Get cmd line args
    args = mainArgs()

    # Check for required programs
    tools: List[str] = [args.lzpath, args.TRFpath]
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
    logging.info('Starting genome mapping workflow.')
    # Log the command line arguments
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
        runtrf=args.maxtandem,
    )
    # Log the paths
    logging.info(f'Output directory: {outdir}')
    logging.info(f'Output alignment file: {outtab}')
    logging.info(f'Output GFF file: {gffout}')
    logging.info(f'Temporary directory: {tempdir}')

    # Get all possible alignment pairs between genomes A and B

    pairs = get_all_pairs(Adir=adir_path, Bdir=bdir_path)
    # Log count of pairs
    logging.info('Number of pairs to align: %d', len(pairs))

    # Get chromosome lengths for GFF header
    logging.info('Calculating chromosome lengths...')
    chrLens = chromlens(seqDir=adir_path)

    # Log chromosome lengths from chrLens dictionary as for key, value pairs
    chrLens_lines = [f'{chr}: {len}' for chr, len in chrLens]
    logging.info('Chromosome lengths:\n' + '\n'.join(chrLens_lines))

    # Generate new alignments if needed (not in recycle mode or output doesn't exist)
    if not args.recycle or not os.path.isfile(outtab):
        if not pairs:
            logging.error(
                'No files to align. Check --adir and --bdir contain \
                  at least one fasta each.'
            )
            sys.exit(1)

        # Compose alignment commands using LASTZ
        logging.info('Generating alignment commands...')
        cmds: List[str] = map_LZ_cmds(
            lzpath=args.lzpath,
            pairs=pairs,
            minIdt=args.minIdt,
            minLen=args.minLen,
            hspthresh=args.hspthresh,
            outfile=outtab,
            verbose=args.verbose,
        )

        # Run alignment commands
        logging.info('Running alignments...')
        run_cmd(cmds, verbose=args.verbose, keeptemp=args.keeptemp)

    # Import alignment results as dataframe
    logging.info(f'Importing alignments from {outtab}')
    alignments = import_Align(
        infile=outtab, prefix=args.prefix, minLen=args.minLen, minIdt=args.minIdt
    )

    # Optional: Filter alignments based on tandem repeat content
    if args.maxtandem:
        logging.info('Filtering alignments by tandem repeat content...')
        alignments = trfFilter(
            alignDF=alignments,
            tempdir=tempdir,
            prefix=args.prefix,
            adir=adir_path,
            TRFpath=args.TRFpath,
            tmatch=args.tmatch,
            tmismatch=args.tmismatch,
            tdelta=args.tdelta,
            tPM=args.tPM,
            tPI=args.tPI,
            tminscore=args.tminscore,
            tmaxperiod=args.tmaxperiod,
            maxtandem=args.maxtandem,
        )

        # Write TRF-filtered alignments if requested
        if args.writeTRF:
            logging.info(f'Writing TRF-filtered alignments to file: {outtab + ".trf"}')
            writetrf(alignDF=alignments, outtab=outtab)

    # Write results to GFF3 format if output file specified
    if gffout:
        logging.info(f'Writing GFF3 output to {gffout}')
        with open(gffout, 'w') as f:
            for x in writeGFFlines(alnDF=alignments, chrlens=chrLens, ftype=args.label):
                f.write(x)

    # Clean up temporary files unless keeptemp flag is set
    if tempdir and os.path.isdir(tempdir) and not args.keeptemp:
        logging.info(f'Removing temporary files in {tempdir}')
        shutil.rmtree(tempdir)

    logging.info('Finished!')
