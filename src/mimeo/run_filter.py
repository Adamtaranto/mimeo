import argparse
import os

import mimeo


def mainArgs():
    parser = argparse.ArgumentParser(
        description='Filter SSR containing sequences from fasta library of repeats.',
        prog='mimeo-filter',
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
    args = parser.parse_args()
    return args


def main():
    # Get cmd line args
    args = mainArgs()
    # Check for required programs.
    tools = [args.TRFpath]
    missing_tools = []
    for tool in tools:
        missing_tools += mimeo.missing_tool(tool)
    if missing_tools:
        print(
            'WARNING: Some tools required by mimeo could not be found: '
            + ', '.join(missing_tools)
        )
        print('You may need to install them to use all features.')
    # Set outfile path
    if not args.outfile:
        outname = os.path.splitext(os.path.basename(args.infile))[0] + '_filtered.fa'
    else:
        outname = args.outfile
    if args.outdir:
        outfile = os.path.join(os.path.abspath(args.outdir), outname)
    else:
        outfile = os.path.join(os.getcwd(), outname)
    # Set abs path to element library
    infile = os.path.abspath(args.infile)
    # Run filter
    mimeo.trfFasta(
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
    print('Finished!')
