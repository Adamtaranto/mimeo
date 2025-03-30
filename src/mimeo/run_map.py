import argparse
import os
import shutil
import sys

import mimeo


def mainArgs():
    parser = argparse.ArgumentParser(
        description='Find all high-identity segments shared between genomes.',
        prog='mimeo-map',
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
    args = parser.parse_args()
    return args


def main():
    # Get cmd line args
    args = mainArgs()
    # Check for required programs.
    tools = [args.lzpath, args.TRFpath]
    missing_tools = []
    for tool in tools:
        missing_tools += mimeo.missing_tool(tool)
    if missing_tools:
        print(
            'WARNING: Some tools required by mimeo could not be found: '
            + ', '.join(missing_tools)
        )
        print('You may need to install them to use all features.')
    # Set output paths
    adir_path, bdir_path, outdir, outtab, gffout, tempdir = mimeo.set_paths(
        adir=args.adir,
        bdir=args.bdir,
        afasta=args.afasta,
        bfasta=args.bfasta,
        outdir=args.outdir,
        outtab=args.outfile,
        gffout=args.gffout,
        runtrf=args.maxtandem,
    )
    # Get all B to A alignment pairs
    pairs = mimeo.get_all_pairs(Adir=adir_path, Bdir=bdir_path)
    # Get chrm lens for GFF header
    chrLens = mimeo.chromlens(seqDir=adir_path)
    # Do not realign if outtab exists AND recycle mode is set
    if not args.recycle or not os.path.isfile(outtab):
        if not pairs:
            print(
                'No files to align. Check --adir and --bdir contain \
                  at least one fasta each.'
            )
            sys.exit(1)
        # Compose alignment commands
        cmds = mimeo.map_LZ_cmds(
            lzpath=args.lzpath,
            pairs=pairs,
            minIdt=args.minIdt,
            minLen=args.minLen,
            hspthresh=args.hspthresh,
            outfile=outtab,
            verbose=args.verbose,
        )
        # Run alignments
        mimeo.run_cmd(cmds, verbose=args.verbose, keeptemp=args.keeptemp)
    # Import alignment as df
    alignments = mimeo.import_Align(
        infile=outtab, prefix=args.prefix, minLen=args.minLen, minIdt=args.minIdt
    )
    # Filter alignments if A-genome location >= x% masked by TRF
    if args.maxtandem:
        alignments = mimeo.trfFilter(
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
        if args.writeTRF:
            mimeo.writetrf(alignDF=alignments, outtab=outtab)
    # Write to GFF3
    if gffout:
        with open(gffout, 'w') as f:
            for x in mimeo.writeGFFlines(
                alnDF=alignments, chrlens=chrLens, ftype=args.label
            ):
                f.write(x)
    if tempdir and os.path.isdir(tempdir) and not args.keeptemp:
        shutil.rmtree(tempdir)
    print('Finished!')
