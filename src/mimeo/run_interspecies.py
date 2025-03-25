import argparse
import os
import shutil

import mimeo


def mainArgs():
    parser = argparse.ArgumentParser(
        description='Cross-species repeat finder. Mimeo-x searches for features which are abundant in an external reference genome.',
        prog='mimeo-x',
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
    args = parser.parse_args()
    return args


def main():
    # Get cmd line args
    args = mainArgs()
    # Check for required programs.
    tools = [args.lzpath, args.bedtools]
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
    )
    # Get file names to align
    pairs = mimeo.get_all_pairs(Adir=adir_path, Bdir=bdir_path)
    # Get A-genome chromosome lengths for coverage calcs
    lenPathA = os.path.join(outdir, 'A_gen_lens.txt')
    # Make chrm length file for bedtool coverage calc
    chrLensA = mimeo.chromlens(seqDir=adir_path, outfile=lenPathA)
    # Compose alignment commands
    cmds = mimeo.xspecies_LZ_cmds(
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
    # Run alignments
    print('Running alignments...')
    mimeo.run_cmd(cmds, verbose=args.verbose, keeptemp=args.keeptemp)
    if tempdir and os.path.isdir(tempdir) and not args.keeptemp:
        shutil.rmtree(tempdir)
    print('Finished!')
