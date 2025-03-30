from collections import Counter
from datetime import datetime, timezone
import glob
import os
from shlex import quote
import shutil
import subprocess
import sys
import tempfile

from Bio import SeqIO
import pandas as pd


def import_pairs(file=None, Adir=None, Bdir=None):
    pairs = []
    with open(file) as f:
        for line in f:
            li = line.strip()
            if not li.startswith('#'):
                A, B = li.split()[:2]
                pairs.append((os.path.join(Adir, A), os.path.join(Bdir, B)))
    return pairs


def get_all_pairs(Adir=None, Bdir=None):
    pairs = []
    if Adir and Bdir:
        for A in glob.glob(os.path.join(Adir, '*')):
            for B in glob.glob(os.path.join(Bdir, '*')):
                pairs.append((A, B))
    elif Adir:
        print('Compose self-genome alignment pairs.')
        for A in glob.glob(os.path.join(Adir, '*')):
            for B in glob.glob(os.path.join(Adir, '*')):
                pairs.append((A, B))
    else:
        print('Need at least one seq directory to compose alignment pairs.')
        sys.exit(1)
    return pairs


def _write_script(cmds, script):
    """Write commands into a bash script"""
    f = open(script, 'w+')
    for cmd in cmds:
        print(cmd, file=f)
    f.close()


def syscall(
    cmd: str | list,
    verbose: bool = False,
    timeout: float | None = None,
    encoding: str = 'utf-8',
    shell: bool = True,
) -> str:
    """Execute a system command with comprehensive error handling.

    This function provides a standardized way to run shell commands while capturing
    output and handling errors consistently. It wraps subprocess.check_output with
    additional error reporting and output management.

    Parameters
    ----------
    cmd : str or list
        Command to execute as string or list of arguments. When shell=True (default),
        this should be a string. When shell=False, this should be a list of strings.
    verbose : bool, optional
        If True, prints the command before execution and its output after completion.
        Default is False.
    timeout : float or None, optional
        Maximum time in seconds for the command to complete before raising a timeout error.
        Default is None (no timeout).
    encoding : str, optional
        Character encoding to use when decoding command output. Default is 'utf-8'.
    shell : bool, optional
        If True, execute command through the system shell. This allows features like
        pipes and redirection but poses security risks with untrusted input. Default is True.

    Returns
    -------
    str
        Decoded command output as string.

    Raises
    ------
    RuntimeError
        If the command returns a non-zero exit code, with details about the failure.
    TimeoutExpired
        If the command execution exceeds the specified timeout.
    """
    # Display the command if in verbose mode
    if verbose:
        print('Running command:', cmd, flush=True)

    try:
        # Execute the command and capture output (including stderr)
        output = subprocess.check_output(
            cmd,
            shell=shell,
            stderr=subprocess.STDOUT,  # Merge stderr into stdout
            timeout=timeout,
        )

        # Decode the binary output using the specified encoding
        decoded_output = output.decode(encoding)

        # Display command output if in verbose mode
        if verbose:
            print(decoded_output, flush=True)

        return decoded_output

    except subprocess.CalledProcessError as error:
        # Handle command execution failures (non-zero exit codes)
        # Print error information to stderr for immediate feedback
        print(
            'The following command failed with exit code',
            error.returncode,
            file=sys.stderr,
            flush=True,
        )
        print(cmd, file=sys.stderr, flush=True)
        print('\nThe output was:\n', file=sys.stderr, flush=True)
        print(error.output.decode(encoding), file=sys.stderr, flush=True)

        # Re-raise as RuntimeError with detailed message while preserving the original error chain
        raise RuntimeError(
            f'Error running command (exit code {error.returncode}): {cmd}'
        ) from error


def run_cmd(cmds, verbose=False, keeptemp=False):
    """Write and excute HMMER script"""
    tmpdir = tempfile.mkdtemp(prefix='tmp.', dir=os.getcwd())
    original_dir = os.getcwd()
    os.chdir(tmpdir)
    script = 'run_jobs.sh'
    _write_script(cmds, script)
    syscall('bash ' + script, verbose=verbose)
    os.chdir(original_dir)
    if not keeptemp:
        shutil.rmtree(tmpdir)


def getTimestring():
    """Return int only string of current datetime with milliseconds."""
    (dt, micro) = datetime.now(timezone.utc).strftime('%Y%m%d%H%M%S.%f').split('.')
    dt = '%s%03d' % (dt, int(micro) // 1000)
    return dt


def splitFasta(infile, outdir, unique=True):
    seen = []
    for rec in SeqIO.parse(infile, 'fasta'):
        if str(rec.id) in seen and unique:
            print('Non-unique name in genome: %s. Quitting.' % str(rec.id))
            sys.exit(1)
        else:
            seen.append(str(rec.id))
        outfile = os.path.join(outdir, rec.id + '.fa')
        with open(outfile, 'w') as handle:
            SeqIO.write(rec, handle, 'fasta')


def isfile(path):
    path = os.path.abspath(path)
    if not os.path.isfile(path):
        print('Input file not found: %s' % path)
        sys.exit(1)
    else:
        return path


def set_paths(
    adir=None,
    bdir=None,
    afasta=None,
    bfasta=None,
    outdir=None,
    outtab=None,
    gffout=None,
    suppresBdir=False,
    runtrf=None,
):
    if not adir:
        # Make temp directory
        tempdir = os.path.join(os.getcwd(), 'temp_' + getTimestring())
        os.makedirs(tempdir)
    elif not bdir and not suppresBdir:
        # Make temp directory
        tempdir = os.path.join(os.getcwd(), 'temp_' + getTimestring())
        os.makedirs(tempdir)
    elif runtrf:
        # Make temp directory
        tempdir = os.path.join(os.getcwd(), 'temp_' + getTimestring())
        os.makedirs(tempdir)
    else:
        tempdir = None
    # Check that A/B genome directories exist
    if adir:
        adir = os.path.abspath(adir)
        if not os.path.isdir(adir):
            print('Creating Adir: %s' % adir)
            os.makedirs(adir)
            if not afasta:
                print('No A-genome fasta file provided. Quitting.')
                sys.exit(1)
    else:
        adir = os.path.join(tempdir, 'A_genome_split')
        os.makedirs(adir)
    if bdir:
        bdir = os.path.abspath(bdir)
        if not os.path.isdir(bdir):
            print('Creating Bdir: %s' % bdir)
            os.makedirs(bdir)
            if not bfasta:
                print('No B-genome fasta file provided. Quitting.')
                sys.exit(1)
    elif not suppresBdir:
        bdir = os.path.join(tempdir, 'B_genome_split')
        os.makedirs(bdir)
    # Split genomes if given
    if afasta:
        if os.path.isfile(afasta):
            splitFasta(afasta, adir)
        else:
            print('A-genome fasta not found at path: %s' % afasta)
    if bfasta:
        if os.path.isfile(bfasta):
            splitFasta(bfasta, bdir)
        elif not suppresBdir:
            print('B-genome fasta not found at path: %s' % bfasta)
    # Set outdir
    if outdir:
        outdir = os.path.abspath(outdir)
        if not os.path.isdir(outdir):
            print('Create output directory: %s' % outdir)
            os.makedirs(outdir)
    else:
        outdir = os.getcwd()
    # Compose path to outfile
    if outtab:
        outtab = os.path.join(outdir, outtab)
        if os.path.isfile(outtab):
            print('Previous alignment found: %s' % outtab)
    if gffout:
        gffout = os.path.join(outdir, gffout)
    # Return paths
    return adir, bdir, outdir, outtab, gffout, tempdir


def checkUniqueID(records):
    """Check that IDs for input elements are unique."""
    seqIDs = [records[x].id for x in range(len(records))]
    IDcounts = Counter(seqIDs)
    duplicates = [k for k, v in IDcounts.items() if v > 1]
    if duplicates:
        print('Input sequence IDs not unique. Quiting.')
        print(duplicates)
        sys.exit(1)
    else:
        pass


def chromlens(seqDir=None, outfile=None):
    """Get chrom lens from fasta dir"""
    records = []
    for f in glob.glob(os.path.join(seqDir, '*')):
        records += list(SeqIO.parse(f, 'fasta'))
    # Check records exist
    if not records:
        print('No sequences found in %s \n Cannot calculate seq lengths.' % seqDir)
        sys.exit(1)
    # Check names are unique
    checkUniqueID(records)
    # Make list
    chrlens = []
    for rec in records:
        chrlens.append((str(rec.id), str(len(rec.seq))))
    # Sort list
    chrlens.sort(key=lambda x: x[0])
    if outfile:
        with open(outfile, 'w') as handle:
            for x, y in chrlens:
                handle.write('\t'.join([x, y]) + '\n')
    return chrlens


def missing_tool(tool_name):
    path = shutil.which(tool_name)
    if path is None:
        return [tool_name]
    else:
        return []


def import_Align(infile=None, prefix=None, minLen=100, minIdt=95):
    """Import LASTZ alignment file to pandas dataframe"""
    hits = []
    with open(infile, 'r') as f:
        for line in f.readlines():
            li = line.strip()
            if not li.startswith('#'):
                li = li.split()
                if int(li[3]) - int(li[2]) >= minLen and float(li[9]) >= minIdt:
                    hits.append(
                        {
                            'tName': li[0],
                            'tStrand': li[1],
                            'tStart': li[2],
                            'tEnd': li[3],
                            'qName': li[4],
                            'qStrand': li[5],
                            'qStart': li[6],
                            'qEnd': li[7],
                            'score': li[8],
                            'pID': li[9],
                            'UID': None,
                        }
                    )
    # Die if no hits
    if not hits:
        print('No alignments found in %s' % infile)
        sys.exit(1)
    # Convert list of dicts into dataframe
    df = pd.DataFrame(hits)
    # Sort hits by Chromosome, location, and strand
    df = df.sort_values(
        ['tName', 'tStart', 'tEnd', 'tStrand'], ascending=[True, True, True, True]
    )
    # Reindex
    df = df.reset_index(drop=True)
    df.index = df.index + 1
    fillLen = len(str(len(df.index)))
    if prefix:
        df['UID'] = str(prefix) + '_' + df.index.astype(str).str.zfill(fillLen)
    else:
        df['UID'] = 'BHit_' + df.index.astype(str).str.zfill(fillLen)
    return df


def trfFilter(
    alignDF=None,
    tempdir=None,
    adir=None,
    prefix=None,
    TRFpath='trf',
    tmatch=2,
    tmismatch=7,
    tdelta=7,
    tPM=80,
    tPI=10,
    tminscore=50,
    tmaxperiod=50,
    maxtandem=40,
):
    alignDF['UID'] = 'TEMPID_' + alignDF.index.astype(str)
    # Reconstitute A-genome fasta from split dir
    seqMaster = {}
    for A in glob.glob(os.path.join(adir, '*')):
        for rec in SeqIO.parse(A, 'fasta'):
            seqMaster[rec.id] = rec
    # Write aligned segments
    AlnFasta = os.path.join(tempdir, 'raw_A_genome_hits.fa')
    with open(AlnFasta, 'w') as handle:
        for _index, row in alignDF.iterrows():
            rec = seqMaster[row['tName']][int(row['tStart']) : int(row['tEnd'])]
            rec.id = row['UID']
            SeqIO.write(rec, handle, 'fasta')
    # Run TRF
    cmds = []
    trf_cmd = ' '.join(
        [
            str(TRFpath),
            quote(AlnFasta),
            str(tmatch),
            str(tmismatch),
            str(tdelta),
            str(tPM),
            str(tPI),
            str(tminscore),
            str(tmaxperiod),
            '-m',
            '-h',
            '-ngs >',
            quote(AlnFasta + '.dat'),
        ]
    )
    # cmds.append(' '.join(["echo 'Run TRF as: ' >&2"]))
    # cmds.append(' '.join(["echo '" + trf_cmd + "' >&2"]))
    cmds.append(trf_cmd)
    maskfile = '.'.join(
        [
            'raw_A_genome_hits.fa',
            str(tmatch),
            str(tmismatch),
            str(tdelta),
            str(tPM),
            str(tPI),
            str(tminscore),
            str(tmaxperiod),
            'mask',
        ]
    )
    cmds.append(' '.join(['mv', maskfile, quote(AlnFasta + '.mask')]))
    alnmasked = AlnFasta + '.mask'
    run_cmd(cmds, verbose=True, keeptemp=False)
    # Make keeplist
    keeplist = []
    for rec in SeqIO.parse(alnmasked, 'fasta'):
        if rec.seq.count('N') / len(rec.seq) * 100 < float(maxtandem):
            keeplist.append(rec.id)
    # Filter align dataframe
    alignments = alignDF.loc[alignDF['UID'].isin(keeplist)].copy()
    # Reindex
    alignments = alignments.sort_values(
        ['tName', 'tStart', 'tEnd', 'tStrand'], ascending=[True, True, True, True]
    )
    alignments = alignments.reset_index(drop=True)
    alignments.index = alignments.index + 1
    fillLen = len(str(len(alignments.index)))
    if prefix:
        alignments['UID'] = (
            str(prefix) + '_' + alignments.index.astype(str).str.zfill(fillLen)
        )
    else:
        alignments['UID'] = 'BHit_' + alignments.index.astype(str).str.zfill(fillLen)
    # Return filtered subset
    return alignments


def trfFasta(
    fasta=None,
    outfile=None,
    TRFpath='trf',
    tmatch=2,
    tmismatch=7,
    tdelta=7,
    tPM=80,
    tPI=10,
    tminscore=50,
    tmaxperiod=50,
    maxtandem=40,
    verbose=True,
    keeptemp=False,
):
    # Run TRF
    cmds = []
    trf_cmd = ' '.join(
        [
            str(TRFpath),
            quote(fasta),
            str(tmatch),
            str(tmismatch),
            str(tdelta),
            str(tPM),
            str(tPI),
            str(tminscore),
            str(tmaxperiod),
            '-m',
            '-h',
            '-ngs >',
            os.path.basename(fasta) + '.dat',
        ]
    )
    cmds.append(' '.join(["echo 'Run TRF as: ' >&2"]))
    # cmds.append(' '.join(["echo '" + trf_cmd + "' >&2"]))
    cmds.append(trf_cmd)
    maskfile = '.'.join(
        [
            os.path.basename(fasta),
            str(tmatch),
            str(tmismatch),
            str(tdelta),
            str(tPM),
            str(tPI),
            str(tminscore),
            str(tmaxperiod),
            'mask',
        ]
    )
    cmds.append(' '.join(['mv', maskfile, quote(fasta + '.mask')]))
    masked = fasta + '.mask'
    run_cmd(cmds, verbose=True, keeptemp=False)
    # Make keeplist
    keeplist = []
    for rec in SeqIO.parse(masked, 'fasta'):
        if rec.seq.count('N') / len(rec.seq) * 100 < float(maxtandem):
            keeplist.append(rec.id)
    # Write seqs that pass filter to new file
    with open(outfile, 'w') as handle:
        for rec in SeqIO.parse(fasta, 'fasta'):
            if rec.id in keeplist:
                SeqIO.write(rec, handle, 'fasta')


def writetrf(alignDF=None, outtab=None):
    # Write alignment dataframe to LZ tab format
    outfile = outtab + '.trf'
    with open(outfile, 'w') as handle:
        handle.write(
            '\t'.join(
                [
                    '#name1',
                    'strand1',
                    'start1',
                    'end1',
                    'name2',
                    'strand2',
                    'start2+',
                    'end2+',
                    'score',
                    'identity' + '\n',
                ]
            )
        )
        for _index, row in alignDF.iterrows():
            handle.write(
                '\t'.join(
                    [
                        row['tName'],
                        row['tStrand'],
                        row['tStart'],
                        row['tEnd'],
                        row['qName'],
                        row['qStrand'],
                        row['qStart'],
                        row['qEnd'],
                        row['score'],
                        row['pID'] + '\n',
                    ]
                )
            )
    return outfile


def writeGFFlines(alnDF=None, chrlens=None, ftype='BHit'):
    yield '##gff-version 3\n'
    if chrlens:
        for name, maxlen in chrlens:
            yield ' '.join(['##sequence-region', str(name), '1', str(maxlen) + '\n'])
    yield '\t'.join(
        [
            '##seqid',
            'source',
            'type',
            'start',
            'end',
            'score',
            'strand',
            'phase',
            'attributes' + '\n',
        ]
    )
    for _index, row in alnDF.iterrows():
        attributes = ';'.join(
            [
                'ID=' + row['UID'],
                'identity=' + str(row['pID']),
                'B_locus='
                + row['qName']
                + '_'
                + row['qStrand']
                + '_'
                + str(row['qStart'])
                + '_'
                + str(row['qEnd']),
            ]
        )
        yield '\t'.join(
            [
                row['tName'],
                'mimeo-map',
                ftype,
                str(row['tStart']),
                str(row['tEnd']),
                str(row['score']),
                row['tStrand'],
                '.',
                attributes + '\n',
            ]
        )


def map_LZ_cmds(
    lzpath: str = 'lastz',
    pairs: list[tuple[str, str]] = None,
    minIdt: float = 95,
    minLen: int = 100,
    hspthresh: int = 3000,
    outfile: str = None,
    verbose: bool = False,
    lastz_format: str = 'general:name1,strand1,start1,end1,length1,name2,strand2,start2+,end2+,length2,score,identity',
    step_size: int = 1,
    strand_mode: str = 'both',
    chain: bool = True,
    gapped: bool = True,
) -> list[str]:
    """Generate commands to map high identity segments between genome files using LASTZ.

    This function creates a list of shell commands that identify regions with high sequence
    identity between pairs of genome files. These regions can represent candidate horizontal
    gene transfer (HGT) regions. The commands run LASTZ alignments and filter/format the results.

    Parameters
    ----------
    lzpath : str, optional
        Path to LASTZ executable, default is 'lastz'
    pairs : list of (str, str) tuples, optional
        List of file path pairs to align. Each tuple contains (target_path, query_path)
    minIdt : float, optional
        Minimum percent identity threshold for reported alignments, default is 95
    minLen : int, optional
        Minimum alignment length to report (bp), default is 100
    hspthresh : int, optional
        High-scoring segment pair threshold for LASTZ sensitivity, default is 3000
    outfile : str, optional
        Path where filtered alignments will be written
    verbose : bool, optional
        Enable verbose output during alignment, default is False
    lastz_format : str, optional
        Format string for LASTZ output fields
    step_size : int, optional
        Step size parameter for LASTZ seed generation, default is 1
    strand_mode : str, optional
        Strand to align: 'plus', 'minus', or 'both', default is 'both'
    chain : bool, optional
        Enable chaining of alignments, default is True
    gapped : bool, optional
        Enable gapped extensions, default is True

    Returns
    -------
    list of str
        Shell commands to run LASTZ alignments and process the results

    Notes
    -----
    The commands perform the following operations for each genome pair:
    1. Run LASTZ alignment with specified parameters
    2. Format output and remove percentage symbols
    3. Filter alignments by minimum identity and length thresholds
    4. Sort results by chromosome, start, and end positions
    """
    # Validate required inputs
    if pairs is None or len(pairs) == 0:
        raise ValueError('No sequence pairs provided for alignment')

    if outfile is None:
        raise ValueError('Output file path is required')

    # Determine platform for sed command compatibility
    import platform

    system = platform.system().lower()
    if system == 'darwin':  # macOS
        sed_command = "sed -i '' -e 's/%//g'"
    else:  # Linux and others
        sed_command = "sed -i 's/%//g'"

    # Set LASTZ verbosity level
    verbosity_level = 1 if verbose else 0

    # Build LASTZ core options
    lastz_options = [
        '--entropy',
        f'--format={lastz_format}',
        '--markend',
        '--gfextend',
    ]

    # Add conditional options based on parameters
    if chain:
        lastz_options.append('--chain')
    if gapped:
        lastz_options.append('--gapped')

    # Add remaining options
    lastz_options.extend(
        [f'--step={step_size}', f'--strand={strand_mode}', f'--hspthresh={hspthresh}']
    )

    # Initialize command list
    cmds = []

    # Write header to output file
    header_command = (
        "echo $'#name1\\tstrand1\\tstart1\\tend1\\tname2\\tstrand2\\tstart2+\\tend2+\\tscore\\tidentity' >"
        f' {quote(outfile)}'
    )
    cmds.append(header_command)

    # Process each genome pair
    # A = target genome
    # B = query genome
    # Loop through each pair of genome files
    # and generate LASTZ commands
    for target_file, query_file in pairs:
        # Extract file information for naming temporary files
        target_name = os.path.splitext(os.path.basename(target_file))[0]
        query_name = os.path.splitext(os.path.basename(query_file))[0]
        temp_outfile = f'temp_{query_name}_onto_{target_name}.tab'

        # Compose LASTZ alignment command
        lastz_cmd = [
            lzpath,
            quote(target_file),
            quote(query_file),
            ' '.join(lastz_options),
            f'--output={temp_outfile}',
            f'--verbosity={verbosity_level}',
        ]
        cmds.append(' '.join(lastz_cmd))

        # Remove percentage symbols from identity values
        cmds.append(f'{sed_command} {temp_outfile}')

        # Filter, format and sort alignment results
        # 1. Extract non-header lines
        # 2. Filter by minimum length
        # 3. Filter by minimum identity
        # 4. Format into tab-delimited columns
        # 5. Sort by chromosome, start, and end positions
        # 6. Append to final output file
        filter_cmd = [
            "awk '!/^#/ { print; }'",
            temp_outfile,
            f'| awk -v minLen={minLen}',
            f"'0+$5 >= minLen {{print ;}}'| awk -v OFS='\\t' -v minIdt={minIdt}",
            "'0+$13 >= minIdt {print $1,$2,$3,$4,$6,$7,$8,$9,$11,$13;}'",
            "| sed 's/ //g'",
            '| sort -k 1,1 -k 3n,4n >>',
            quote(outfile),
        ]
        cmds.append(' '.join(filter_cmd))

        # Remove temporary file
        cmds.append(f'rm {temp_outfile}')

    return cmds


def xspecies_LZ_cmds(
    lzpath='lastz',
    bdtlsPath='bedtools',
    Adir=None,
    Bdir=None,
    pairs=None,
    outtab=None,
    outgff=None,
    minIdt=60,
    minLen=100,
    hspthresh=3000,
    minCov=5,
    AchrmLens=None,
    reuseTab=False,
    label='B_repeats',
    prefix=None,
    verbose=False,
):
    """
    Screen genome A (target) for features which are high copy in
    genome B (query).
    """
    if verbose:
        verb = 1
    else:
        verb = 0
    cmds = []
    if not reuseTab or not os.path.isfile(outtab):
        cmds.append(
            ' '.join(
                [
                    "echo $'#name1\\tstrand1\\tstart1\\tend1\\tname2\\tstrand2\\tstart2+\\tend2+\\tscore\\tidentity' >",
                    quote(outtab),
                ]
            )
        )
        # loop writing B to A mapping cmds
        for A, B in pairs:
            t_file = A
            t_name = os.path.splitext(os.path.basename(A))[0]
            q_file = B
            q_name = os.path.splitext(os.path.basename(B))[0]
            temp_outfile = '_'.join(['temp', q_name, 'onto', t_name, '.tab'])
            # Compose LASTZ command
            cmds.append(
                ' '.join(
                    [
                        lzpath,
                        quote(t_file),
                        quote(q_file),
                        '--entropy --format=general:name1,strand1,start1,end1,length1,name2,strand2,start2+,end2+,length2,score,identity --markend --gfextend --chain --gapped --step=1 --strand=both --hspthresh='
                        + str(hspthresh),
                        '--output=' + temp_outfile,
                        '--verbosity=' + str(verb),
                    ]
                )
            )
            # Scrub % symbols
            cmds.append(' '.join(["sed -i '' -e 's/%//g'", temp_outfile]))
            # Filter Inter_Chrome targets to min len $minLen [100], min identity $minIdt [90]
            # New Header = name1,strand1,start1,end1,name2,strand2,start2+,end2+,score,identity
            # Sort filtered file by chrom, start, stop
            cmds.append(
                ' '.join(
                    [
                        "awk '!/^#/ { print; }'",
                        temp_outfile,
                        '| awk -v minLen=' + str(minLen),
                        "'0+$5 >= minLen {print ;}' | awk -v OFS='\\t' -v minIdt="
                        + str(minIdt),
                        "'0+$13 >= minIdt {print $1,$2,$3,$4,$6,$7,$8,$9,$11,$13;}' | sed 's/ //g' | sort -k 1,1 -k 3n,4n >>",
                        quote(outtab),
                    ]
                )
            )
            cmds.append(' '.join(['rm', temp_outfile]))
    # Set temp output files
    temp_bed = 'temp.bed'
    temp_bed_sorted = 'temp_sorted.bed'
    # Port filtered hits to bed format: "name1,start1,end1"
    # cmds.append(' '.join(["echo 'Port filtered hits to bed format: name1,start1,end1' >&2"]))
    cmds.append(
        ' '.join(["awk -v OFS='\\t' '!/^#/ {print $1,$3,$4;}'", outtab, '>', temp_bed])
    )
    # Remove spaces
    cmds.append(' '.join(["sed -i '' -e 's/ //g'", temp_bed]))
    # Sort filtered bed file by chrom, start, stop
    cmds.append(' '.join(['sort -k 1,1 -k 2n,3n', temp_bed, '>', temp_bed_sorted]))
    # Generate non-zero coverage scores for target genome regions, filter for min coverage of x
    cmds.append(
        ' '.join(
            [
                "echo 'Generate non-zero coverage scores for target genome regions, filter for min coverage of x' >&2"
            ]
        )
    )
    cmds.append(
        ' '.join(
            [
                bdtlsPath,
                'genomecov -bg -i',
                temp_bed_sorted,
                '-g',
                AchrmLens,
                "| awk -v OFS='\\t' -v cov=" + str(minCov),
                "'0+$4 >= cov {print ;}' >",
                temp_bed,
            ]
        )
    )
    # Re-sort
    cmds.append(' '.join(['sort -k 1,1 -k 2n,3n', temp_bed, '>', temp_bed_sorted]))
    # Merge end-to-end annotations
    cmds.append(' '.join([bdtlsPath, 'merge -i', temp_bed_sorted, '>', temp_bed]))
    # Write GFF header
    cmds.append(' '.join(["echo 'Writing GFF' >&2"]))
    cmds.append(
        ' '.join(
            [
                "echo $'##gff-version 3\\n#seqid\\tsource\\ttype\\tstart\\tend\\tscore\\tstrand\\tphase\\tattributes' >",
                quote(outgff),
            ]
        )
    )
    # Filter orphan fragments < xx && Create GFF3 file
    cmds.append(
        ' '.join(
            [
                'awk -v minLen=' + str(minLen),
                "'{ if($3 - $2 >= minLen) print ;}'",
                temp_bed,
                '| awk -v OFS=\'\\t\' \'BEGIN{i=0}{i++;}{j= sprintf("%05d", i)}{print $1,"mimeo","'
                + str(label)
                + '",$2,$3,".","+",".","ID='
                + str(prefix)
                + '_"j ;}\' >>',
                quote(outgff),
            ]
        )
    )
    return cmds


def self_LZ_cmds(
    lzpath='lastz',
    bdtlsPath='bedtools',
    splitSelf=False,
    Adir=None,
    Bdir=None,
    pairs=None,
    outtab=None,
    outgff=None,
    minIdt=60,
    minLen=100,
    hspthresh=3000,
    minCov=3,
    intraCov=5,
    AchrmLens=None,
    reuseTab=False,
    label='Self_repeats',
    prefix=None,
    verbose=False,
):
    """
    Align genome to itself. Opt: Process intra-chrom alignments separately to
    avoid local SSRs.
    """
    if verbose:
        verb = 1
    else:
        verb = 0
    cmds = []
    if splitSelf:
        outtab_intra = outtab + '_intra.tab'
    if not reuseTab or not os.path.isfile(outtab):
        # Write alignment out file headers
        cmds.append(
            ' '.join(
                [
                    "echo $'#name1\tstrand1\tstart1\tend1\tname2\tstrand2\tstart2+\tend2+\tscore\tidentity' >",
                    quote(outtab),
                ]
            )
        )
        if splitSelf:
            outtab_intra = outtab + '_intra.tab'
            cmds.append(
                ' '.join(
                    [
                        "echo $'#name1\tstrand1\tstart1\tend1\tname2\tstrand2\tstart2+\tend2+\tscore\tidentity' >",
                        quote(outtab_intra),
                    ]
                )
            )
        # loop writing B to A mapping cmds
        for A, B in pairs:
            if A != B or not splitSelf:
                t_file = A
                t_name = os.path.splitext(os.path.basename(A))[0]
                q_file = B
                q_name = os.path.splitext(os.path.basename(B))[0]
                temp_outfile = '_'.join(['temp', q_name, 'onto', t_name, '.tab'])
                # Compose LASTZ command
                cmds.append(
                    ' '.join(
                        [
                            lzpath,
                            quote(t_file),
                            quote(q_file),
                            '--entropy --format=general:name1,strand1,start1,end1,length1,name2,strand2,start2+,end2+,length2,score,identity --markend --gfextend --chain --gapped --step=1 --strand=both --hspthresh='
                            + str(hspthresh),
                            '--output=' + temp_outfile,
                            '--verbosity=' + str(verb),
                        ]
                    )
                )
                # Scrub % symbols
                cmds.append(' '.join(["sed -i '' -e 's/%//g'", temp_outfile]))
                # Filter Inter_Chrome targets to min len $minLen [100], min identity $minIdt [90]
                # New Header = name1,strand1,start1,end1,name2,strand2,start2+,end2+,score,identity
                # Sort filtered file by chrom, start, stop
                cmds.append(
                    ' '.join(
                        [
                            "awk '!/^#/ { print; }'",
                            temp_outfile,
                            '| awk -v minLen=' + str(minLen),
                            "'0+$5 >= minLen {print ;}' | awk -v OFS='\\t' -v minIdt="
                            + str(minIdt),
                            "'0+$13 >= minIdt {print $1,$2,$3,$4,$6,$7,$8,$9,$11,$13;}' | sed 's/ //g' | sort -k 1,1 -k 3n,4n >>",
                            quote(outtab),
                        ]
                    )
                )
                cmds.append(' '.join(['rm', temp_outfile]))
            else:
                # Align to self with diff settings
                t_file = A
                t_name = os.path.splitext(os.path.basename(A))[0]
                q_file = B
                q_name = os.path.splitext(os.path.basename(B))[0]
                temp_outfile = '_'.join(['temp', q_name, 'onto', t_name, '.tab'])
                # Compose LASTZ command
                cmds.append(
                    ' '.join(
                        [
                            lzpath,
                            quote(t_file),
                            quote(q_file),
                            '--entropy --format=general:name1,strand1,start1,end1,length1,name2,strand2,start2+,end2+,length2,score,identity --markend --gfextend --chain --gapped --step=1 --strand=both --hspthresh='
                            + str(hspthresh),
                            '--output=' + temp_outfile,
                            '--verbosity=' + str(verb),
                        ]
                    )
                )
                # Scrub % symbols
                cmds.append(' '.join(["sed -i '' -e 's/%//g'", temp_outfile]))
                # Filter Inter_Chrome targets to min len $minLen [100], min identity $minIdt [90]
                # New Header = name1,strand1,start1,end1,name2,strand2,start2+,end2+,score,identity
                # Sort filtered file by chrom, start, stop
                cmds.append(
                    ' '.join(
                        [
                            "awk '!/^#/ { print; }'",
                            temp_outfile,
                            '| awk -v minLen=' + str(minLen),
                            "'0+$5 >= minLen {print ;}' | awk -v OFS='\\t' -v minIdt="
                            + str(minIdt),
                            "'0+$13 >= minIdt {print $1,$2,$3,$4,$6,$7,$8,$9,$11,$13;}' | sed 's/ //g' | sort -k 1,1 -k 3n,4n >>",
                            quote(outtab_intra),
                        ]
                    )
                )
                cmds.append(' '.join(['rm', temp_outfile]))
    # Coverage filtering for BETWEEN chromosome hits (or all if not in selfSplit mode)
    cmds.append(
        ' '.join(
            [
                "echo 'Coverage filtering for BETWEEN chromosome hits (or all if not in selfSplit mode)' >&2"
            ]
        )
    )
    # Set temp output files
    temp_bed = 'temp.bed'
    temp_bed_sorted = 'temp_sorted.bed'
    # Port filtered hits to bed format: "name1,start1,end1"
    cmds.append(
        ' '.join(["awk -v OFS='\\t' '!/^#/ {print $1,$3,$4;}'", outtab, '>', temp_bed])
    )
    # Remove spaces
    cmds.append(' '.join(["sed -i '' -e 's/ //g'", temp_bed]))
    # Sort filtered bed file by chrom, start, stop
    cmds.append(' '.join(['sort -k 1,1 -k 2n,3n', temp_bed, '>', temp_bed_sorted]))
    # Generate non-zero coverage scores for target genome regions,
    # filter for min coverage of x
    cmds.append(
        ' '.join(
            [
                bdtlsPath,
                'genomecov -bg -i',
                temp_bed_sorted,
                '-g',
                AchrmLens,
                "| awk -v OFS='\\t' -v cov=" + str(minCov),
                "'0+$4 >= cov {print ;}' >",
                temp_bed,
            ]
        )
    )
    # Re-sort
    cmds.append(' '.join(['sort -k 1,1 -k 2n,3n', temp_bed, '>', temp_bed_sorted]))
    # Merge end-to-end annotations
    cmds.append(' '.join([bdtlsPath, 'merge -i', temp_bed_sorted, '>', temp_bed]))
    # Initialise coverage GFF file
    cmds.append(
        ' '.join(
            [
                "echo $'##gff-version 3\n#seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes' >",
                quote(outgff),
            ]
        )
    )
    # Filter orphan fragments < xx && Create GFF3 file
    cmds.append(
        ' '.join(
            [
                'awk -v minLen=' + str(minLen),
                "'{ if($3 - $2 >= minLen) print ;}'",
                temp_bed,
                '| awk -v OFS=\'\\t\' \'BEGIN{i=0}{i++;}{j= sprintf("%05d", i)}{print $1,"mimeo-self","'
                + str(label)
                + '",$2,$3,".","+",".","ID='
                + str(prefix)
                + '_"j ;}\' >>',
                quote(outgff),
            ]
        )
    )
    # Coverage filtering for WITHIN chromosome hits
    if splitSelf:
        if reuseTab and not os.path.isfile(outtab_intra) and os.path.isfile(outtab):
            print(
                "Warning: Could not find intra-chrom results file: %s \nRe-run in '--strictSelf' mode if required."
                % outtab_intra
            )
        else:
            cmds.append(
                ' '.join(
                    [
                        "echo 'Applying separate coverage filtering for WITHIN chromosome hits' >&2"
                    ]
                )
            )
            temp_bed_intra = 'temp.bed'
            temp_bed_sorted_intra = 'temp_sorted.bed'
            # Port filtered hits to bed format: "name1,start1,end1"
            cmds.append(
                ' '.join(
                    [
                        "awk -v OFS='\\t' '!/^#/ {print $1,$3,$4;}'",
                        outtab_intra,
                        '>',
                        temp_bed_intra,
                    ]
                )
            )
            # Remove spaces
            cmds.append(' '.join(["sed -i '' -e 's/ //g'", temp_bed_intra]))
            # Sort filtered bed file by chrom, start, stop
            cmds.append(
                ' '.join(
                    ['sort -k 1,1 -k 2n,3n', temp_bed_intra, '>', temp_bed_sorted_intra]
                )
            )
            # Generate non-zero coverage scores for target genome regions,
            # filter for min coverage of x
            # Note: Recycling $temp_bed_intra
            cmds.append(
                ' '.join(
                    [
                        bdtlsPath,
                        'genomecov -bg -i',
                        temp_bed_sorted_intra,
                        '-g',
                        AchrmLens,
                        "| awk -v OFS='\\t' -v cov=" + str(intraCov),
                        "'0+$4 >= cov {print ;}' >",
                        temp_bed_intra,
                    ]
                )
            )
            # Re-sort
            cmds.append(
                ' '.join(
                    ['sort -k 1,1 -k 2n,3n', temp_bed_intra, '>', temp_bed_sorted_intra]
                )
            )
            # Merge end-to-end annotations
            cmds.append(
                ' '.join(
                    [bdtlsPath, 'merge -i', temp_bed_sorted_intra, '>', temp_bed_intra]
                )
            )
            # Filter orphan fragments < xx && Create GFF3 file
            cmds.append(
                ' '.join(
                    [
                        'awk -v minLen=' + str(minLen),
                        "'{ if($3 - $2 >= minLen) print ;}'",
                        temp_bed_intra,
                        '| awk -v OFS=\'\\t\' \'BEGIN{i=0}{i++;}{j= sprintf("%05d", i)}{print $1,"mimeo-self","'
                        + str(label)
                        + '_intra'
                        + '",$2,$3,".","+",".","ID='
                        + str(prefix)
                        + '_"j ;}\' >>',
                        quote(outgff),
                    ]
                )
            )
    # Return cmds list
    return cmds
