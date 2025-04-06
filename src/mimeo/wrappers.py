"""
Wrapper functions for genomic alignment and repeat identification tools.

This module provides wrapper functions for external tools like LASTZ and TRF (Tandem
Repeats Finder), as well as utilities for processing, filtering, and formatting
genomic alignment data. It includes functions for:

- Reading and writing alignment files in various formats
- Generating shell commands for genomic alignments
- Filtering alignments based on identity and length thresholds
- Processing tandem repeat data
- Converting between alignment formats (table, BED, GFF3)
- Identifying repetitive elements through self-alignments and cross-species comparisons

These functions support the core analysis pipeline in the mimeo package for
horizontal gene transfer and repetitive element detection.
"""

import glob
import logging
import os
import platform
from shlex import quote
import sys
from typing import List, Tuple

from Bio import SeqIO
import pandas as pd

from .utils import run_cmd


def import_Align(
    infile: str = None, prefix: str = None, minLen: int = 100, minIdt: float = 95
) -> pd.DataFrame:
    """
    Import LASTZ alignment file to pandas dataframe.

    Reads a LASTZ alignment file in tabular format, filters alignments based on
    minimum length and identity thresholds, and organizes the data into a pandas
    DataFrame with unique identifiers.

    Parameters
    ----------
    infile : str
        Path to LASTZ alignment file in tabular format.
    prefix : str, optional
        Prefix to use when generating unique IDs for alignments.
    minLen : int, optional
        Minimum alignment length to include, default is 100.
    minIdt : float, optional
        Minimum percent identity to include, default is 95.

    Returns
    -------
    pd.DataFrame
        DataFrame containing filtered alignment data with columns:
        tName, tStrand, tStart, tEnd, qName, qStrand, qStart, qEnd, score, pID, UID.

    Raises
    ------
    SystemExit
        If no alignments are found or none pass the filtering criteria.
    """
    # Initialize empty list to store alignment hits
    hits = []

    # Read and parse alignment file
    with open(infile, 'r') as f:
        for line in f.readlines():
            li = line.strip()
            # Skip header lines
            if not li.startswith('#'):
                li = li.split()
                # Apply length and identity filters
                if int(li[3]) - int(li[2]) >= minLen and float(li[9]) >= minIdt:
                    hits.append(
                        {
                            'tName': li[0],  # Target name
                            'tStrand': li[1],  # Target strand
                            'tStart': li[2],  # Target start position
                            'tEnd': li[3],  # Target end position
                            'qName': li[4],  # Query name
                            'qStrand': li[5],  # Query strand
                            'qStart': li[6],  # Query start position
                            'qEnd': li[7],  # Query end position
                            'score': li[8],  # Alignment score
                            'pID': li[9],  # Percent identity
                            'UID': None,  # Will be filled with unique ID
                        }
                    )

    # Exit if no hits were found
    if not hits:
        logging.warning('No alignments found in %s' % infile)
        sys.exit(1)

    # Convert list of dicts into dataframe
    df = pd.DataFrame(hits)

    # Sort hits by Chromosome, location, and strand
    df = df.sort_values(
        ['tName', 'tStart', 'tEnd', 'tStrand'], ascending=[True, True, True, True]
    )

    # Reindex and create unique IDs
    df = df.reset_index(drop=True)
    df.index = df.index + 1
    fillLen = len(str(len(df.index)))

    # Create unique identifiers with provided prefix or default
    if prefix:
        df['UID'] = str(prefix) + '_' + df.index.astype(str).str.zfill(fillLen)
    else:
        df['UID'] = 'BHit_' + df.index.astype(str).str.zfill(fillLen)

    return df


def trfFilter(
    alignDF: pd.DataFrame = None,
    tempdir: str = None,
    adir: str = None,
    prefix: str = None,
    TRFpath: str = 'trf',
    tmatch: int = 2,
    tmismatch: int = 7,
    tdelta: int = 7,
    tPM: int = 80,
    tPI: int = 10,
    tminscore: int = 50,
    tmaxperiod: int = 50,
    maxtandem: float = 40,
) -> pd.DataFrame:
    """
    Filter alignments based on tandem repeat content using TRF.

    Extracts sequences for each alignment from genome A, runs Tandem Repeats Finder (TRF)
    on them, and filters out alignments where the tandem repeat content exceeds a specified
    threshold (maxtandem).

    Parameters
    ----------
    alignDF : pd.DataFrame
        DataFrame containing alignment data with columns as created by import_Align.
    tempdir : str
        Path to directory for temporary files.
    adir : str
        Path to directory containing A genome sequence files.
    prefix : str, optional
        Prefix to use when regenerating unique IDs for filtered alignments.
    TRFpath : str, optional
        Path to TRF executable, default is 'trf'.
    tmatch : int, optional
        TRF matching weight parameter, default is 2.
    tmismatch : int, optional
        TRF mismatch penalty parameter, default is 7.
    tdelta : int, optional
        TRF indel penalty parameter, default is 7.
    tPM : int, optional
        TRF match probability parameter, default is 80.
    tPI : int, optional
        TRF indel probability parameter, default is 10.
    tminscore : int, optional
        TRF minimum score to report parameter, default is 50.
    tmaxperiod : int, optional
        TRF maximum period size parameter, default is 50.
    maxtandem : float, optional
        Maximum percentage of sequence that can be tandem repeats, default is 40.

    Returns
    -------
    pd.DataFrame
        DataFrame containing filtered alignments that pass the tandem repeat threshold.
    """
    # Assign temporary unique IDs for processing
    alignDF['UID'] = 'TEMPID_' + alignDF.index.astype(str)

    # Reconstitute A-genome fasta from split dir
    seqMaster = {}
    for A in glob.glob(os.path.join(adir, '*')):
        for rec in SeqIO.parse(A, 'fasta'):
            seqMaster[rec.id] = rec

    # Write aligned segments to a FASTA file for TRF processing
    AlnFasta = os.path.join(tempdir, 'raw_A_genome_hits.fa')
    with open(AlnFasta, 'w') as handle:
        for _index, row in alignDF.iterrows():
            # Extract sequence for this alignment from master sequence
            rec = seqMaster[row['tName']][int(row['tStart']) : int(row['tEnd'])]
            rec.id = row['UID']
            SeqIO.write(rec, handle, 'fasta')

    # Prepare TRF command
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
            '-m',  # Masked sequence output
            '-h',  # Suppress HTML output
            '-ngs >',  # NGS mode with stdout redirect
            quote(AlnFasta + '.dat'),
        ]
    )
    cmds.append(trf_cmd)

    # Move mask file to expected location
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

    # Run TRF and process the output
    run_cmd(cmds, verbose=True, keeptemp=False)

    # Create list of sequence IDs that pass the tandem repeat filter
    keeplist = []
    for rec in SeqIO.parse(alnmasked, 'fasta'):
        # Calculate percentage of sequence masked as tandem repeats (N characters)
        if rec.seq.count('N') / len(rec.seq) * 100 < float(maxtandem):
            keeplist.append(rec.id)

    # Filter alignment dataframe to keep only sequences that pass the filter
    alignments = alignDF.loc[alignDF['UID'].isin(keeplist)].copy()

    # Sort and reindex the filtered alignments
    alignments = alignments.sort_values(
        ['tName', 'tStart', 'tEnd', 'tStrand'], ascending=[True, True, True, True]
    )
    alignments = alignments.reset_index(drop=True)
    alignments.index = alignments.index + 1

    # Generate new unique IDs with proper zero-padding
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
    fasta: str = None,
    outfile: str = None,
    TRFpath: str = 'trf',
    tmatch: int = 2,
    tmismatch: int = 7,
    tdelta: int = 7,
    tPM: int = 80,
    tPI: int = 10,
    tminscore: int = 50,
    tmaxperiod: int = 50,
    maxtandem: float = 40,
    verbose: bool = True,
    keeptemp: bool = False,
) -> None:
    """
    Filter sequences in a FASTA file based on tandem repeat content.

    Processes a FASTA file using Tandem Repeats Finder (TRF), then filters out
    sequences where the tandem repeat content exceeds a specified threshold.

    Parameters
    ----------
    fasta : str
        Path to input FASTA file.
    outfile : str
        Path to output filtered FASTA file.
    TRFpath : str, optional
        Path to TRF executable, default is 'trf'.
    tmatch : int, optional
        TRF matching weight parameter, default is 2.
    tmismatch : int, optional
        TRF mismatch penalty parameter, default is 7.
    tdelta : int, optional
        TRF indel penalty parameter, default is 7.
    tPM : int, optional
        TRF match probability parameter, default is 80.
    tPI : int, optional
        TRF indel probability parameter, default is 10.
    tminscore : int, optional
        TRF minimum score to report parameter, default is 50.
    tmaxperiod : int, optional
        TRF maximum period size parameter, default is 50.
    maxtandem : float, optional
        Maximum percentage of sequence that can be tandem repeats, default is 40.
    verbose : bool, optional
        If True, display verbose output, default is True.
    keeptemp : bool, optional
        If True, keep temporary files, default is False.

    Returns
    -------
    None
        Writes filtered sequences to output file.
    """
    # Run TRF
    cmds = []

    # Construct TRF command
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
            '-m',  # Masked sequence output
            '-h',  # Suppress HTML output
            '-ngs >',  # NGS mode with stdout redirect
            os.path.basename(fasta) + '.dat',
        ]
    )

    # Add commands to echo TRF command and run it
    cmds.append(' '.join(["echo 'Run TRF as: ' >&2"]))
    cmds.append(trf_cmd)

    # Move mask file to expected location
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

    # Execute commands
    run_cmd(cmds, verbose=True, keeptemp=False)

    # Create list of sequence IDs that pass the tandem repeat filter
    keeplist = []
    for rec in SeqIO.parse(masked, 'fasta'):
        # Calculate percentage of sequence masked as tandem repeats (N characters)
        if rec.seq.count('N') / len(rec.seq) * 100 < float(maxtandem):
            keeplist.append(rec.id)

    # Write sequences that pass filter to new file
    with open(outfile, 'w') as handle:
        for rec in SeqIO.parse(fasta, 'fasta'):
            if rec.id in keeplist:
                SeqIO.write(rec, handle, 'fasta')


def writetrf(alignDF: pd.DataFrame = None, outtab: str = None) -> str:
    """
    Write alignment dataframe to LASTZ tab format.

    Converts a dataframe of alignments (as processed by trfFilter) back to the
    LASTZ tabular format and writes it to a file with '.trf' suffix.

    Parameters
    ----------
    alignDF : pd.DataFrame
        DataFrame containing alignment data.
    outtab : str
        Base path for output file (will append '.trf').

    Returns
    -------
    str
        Path to the output file.
    """
    # Define output file path with .trf extension
    outfile = outtab + '.trf'

    # Write alignments in LASTZ tabular format
    with open(outfile, 'w') as handle:
        # Write header
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

        # Write each alignment as a tab-delimited row
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


def writeGFFlines(
    alnDF: pd.DataFrame = None,
    chrlens: List[Tuple[str, str]] = None,
    ftype: str = 'BHit',
) -> str:
    """
    Generate GFF3 format lines from alignment data.

    Converts alignment data to GFF3 format, including optional sequence-region
    directives based on chromosome lengths.

    Parameters
    ----------
    alnDF : pd.DataFrame
        DataFrame containing alignment data.
    chrlens : List[Tuple[str, str]], optional
        List of (chromosome_name, length) tuples for sequence-region directives.
    ftype : str, optional
        Feature type to use in GFF3, default is 'BHit'.

    Yields
    ------
    str
        Lines of GFF3 format text, one at a time.
    """
    # Generate GFF3 version line
    yield '##gff-version 3\n'

    # Add sequence-region directives if chromosome lengths are provided
    if chrlens:
        for name, maxlen in chrlens:
            yield ' '.join(['##sequence-region', str(name), '1', str(maxlen) + '\n'])

    # Add column headers as a comment
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

    # Write each alignment as a GFF feature
    for _index, row in alnDF.iterrows():
        # Construct attributes field
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

        # Yield tab-delimited GFF3 feature line
        yield '\t'.join(
            [
                row['tName'],  # seqid
                'mimeo-map',  # source
                ftype,  # type
                str(row['tStart']),  # start
                str(row['tEnd']),  # end
                str(row['score']),  # score
                row['tStrand'],  # strand
                '.',  # phase
                attributes + '\n',  # attributes
            ]
        )


def map_LZ_cmds(
    lzpath: str = 'lastz',
    pairs: List[Tuple[str, str]] = None,
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
) -> List[str]:
    """
    Generate commands to map high identity segments between genome files using LASTZ.

    This function creates a list of shell commands that identify regions with high sequence
    identity between pairs of genome files. These regions can represent candidate horizontal
    gene transfer (HGT) regions. The commands run LASTZ alignments and filter/format the results.

    Parameters
    ----------
    lzpath : str, optional
        Path to LASTZ executable, default is 'lastz'.
    pairs : list of (str, str) tuples, optional
        List of file path pairs to align. Each tuple contains (target_path, query_path).
    minIdt : float, optional
        Minimum percent identity threshold for reported alignments, default is 95.
    minLen : int, optional
        Minimum alignment length to report (bp), default is 100.
    hspthresh : int, optional
        High-scoring segment pair threshold for LASTZ sensitivity, default is 3000.
    outfile : str, optional
        Path where filtered alignments will be written.
    verbose : bool, optional
        Enable verbose output during alignment, default is False.
    lastz_format : str, optional
        Format string for LASTZ output fields.
    step_size : int, optional
        Step size parameter for LASTZ seed generation, default is 1.
    strand_mode : str, optional
        Strand to align: 'plus', 'minus', or 'both', default is 'both'.
    chain : bool, optional
        Enable chaining of alignments, default is True.
    gapped : bool, optional
        Enable gapped extensions, default is True.

    Returns
    -------
    list of str
        Shell commands to run LASTZ alignments and process the results.

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

    logging.info(f'Using LASTZ at {lzpath}')
    logging.info(f'LASTZ format: {lastz_format}')

    # Determine platform for sed command compatibility
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
    lzpath: str = 'lastz',
    bdtlsPath: str = 'bedtools',
    Adir: str = None,
    Bdir: str = None,
    pairs: List[Tuple[str, str]] = None,
    outtab: str = None,
    outgff: str = None,
    minIdt: float = 60,
    minLen: int = 100,
    hspthresh: int = 3000,
    minCov: int = 5,
    AchrmLens: str = None,
    reuseTab: bool = False,
    label: str = 'B_repeats',
    prefix: str = None,
    verbose: bool = False,
) -> List[str]:
    """
    Generate commands to identify B-genome regions with high coverage in A-genome.

    Screen genome A (target) for features which are high copy in genome B (query).
    This function creates shell commands that align B genome segments to A genome,
    then identify regions in A that are covered by multiple B segments, which may
    represent repetitive elements.

    Parameters
    ----------
    lzpath : str, optional
        Path to LASTZ executable, default is 'lastz'.
    bdtlsPath : str, optional
        Path to bedtools executable, default is 'bedtools'.
    Adir : str
        Directory containing A genome sequence files.
    Bdir : str
        Directory containing B genome sequence files.
    pairs : List[Tuple[str, str]]
        List of (A path, B path) tuples for alignment pairs.
    outtab : str
        Path for output alignment table.
    outgff : str
        Path for output GFF annotation file.
    minIdt : float, optional
        Minimum percent identity threshold, default is 60.
    minLen : int, optional
        Minimum alignment length to report (bp), default is 100.
    hspthresh : int, optional
        High-scoring segment pair threshold for LASTZ, default is 3000.
    minCov : int, optional
        Minimum depth of B-genome segments covering an A-genome region, default is 5.
    AchrmLens : str
        Path to file containing A genome chromosome lengths.
    reuseTab : bool, optional
        If True, use existing alignment file if found, default is False.
    label : str, optional
        Feature type label for GFF output, default is 'B_repeats'.
    prefix : str, optional
        Prefix for feature IDs, default is None.
    verbose : bool, optional
        Enable verbose output during alignment, default is False.

    Returns
    -------
    List[str]
        Shell commands to run the cross-species repeat identification workflow.
    """
    # Set verbosity level based on parameter
    if verbose:
        verb = 1
    else:
        verb = 0

    # Initialize command list
    cmds = []

    # Determine platform for sed command compatibility
    system = platform.system().lower()
    if system == 'darwin':  # macOS
        sed_command = "sed -i '' -e 's/%//g'"
    else:  # Linux and others
        sed_command = "sed -i 's/%//g'"

    # Generate alignments if needed (not reusing existing file)
    if not reuseTab or not os.path.isfile(outtab):
        # Write header to output file
        cmds.append(
            ' '.join(
                [
                    "echo $'#name1\\tstrand1\\tstart1\\tend1\\tname2\\tstrand2\\tstart2+\\tend2+\\tscore\\tidentity' >",
                    quote(outtab),
                ]
            )
        )

        # Process each genome pair
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

            # Remove percentage symbols
            cmds.append(' '.join([sed_command, temp_outfile]))

            # Filter alignments by min length and min identity
            # Format output and sort
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

            # Remove temporary file
            cmds.append(' '.join(['rm', temp_outfile]))

    # Set temp output files
    temp_bed = 'temp.bed'
    temp_bed_sorted = 'temp_sorted.bed'

    # Convert alignments to BED format (target coordinates only)
    cmds.append(
        ' '.join(["awk -v OFS='\\t' '!/^#/ {print $1,$3,$4;}'", outtab, '>', temp_bed])
    )

    # Remove spaces from BED file
    cmds.append(' '.join([sed_command, temp_bed]))

    # Sort BED file by chromosome, start, stop
    cmds.append(' '.join(['sort -k 1,1 -k 2n,3n', temp_bed, '>', temp_bed_sorted]))

    # Calculate coverage and filter by minimum coverage
    cmds.append(
        ' '.join(
            [
                "echo 'Generate non-zero coverage scores for target genome regions, filter for min coverage of x' >&2"
            ]
        )
    )

    # Use bedtools genomecov to calculate coverage
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

    # Re-sort filtered coverage results
    cmds.append(' '.join(['sort -k 1,1 -k 2n,3n', temp_bed, '>', temp_bed_sorted]))

    # Merge adjacent regions
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

    # Filter by minimum length and convert to GFF format
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
    lzpath: str = 'lastz',
    bdtlsPath: str = 'bedtools',
    splitSelf: bool = False,
    Adir: str = None,
    Bdir: str = None,
    pairs: List[Tuple[str, str]] = None,
    outtab: str = None,
    outgff: str = None,
    minIdt: float = 60,
    minLen: int = 100,
    hspthresh: int = 3000,
    minCov: int = 3,
    intraCov: int = 5,
    AchrmLens: str = None,
    reuseTab: bool = False,
    label: str = 'Self_repeats',
    prefix: str = None,
    verbose: bool = False,
) -> List[str]:
    """
    Generate commands for identifying repeats through genome self-alignment.

    Creates shell commands to align a genome to itself and identify regions that are
    present in multiple copies based on coverage. Optionally processes intra-chromosome
    alignments separately to handle local tandem duplications differently.

    Parameters
    ----------
    lzpath : str, optional
        Path to LASTZ executable, default is 'lastz'.
    bdtlsPath : str, optional
        Path to bedtools executable, default is 'bedtools'.
    splitSelf : bool, optional
        If True, process intra-chromosome alignments separately, default is False.
    Adir : str
        Directory containing genome sequence files.
    Bdir : str
        Directory containing genome sequence files (same as Adir for self-alignment).
    pairs : List[Tuple[str, str]]
        List of (file1, file2) tuples for alignment pairs.
    outtab : str
        Path for output alignment table.
    outgff : str
        Path for output GFF annotation file.
    minIdt : float, optional
        Minimum percent identity threshold, default is 60.
    minLen : int, optional
        Minimum alignment length to report (bp), default is 100.
    hspthresh : int, optional
        High-scoring segment pair threshold for LASTZ, default is 3000.
    minCov : int, optional
        Minimum depth of aligned segments to report repeat feature, default is 3.
    intraCov : int, optional
        Minimum depth for intra-chromosome alignments (when splitSelf=True), default is 5.
    AchrmLens : str
        Path to file containing chromosome lengths.
    reuseTab : bool, optional
        If True, use existing alignment file if found, default is False.
    label : str, optional
        Feature type label for GFF output, default is 'Self_repeats'.
    prefix : str, optional
        Prefix for feature IDs, default is None.
    verbose : bool, optional
        Enable verbose output during alignment, default is False.

    Returns
    -------
    List[str]
        Shell commands to run the self-alignment repeat identification workflow.
    """
    # Set verbosity level based on parameter
    if verbose:
        verb = 1
    else:
        verb = 0

    # Initialize command list
    cmds = []

    # Determine platform for sed command compatibility
    system = platform.system().lower()
    if system == 'darwin':  # macOS
        sed_command = "sed -i '' -e 's/%//g'"
    else:  # Linux and others
        sed_command = "sed -i 's/%//g'"

    # Set up output file for intra-chromosomal alignments if needed
    if splitSelf:
        outtab_intra = outtab + '_intra.tab'

    # Generate alignments if needed (not reusing existing file)
    if not reuseTab or not os.path.isfile(outtab):
        # Write alignment out file headers
        cmds.append(
            ' '.join(
                [
                    "echo $'#name1\\tstrand1\\tstart1\\tend1\\tname2\\tstrand2\\tstart2+\\tend2+\\tscore\\tidentity' >",
                    quote(outtab),
                ]
            )
        )

        # Create separate file for intra-chromosomal alignments if using split mode
        if splitSelf:
            outtab_intra = outtab + '_intra.tab'
            cmds.append(
                ' '.join(
                    [
                        "echo $'#name1\\tstrand1\\tstart1\\tend1\\tname2\\tstrand2\\tstart2+\\tend2+\\tscore\\tidentity' >",
                        quote(outtab_intra),
                    ]
                )
            )

        # Process each sequence pair
        for A, B in pairs:
            if A != B or not splitSelf:
                # Regular alignment (inter-chromosomal or non-split mode)
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

                # Remove percentage symbols
                cmds.append(' '.join([sed_command, temp_outfile]))

                # Filter alignments by min length and min identity
                # Format output and sort
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

                # Remove temporary file
                cmds.append(' '.join(['rm', temp_outfile]))

            else:
                # Special handling for intra-chromosomal alignments in split mode
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

                # Remove percentage symbols
                cmds.append(' '.join([sed_command, temp_outfile]))

                # Filter alignments by min length and min identity
                # Format output and sort
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

                # Remove temporary file
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

    # Convert alignments to BED format (target coordinates only)
    cmds.append(
        ' '.join(["awk -v OFS='\\t' '!/^#/ {print $1,$3,$4;}'", outtab, '>', temp_bed])
    )

    # Remove spaces from BED file
    cmds.append(' '.join([sed_command, temp_bed]))

    # Sort BED file by chromosome, start, stop
    cmds.append(' '.join(['sort -k 1,1 -k 2n,3n', temp_bed, '>', temp_bed_sorted]))

    # Calculate coverage and filter by minimum coverage
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

    # Re-sort filtered coverage results
    cmds.append(' '.join(['sort -k 1,1 -k 2n,3n', temp_bed, '>', temp_bed_sorted]))

    # Merge adjacent regions
    cmds.append(' '.join([bdtlsPath, 'merge -i', temp_bed_sorted, '>', temp_bed]))

    # Initialize GFF output file
    cmds.append(
        ' '.join(
            [
                "echo $'##gff-version 3\\n#seqid\\tsource\\ttype\\tstart\\tend\\tscore\\tstrand\\tphase\\tattributes' >",
                quote(outgff),
            ]
        )
    )

    # Filter by minimum length and convert to GFF format
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

    # Process intra-chromosomal alignments separately if in split mode
    if splitSelf:
        if reuseTab and not os.path.isfile(outtab_intra) and os.path.isfile(outtab):
            logging.warning(
                "Warning: Could not find intra-chrom results file: %s \nRe-run in '--strictSelf' mode if required."
                % outtab_intra
            )
        else:
            # Apply separate coverage filtering for WITHIN chromosome hits
            cmds.append(
                ' '.join(
                    [
                        "echo 'Applying separate coverage filtering for WITHIN chromosome hits' >&2"
                    ]
                )
            )

            # Set up temporary files for intra-chromosomal processing
            temp_bed_intra = 'temp.bed'
            temp_bed_sorted_intra = 'temp_sorted.bed'

            # Convert intra-chromosomal alignments to BED format
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

            # Remove spaces from BED file
            cmds.append(' '.join([sed_command, temp_bed_intra]))

            # Sort BED file
            cmds.append(
                ' '.join(
                    ['sort -k 1,1 -k 2n,3n', temp_bed_intra, '>', temp_bed_sorted_intra]
                )
            )

            # Calculate coverage using higher threshold for intra-chromosomal alignments
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

            # Re-sort filtered coverage results
            cmds.append(
                ' '.join(
                    ['sort -k 1,1 -k 2n,3n', temp_bed_intra, '>', temp_bed_sorted_intra]
                )
            )

            # Merge adjacent regions
            cmds.append(
                ' '.join(
                    [bdtlsPath, 'merge -i', temp_bed_sorted_intra, '>', temp_bed_intra]
                )
            )

            # Filter by minimum length and convert to GFF format with intra-specific label
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

    # Return the complete list of commands
    return cmds
