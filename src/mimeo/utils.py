"""
Utility functions for genome and sequence analysis in the mimeo package.

This module provides helper functions for common operations used throughout the mimeo
package, including:
- Handling sequence files and directories
- Processing and validating FASTA formats
- Managing sequence alignments and pairs
- Executing system commands with error handling
- Calculating sequence properties
- Handling temporary files and directories

These functions support the core operations of the mimeo package for horizontal
gene transfer analysis and genomic sequence comparison.
"""

from collections import Counter
from datetime import datetime, timezone
import glob
import os
import shutil
import subprocess
import sys
import tempfile
from typing import List, Optional, Tuple, Union
import logging

from Bio import SeqIO


def import_pairs(
    file: str = None, Adir: str = None, Bdir: str = None
) -> List[Tuple[str, str]]:
    """
    Import sequence pairs from a file for alignment.

    Reads a tab-delimited file containing pairs of sequence files to align.
    The file should have two columns representing files from genome A and B.
    Lines starting with # are treated as comments and ignored.

    Parameters
    ----------
    file : str
        Path to tab-delimited file listing sequence pairs.
    Adir : str
        Directory containing genome A sequences.
    Bdir : str
        Directory containing genome B sequences.

    Returns
    -------
    List[Tuple[str, str]]
        List of (A path, B path) tuples for each sequence pair.
    """
    pairs = []
    with open(file) as f:
        for line in f:
            li = line.strip()
            if not li.startswith('#'):
                A, B = li.split()[:2]
                pairs.append((os.path.join(Adir, A), os.path.join(Bdir, B)))
    return pairs


def get_all_pairs(
    Adir: Optional[str] = None, Bdir: Optional[str] = None
) -> List[Tuple[str, str]]:
    """
    Generate all possible sequence file pairs between two directories.

    Creates a list of all possible combinations of files from directory A and B.
    If only A is provided, creates self-alignment pairs from all files in A.

    Parameters
    ----------
    Adir : str, optional
        Directory containing genome A sequences.
    Bdir : str, optional
        Directory containing genome B sequences.

    Returns
    -------
    List[Tuple[str, str]]
        List of (A path, B path) tuples for each sequence pair.

    Raises
    ------
    SystemExit
        If no directories are provided.
    """
    pairs = []
    if Adir and Bdir:
        # Cross-species alignment: pair each file from A with each file from B
        for A in glob.glob(os.path.join(Adir, '*')):
            for B in glob.glob(os.path.join(Bdir, '*')):
                pairs.append((A, B))
    elif Adir:
        # Self-alignment: pair each file from A with each file from A
        logging.info('Compose self-genome alignment pairs.')
        for A in glob.glob(os.path.join(Adir, '*')):
            for B in glob.glob(os.path.join(Adir, '*')):
                pairs.append((A, B))
    else:
        logging.error('Need at least one seq directory to compose alignment pairs.')
        sys.exit(1)
    return pairs


def _write_script(cmds: List[str], script: str) -> None:
    """
    Write commands into a bash script file.

    Parameters
    ----------
    cmds : List[str]
        List of shell commands to write.
    script : str
        Filename to write the shell script to.

    Returns
    -------
    None
        Creates a shell script with the provided commands.
    """
    f = open(script, 'w+')
    for cmd in cmds:
        print(cmd, file=f)
    f.close()


def syscall(
    cmd: Union[str, list],
    verbose: bool = False,
    timeout: Optional[float] = None,
    encoding: str = 'utf-8',
    shell: bool = True,
) -> str:
    """
    Execute a system command with comprehensive error handling.

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
    logging.debug(f'Running command: {cmd}')

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
        logging.debug(f'{decoded_output}')

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


def run_cmd(cmds: List[str], verbose: bool = False, keeptemp: bool = False) -> None:
    """
    Write commands to a script and execute it in a temporary directory.

    Creates a temporary directory, writes commands to a shell script,
    executes the script, and optionally cleans up the temporary directory.

    Parameters
    ----------
    cmds : List[str]
        List of shell commands to execute.
    verbose : bool, optional
        If True, output from command execution will be displayed.
    keeptemp : bool, optional
        If True, temporary directory will not be removed after execution.

    Returns
    -------
    None
        Executes the commands in the temporary directory and handles cleanup.
        If keeptemp is True, the temporary directory will be retained.
    """
    # Create temporary directory for execution
    tmpdir = tempfile.mkdtemp(prefix='tmp.', dir=os.getcwd())
    original_dir = os.getcwd()

    # Change to temporary directory
    os.chdir(tmpdir)

    # Write commands to script file
    script = 'run_jobs.sh'
    _write_script(cmds, script)

    # Execute the script
    syscall('bash ' + script, verbose=verbose)

    # Return to original directory
    os.chdir(original_dir)

    # Clean up temporary directory unless keeptemp flag is set
    if not keeptemp:
        shutil.rmtree(tmpdir)


def getTimestring() -> str:
    """
    Return integer-only string of current datetime with milliseconds.

    Generates a timestamp string in format YYYYMMDDHHMMSSMMM where MMM represents
    milliseconds. Useful for creating unique file/directory names.

    Returns
    -------
    str
        Current time formatted as YYYYMMDDHHMMSSMMM.
    """
    (dt, micro) = datetime.now(timezone.utc).strftime('%Y%m%d%H%M%S.%f').split('.')
    dt = '%s%03d' % (dt, int(micro) // 1000)
    return dt


def splitFasta(infile: str, outdir: str, unique: bool = True) -> None:
    """
    Split a multiFASTA file into individual FASTA files.

    Creates one FASTA file per sequence in the input file, with each output file
    named after the sequence ID and placed in the specified directory.

    Parameters
    ----------
    infile : str
        Path to input multiFASTA file.
    outdir : str
        Directory to write individual FASTA files to.
    unique : bool, optional
        If True, verifies that all sequence IDs are unique and exits if duplicates found.

    Returns
    -------
    None
        Creates individual FASTA files in the specified output directory.

    Raises
    ------
    SystemExit
        If unique=True and non-unique sequence IDs are found.
    """
    seen = []
    for rec in SeqIO.parse(infile, 'fasta'):
        if str(rec.id) in seen and unique:
            logging.error('Non-unique name in genome: %s. Quitting.' % str(rec.id))
            sys.exit(1)
        else:
            seen.append(str(rec.id))
        outfile = os.path.join(outdir, rec.id + '.fa')
        with open(outfile, 'w') as handle:
            SeqIO.write(rec, handle, 'fasta')


def isfile(path: str) -> str:
    """
    Check if a file exists and return its absolute path.

    Parameters
    ----------
    path : str
        Path to file to check.

    Returns
    -------
    str
        Absolute path to the file.

    Raises
    ------
    SystemExit
        If file does not exist.
    """
    path = os.path.abspath(path)
    if not os.path.isfile(path):
        logging.error('Input file not found: %s' % path)
        sys.exit(1)
    else:
        return path


def set_paths(
    adir: Optional[str] = None,
    bdir: Optional[str] = None,
    afasta: Optional[str] = None,
    bfasta: Optional[str] = None,
    outdir: Optional[str] = None,
    outtab: Optional[str] = None,
    gffout: Optional[str] = None,
    suppresBdir: bool = False,
    runtrf: Optional[float] = None,
) -> Tuple[str, Optional[str], str, Optional[str], Optional[str], Optional[str]]:
    """
    Set up directory and file paths for mimeo analyses.

    Handles the creation of necessary directories, splitting of multiFASTA files,
    and setting output file paths. Creates temporary directories as needed.

    Parameters
    ----------
    adir : str, optional
        Directory containing genome A sequence files.
    bdir : str, optional
        Directory containing genome B sequence files.
    afasta : str, optional
        Path to genome A multiFASTA file.
    bfasta : str, optional
        Path to genome B multiFASTA file.
    outdir : str, optional
        Directory for output files.
    outtab : str, optional
        Name of alignment output table file.
    gffout : str, optional
        Name of GFF output file.
    suppresBdir : bool, optional
        If True, genome B directory is not required.
    runtrf : float, optional
        If specified, creates a temporary directory for TRF analysis.

    Returns
    -------
    Tuple[str, Optional[str], str, Optional[str], Optional[str], Optional[str]]
        Tuple containing:
        - adir: Path to genome A directory
        - bdir: Path to genome B directory (or None if suppressed)
        - outdir: Path to output directory
        - outtab: Path to alignment output file (or None)
        - gffout: Path to GFF output file (or None)
        - tempdir: Path to temporary directory (or None if not needed).

    Raises
    ------
    SystemExit
        If required genome files are not provided.
    """
    # Determine if we need a temporary directory
    if not adir:
        # Make temp directory if A directory not specified
        tempdir = os.path.join(os.getcwd(), 'temp_' + getTimestring())
        os.makedirs(tempdir)
    elif not bdir and not suppresBdir:
        # Make temp directory if B directory not specified (and not suppressed)
        tempdir = os.path.join(os.getcwd(), 'temp_' + getTimestring())
        os.makedirs(tempdir)
    elif runtrf:
        # Make temp directory for TRF analysis
        tempdir = os.path.join(os.getcwd(), 'temp_' + getTimestring())
        os.makedirs(tempdir)
    else:
        tempdir = None

    # Setup and validate A genome directory
    if adir:
        adir = os.path.abspath(adir)
        if not os.path.isdir(adir):
            logging.info('Creating Adir: %s' % adir)
            os.makedirs(adir)
            if not afasta:
                logging.error('No A-genome fasta file provided. Quitting.')
                sys.exit(1)
    else:
        # Create A genome directory in temp location
        adir = os.path.join(tempdir, 'A_genome_split')
        os.makedirs(adir)

    # Setup and validate B genome directory (if needed)
    if bdir:
        bdir = os.path.abspath(bdir)
        if not os.path.isdir(bdir):
            logging.info('Creating Bdir: %s' % bdir)
            os.makedirs(bdir)
            if not bfasta:
                logging.error('No B-genome fasta file provided. Quitting.')
                sys.exit(1)
    elif not suppresBdir:
        # Create B genome directory in temp location
        bdir = os.path.join(tempdir, 'B_genome_split')
        os.makedirs(bdir)

    # Split genome multiFASTA files into individual files if provided
    if afasta:
        if os.path.isfile(afasta):
            splitFasta(afasta, adir)
        else:
            logging.error('A-genome fasta not found at path: %s' % afasta)

    if bfasta:
        if os.path.isfile(bfasta):
            splitFasta(bfasta, bdir)
        elif not suppresBdir:
            logging.error('B-genome fasta not found at path: %s' % bfasta)

    # Set up output directory
    if outdir:
        outdir = os.path.abspath(outdir)
        if not os.path.isdir(outdir):
            logging.info('Create output directory: %s' % outdir)
            os.makedirs(outdir)
    else:
        outdir = os.getcwd()

    # Set up output file paths
    if outtab:
        outtab = os.path.join(outdir, outtab)
        if os.path.isfile(outtab):
            logging.info('Previous alignment found: %s' % outtab)

    if gffout:
        gffout = os.path.join(outdir, gffout)

    # Return all paths
    return adir, bdir, outdir, outtab, gffout, tempdir


def checkUniqueID(records: List) -> None:
    """
    Check that IDs for sequence records are unique.

    Parameters
    ----------
    records : List
        List of BioPython SeqRecord objects.

    Returns
    -------
    None
        Prints a message if duplicate IDs are found and exits the program.
        If all IDs are unique, the function does nothing.

    Raises
    ------
    SystemExit
        If duplicate sequence IDs are found.
    """
    seqIDs = [records[x].id for x in range(len(records))]
    IDcounts = Counter(seqIDs)
    duplicates = [k for k, v in IDcounts.items() if v > 1]
    if duplicates:
        logging.error(f'Input sequence IDs not unique:\n{duplicates}\n\nQuitting.')
        sys.exit(1)
    else:
        pass


def chromlens(
    seqDir: str = None, outfile: Optional[str] = None
) -> List[Tuple[str, str]]:
    """
    Get chromosome lengths from a directory of FASTA files.

    Reads all FASTA files in the specified directory and calculates the length
    of each sequence. Optionally writes the results to a tab-delimited file.

    Parameters
    ----------
    seqDir : str
        Directory containing FASTA files.
    outfile : str, optional
        If provided, writes lengths to this file.

    Returns
    -------
    List[Tuple[str, str]]
        List of (chromosome_id, length) tuples.

    Raises
    ------
    SystemExit
        If no sequences found or sequence IDs are not unique.
    """
    # Read all sequence records from files in directory
    records = []
    for f in glob.glob(os.path.join(seqDir, '*')):
        records += list(SeqIO.parse(f, 'fasta'))

    # Check that records were found
    if not records:
        logging.error(
            'No sequences found in %s \n Cannot calculate seq lengths.' % seqDir
        )
        sys.exit(1)

    # Check that sequence IDs are unique
    checkUniqueID(records)

    # Calculate lengths for each sequence
    chrlens = []
    for rec in records:
        chrlens.append((str(rec.id), str(len(rec.seq))))

    # Sort by chromosome ID
    chrlens.sort(key=lambda x: x[0])

    # Write to output file if specified
    if outfile:
        with open(outfile, 'w') as handle:
            for x, y in chrlens:
                handle.write('\t'.join([x, y]) + '\n')

    return chrlens


def missing_tool(tool_name: str) -> List[str]:
    """
    Check if a tool is available in the system PATH.

    Parameters
    ----------
    tool_name : str
        Name of the tool executable to check.

    Returns
    -------
    List[str]
        Empty list if tool is found, list containing tool_name if not found.
    """
    path = shutil.which(tool_name)
    if path is None:
        return [tool_name]
    else:
        return []
