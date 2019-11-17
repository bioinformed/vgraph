"""Smarter file opener that understands compression."""

import io
import os
from subprocess import Popen, PIPE
from pathlib import Path

import pysam


def spawn_compressor(exe, filename, mode, buffering=-1, encoding=None, errors=None, newline=None):
    """Spawn a subprocess to run a (de)compressor like gzip or bzip2.

    Args:
        exe (str): executable file name (must be fully qualified or found on the current PATH)
        filename (str or file object): file name or file object
        mode (str): determine whether the file objects should be opened for input or output, either 'w' or 'r'
        buffering (int): buffering mode and size.  0=unbuffered, 1=linebuffered, >1 buffer size, -1 default buffering
        encoding: encoding method.  Default='utf-8'
        errors: how encoding and decoding errors should be handled.  Default=None (See:https://docs.python.org/3/library/functions.html#open)

    Returns:
        file object to read from or write to

    Examples:
        spawn_compressor(exe, filename, mode, buffering=-1) -> input/output file object

    """
    if 'w' in mode:
        out = open(filename, mode)
        cmd = [exe, '-c']
        f = Popen(cmd, stdin=PIPE, stdout=out, bufsize=buffering).stdin
    else:
        cmd = [exe, '-d', '-c', filename]
        f = Popen(cmd, stdout=PIPE, universal_newlines='U' in mode, bufsize=buffering).stdout

    if 'b' not in mode and 'U' not in mode:
        f = io.TextIOWrapper(f, encoding, errors, newline)

    return f


COMPRESSED_SUFFIXES = {'gz': 'gzip', 'Z': 'gzip', 'bgz': 'gzip',
                       'bz': 'bzip2', 'bz2': 'bzip2'}


def compressed_filename(filename):
    """Determine if the input file is in or needs to be in compressed format.

    Args:
        filename (str or file object): file name or file object

    Returns:
        compression format, if compressed, otherwise an empty string

    Examples:
        >>> compressed_filename('subjects.sdat')
        ''
        >>> compressed_filename('subjects.sdat.gz')
        'gzip'
        >>> compressed_filename('../subjects.sdat')
        ''
        >>> compressed_filename('../subjects.sdat.gz')
        'gzip'
        >>> compressed_filename('../subjects.sdat.bz2')
        'bzip2'

    """
    if not isinstance(filename, str):
        return ''

    filename = Path(filename).name
    parts = filename.split('.')
    ext = parts[-1] if parts else ''
    return COMPRESSED_SUFFIXES.get(ext, '')


def smartfile(filename, mode='r', buffering=-1, encoding=None, errors=None, newline=None):
    """Return a file object in the correct compressed format as specified.

    Args:
        filename (str or file object): file name or file object
        mode (str): determine whether the file objects should be opened for reading or writing, either 'w' or 'r'.  Default='r'
        buffering (int): buffering mode and size.  0=unbuffered, 1=linebuffered, >1 buffer size, -1 default buffering (default)

    Returns:
        file object to read from or write to

    """
    # Pass non-string filename objects back as file-objects (e.g. sys.stdin or sys.stdout)
    if not isinstance(filename, str):
        return filename

    filename = os.path.expanduser(filename)

    if 'r' in mode and not filename.startswith('s3:') and not Path(filename).exists():
        raise OSError('file does not exist')

    comp = compressed_filename(filename)

    if filename.startswith('s3:'):
        f = pysam.HFile(filename, mode)
        if 'b' not in mode:
            f = io.TextIOWrapper(f, encoding, errors, newline)
    elif not comp:
        f = open(filename, mode, buffering, encoding, errors, newline)
    elif comp == 'gzip':
        try:
            f = spawn_compressor('gzip', filename, mode, buffering, encoding, errors, newline)
        except (ImportError, ValueError, OSError):
            import gzip
            f = gzip.GzipFile(filename, mode, encoding, errors, newline)
    elif comp == 'bzip2':
        try:
            f = spawn_compressor('bzip2', filename, mode, buffering, encoding, errors, newline)
        except (ImportError, ValueError, OSError):
            import bz2
            f = bz2.BZ2File(filename, mode, encoding, errors, newline)
    else:
        raise ValueError('Unknown compression scheme: %s' % comp)

    return f
