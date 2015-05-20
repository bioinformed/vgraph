# -*- coding: utf-8 -*-

__abstract__ = 'smart file related utility functions'
__copyright__ = 'Copyright (c) 2015, Kevin Jacobs'

# Original copyright & software license
#
# Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health &
# Human Services. Funded by NCI under Contract N01-CO-12400.
#
# BioInformed LLC and the U.S. Department of Health & Human Services
# (“Copyright Holders”) hereby grant to the public a perpetual, irrevocable,
# royalty-free non-exclusive, worldwide license to use, copy, distribute,
# display and prepare derivatives of the software in source or binary forms,
# with or without modification, subject only to the following conditions:
#
#   o Redistributions of source code must retain the above copyright notice,
#     this list of conditions and the following disclaimer.
#
#   o Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#
#   o Neither the name of BioInformed LLC nor the names of its contributors
#     nor the U.S. Department of Health & Human Services may be used to
#     endorse or promote products derived from this software without specific
#     prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import os
import sys
import argparse


__all__ = ['smartfile', 'namefile', 'hyphen', 'SmartfileArg']


COMPRESSED_SUFFIXES = {'gz': 'gzip', 'Z': 'gzip',
                       'bz': 'bzip2', 'bz2': 'bzip2',
                       }


def is_str(s):
    '''
    is_str(s) -> bool

    Return whether s is some flavor of string

    :param s: an object
    :rtype:   bool

    >>> is_str(True)
    False
    >>> is_str(None)
    False
    >>> is_str('abc')
    True
    >>> is_str(u'abc')
    True
    '''
    return isinstance(s, basestring)


def spawn_compressor(exe, filename, mode, bufsize=-1):
    '''
    spawn_compressor(exe, filename, mode, bufsize=-1) -> input/output file object

    Spawn a subprocess to run a (de)compressor like gzip or bzip2.

    :param       exe: executable file name (must be fully qualified or found on the current PATH)
    :type        exe: str
    :param  filename: file name or file object
    :type   filename: str or file object
    :param      mode: determine whether the file objects should be opened for input or output,
                      either 'w' or 'r'
    :type       mode: str
    :param   bufsize: buffering mode and size.  0=unbuffered, 1=linebuffered, >1 buffer size,
                      -1 default buffering (default)
    :type    bufsize: int
    :return         : file object to read from or write to
    :rtype          : file object
    '''
    if os.environ.get('LOCUS_NOSPAWN'):
        raise OSError('Spawning external processes disabled by LOCUS_NOSPAWN')

    from subprocess import Popen, PIPE

    if 'w' in mode:
        out = file(filename, mode)
        cmd = [exe, '-c']
        f = Popen(cmd, stdin=PIPE, stdout=out, bufsize=bufsize).stdin
    else:
        cmd = [exe, '-d', '-c', filename]
        f = Popen(cmd, stdout=PIPE, universal_newlines='U' in mode, bufsize=bufsize).stdout

    return f


def hyphen(filename, defaultfile):
    '''
    hyphen(filename,defaultfile)

    Allow special handling when '-' is specified as a filename.  Returns a
    default file object if input filename is '-', otherwise the specified file
    name object is returned.

    :param  defaultfile: default file name or file object
    :type   defaultfile: str or file object
    :param     filename: file name or file object
    :type      filename: str or file object
    :return            : file name or file object to read from or write to
    :rtype             : str or file object
    '''
    if not is_str(filename):
        return filename

    if filename == '-':
        return defaultfile

    return filename


def namefile(filething):
    '''
    namefile(filething) -> str

    Return a human-comprehensible name for the file or file object provided.
    Recognizes file objects with the 'name' attribute, including sys.stdin,
    sys.stdout, and sys.stderr.

    :param  filething: file or filename
    :type   filething: file object or string
    :return:           human-comprehensible file name

    >>> namefile(file('/dev/null'))
    '/dev/null'
    >>> namefile(sys.stdin)
    '<stdin>'
    >>> namefile(sys.stderr)
    '<stderr>'
    >>> namefile('/dev/null')
    '/dev/null'
    '''
    if is_str(filething):
        return filething
    elif getattr(filething, 'name', None) is not None:
        return filething.name
    else:
        return repr(filething)


_smartfile_errors = (ImportError, ValueError, OSError)


def compressed_filename(filename):
    '''
    Determine if the input file is in or needs to be in compressed format

    :param  filename: file name or file object
    :type   filename: str or file object
    :return         : compression format, if compressed, otherwise an empty string
    :rtype          : str

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
    '''
    if not is_str(filename):
        return ''

    filename = os.path.expanduser(filename)
    parts = os.path.basename(filename).split('.')
    ext = parts[-1] if parts else ''
    return COMPRESSED_SUFFIXES.get(ext, '')


def smartfile(filename, mode='r', bufsize=-1):
    '''
    Return a file object in the correct compressed format as specified, which
    is ready to read from or write to

    :param  filename: file name or file object
    :type   filename: str or file object
    :param      mode: determine whether the file objects should be opened for reading or writing,
                      either 'w' or 'r'
    :type       mode: str
    :return         : file object to read from or write to
    :rtype          : file object
    '''
    # Pass non-string filename objects back as file-objects (e.g. sys.stdin or sys.stdout)
    if not is_str(filename):
        return filename

    filename = os.path.expanduser(filename)

    if not os.path.exists(filename):
        raise IOError('file does not exist')

    comp = compressed_filename(filename)

    if not comp:
        f = file(filename, mode)
    elif comp == 'gzip':
        try:
            f = spawn_compressor(os.environ.get('LOCUS_GZIP', 'gzip'), filename, mode, bufsize=bufsize)
        except _smartfile_errors:
            import gzip
            f = gzip.GzipFile(filename, mode)
    elif comp == 'bzip2':
        try:
            f = spawn_compressor(os.environ.get('LOCUS_BZIP2', 'bzip2'), filename, mode, bufsize=bufsize)
        except _smartfile_errors:
            import bz2
            f = bz2.BZ2File(filename, mode, buffering=max(0, bufsize))
    else:
        raise ValueError('Unknown compression scheme: %s' % comp)

    return f


class SmartfileArg(object):
    '''Factory for creating file object types

    Instances of SmartfileArg are typically passed as type= arguments to the
    ArgumentParser add_argument() method.

    Keyword Arguments:
        - mode -- A string indicating how the file is to be opened. Accepts the
            same values as the builtin open() function.
        - bufsize -- The file's desired buffer size. Accepts the same values as
            the builtin open() function.
    '''
    def __init__(self, mode='r', bufsize=-1):
        self._mode = mode
        self._bufsize = bufsize

    def __call__(self, value):
        # the special argument "-" means sys.std{in,out}
        if value == '-':
            if 'r' in self._mode:
                return sys.stdin
            elif 'w' in self._mode or 'a' in self._mode:
                return sys.stdout
            else:
                msg = 'argument "-" with mode %r' % self._mode
                raise ValueError(msg)

        # all other arguments are used as file names
        try:
            return smartfile(value, self._mode, self._bufsize)
        except IOError as e:
            message = "can't open '%s': %s"
            raise argparse.ArgumentTypeError(message % (value, e))

    def __repr__(self):
        args = self._mode, self._bufsize
        args_str = ', '.join(repr(arg) for arg in args if arg != -1)
        return '%s(%s)' % (type(self).__name__, args_str)


def _test():
    import doctest
    return doctest.testmod()


if __name__ == '__main__':
    _test()
