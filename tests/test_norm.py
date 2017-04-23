import doctest

from vgraph.norm import trim_common_prefixes, trim_common_suffixes

def test_trim_common_prefixes():
    '''
    >>> trim_common_prefixes([])
    (0, [])
    >>> trim_common_prefixes(['AA'])
    (0, ['AA'])
    >>> trim_common_prefixes(['A','AB'])
    (1, ['', 'B'])
    >>> trim_common_prefixes(['A','BA'])
    (0, ['A', 'BA'])
    >>> trim_common_prefixes(['A','AB','AAAB'])
    (1, ['', 'B', 'AAB'])
    >>> trim_common_prefixes(['AB','ABA','ABAA'])
    (2, ['', 'A', 'AA'])
    >>> trim_common_prefixes(['AB','ABA','ABAA','C'])
    (0, ['AB', 'ABA', 'ABAA', 'C'])
    >>> trim_common_prefixes(['AAAAAAAAA','AAAAAAAAA','AAAAAAAAA','AAAAAAAAA'])
    (9, ['', '', '', ''])
    '''


def test_trim_common_suffixes():
    '''
    >>> trim_common_suffixes([])
    (0, [])
    >>> trim_common_suffixes(['AA'])
    (0, ['AA'])
    >>> trim_common_suffixes(['A','BA'])
    (1, ['', 'B'])
    >>> trim_common_suffixes(['A','AB'])
    (0, ['A', 'AB'])
    >>> trim_common_suffixes(['A','BA','BAAA'])
    (1, ['', 'B', 'BAA'])
    >>> trim_common_suffixes(['BA','ABA','AABA'])
    (2, ['', 'A', 'AA'])
    >>> trim_common_suffixes(['AB','ABA','ABAA','C'])
    (0, ['AB', 'ABA', 'ABAA', 'C'])
    >>> trim_common_suffixes(['AAAAAAAAA','AAAAAAAAA','AAAAAAAAA','AAAAAAAAA'])
    (9, ['', '', '', ''])
    '''


if __name__ == '__main__':
    doctest.testmod()
