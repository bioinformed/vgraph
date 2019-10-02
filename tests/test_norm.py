"""Test variant normalization."""

from vgraph.norm import normalize_seq, trim_common_prefixes, trim_common_suffixes


def test_normalize_seq():
    """Test sequence normalization."""
    assert normalize_seq('') == ''
    assert normalize_seq('ACGT') == 'ACGT'
    assert normalize_seq('ACGTNacgtNRYSWKMBDHV') == 'ACGTNACGTNNNNNNNNNNN'


def test_trim_common_prefixes():
    """Test trim_common_prefixes."""
    assert trim_common_prefixes([]) == (0, [])
    assert trim_common_prefixes(['AA']) == (0, ['AA'])
    assert trim_common_prefixes(['A', 'AB']) == (1, ['', 'B'])
    assert trim_common_prefixes(['A', 'BA']) == (0, ['A', 'BA'])
    assert trim_common_prefixes(['A', 'AB', 'AAAB']) == (1, ['', 'B', 'AAB'])
    assert trim_common_prefixes(['AB', 'ABA', 'ABAA']) == (2, ['', 'A', 'AA'])
    assert trim_common_prefixes(['AB', 'ABA', 'ABAA', 'C']) == (0, ['AB', 'ABA', 'ABAA', 'C'])
    assert trim_common_prefixes(['AAAAAAAAA', 'AAAAAAAAA', 'AAAAAAAAA', 'AAAAAAAAA']) == (9, ['', '', '', ''])


def test_trim_common_suffixes():
    """Test trim_common_suffixes."""
    assert trim_common_suffixes([]) == (0, [])
    assert trim_common_suffixes(['AA']) == (0, ['AA'])
    assert trim_common_suffixes(['A', 'BA']) == (1, ['', 'B'])
    assert trim_common_suffixes(['A', 'AB']) == (0, ['A', 'AB'])
    assert trim_common_suffixes(['A', 'BA', 'BAAA']) == (1, ['', 'B', 'BAA'])
    assert trim_common_suffixes(['BA', 'ABA', 'AABA']) == (2, ['', 'A', 'AA'])
    assert trim_common_suffixes(['AB', 'ABA', 'ABAA', 'C']) == (0, ['AB', 'ABA', 'ABAA', 'C'])
    assert trim_common_suffixes(['AAAAAAAAA', 'AAAAAAAAA', 'AAAAAAAAA', 'AAAAAAAAA']) == (9, ['', '', '', ''])
