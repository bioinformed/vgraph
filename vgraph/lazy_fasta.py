"""Read data from a pysam.Fastafile using lazy queries and a block-level LRU cache."""


from vgraph.lru_cache import LRUCache


class LazyFastaContig:
    """Read data from a pysam.Fastafile using lazy queries and a block-level LRU cache.

    Args:
        fasta (pysam.Fastafile): input FASTA
        contig (str): contig name
        block_size (int): cache block size in bytes.  Default is 10240 (10 kb).
        cache_size (int): cache size in blocks.  Default is 64.

    """
    def __init__(self, fasta, contig, block_size=10240, cache_size=64):
        """Build a new LazyFastaContig."""
        self.fasta = fasta
        self.contig = contig
        self.reference_length = fasta.get_reference_length(contig)
        self.block_size = block_size
        self.cache = LRUCache(cache_size)

    def __len__(self):
        """Return the contig length."""
        return self.reference_length

    def __getitem__(self, index):
        """Retrieve a slice of the FASTA sequence.

        Args:
            index (int or slice): position or slice to retrieve

        Returns:
            str: FASTA sequence for contig and index

        """
        if isinstance(index, int):
            index = slice(index, index + 1)
        elif not isinstance(index, slice):
            raise TypeError('int or slice required')

        if index.start < 0 or index.stop > self.reference_length:
            raise ValueError('indices out of range')

        result = ''.join(self._get_blocks(index))

        # Use as sanity check when debugging -- will be removed after verification
        if 0:
            expected = self.fasta.fetch(self.contig, index.start, index.stop)
            assert result == expected

        if index.step != 1:
            result = result[::index.step]

        return result

    def _get_blocks(self, index):
        """Yield blocks for specified index slice.

        Uses the LRU cache to avoid unnecessary lookups.

        Args:
            index (slice): slice to retrieve

        Yields:
            str: blocks of FASTA results, trimmed to bounds specifed in index

        """
        block_size  = self.block_size
        start_block = (index.start // block_size)
        stop_block  = ((index.stop - 1)  // block_size)

        for block_index in range(start_block, stop_block + 1):
            block = self._get_block(block_index)
            block_start = block_index * block_size
            block_stop = block_start + len(block)

            # Trim right if needed
            if block_stop > index.stop:
                # Uses negative indexing to trim end
                block = block[:index.stop - block_stop]

            # Trim left if needed
            if block_start < index.start:
                block = block[index.start - block_start:]

            yield block

    def _get_block(self, block_index):
        """Get a single block.

        Block is preferably retrieved from the cache, if it is present.  Otherwise, it is queried
        from the underlying FASTA file.

        Args:
            block_index (int): block index correspding to FASTA index from
                               `block_index*block_size` to `(block_index+1)*block_size.`

        Returns:
            str: FASTA sequence for block of length block_size.

        """
        block = self.cache.get(block_index)

        if block is None:
            block = self.fasta.fetch(self.contig, block_index * self.block_size, (block_index + 1) * self.block_size)
            self.cache[block_index] = block

        return block
