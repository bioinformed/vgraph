from collections import OrderedDict


class LRUCache:
    """Bounded cache with Least Recently Used (LRU) replacement policy

    Attrs:
        capacity (int): Cache capacity

    Examples:
        >>> cache = LRUCache(2)
        >>> 'a' not in cache
        True
        >>> cache.get('a') is None
        True
        >>> cache['a'] = 123
        >>> cache['b'] = 456
        >>> cache['a']
        123
        >>> cache['c'] = 555
        >>> cache.get('b') is None
        True

    """
    def __init__(self, capacity):
        self._capacity = capacity
        self.cache = OrderedDict()

    def __len__(self):
        """Return size of the cache.

        Returns:
            int: cache size

        """
        return len(self.cache)

    def __getitem__(self, key):
        """Get value associated with key, raising KeyError if not present."""
        if key not in self.cache:
            raise KeyError(key)

        value = self.cache[key]
        self.cache.move_to_end(key)
        return value

    def __setitem__(self, key, value):
        """Set value for key."""
        if key in self.cache:
            self.cache.move_to_end(key)
        else:
            if len(self.cache) >= self._capacity:
                 self.cache.popitem(last=False)
            self.cache[key] = value

    def __delitem__(self, key):
        """Delete key, raising KeyError if not present."""
        del self.cache[key]

    def __contains__(self, key):
        """Test if cache contains key."""
        return key in self.cache

    def get(self, key, default=None):
        """Get value associated with key, returning default if not present."""
        if key not in self.cache:
            return default

        value = self.cache[key]
        self.cache.move_to_end(key)
        return value

    def clear(self):
        """Clear all items in the cache."""
        self.cache.clear()

    @property
    def capacity(self):
        """Return the capacity of the cache."""
        return self._capacity

    @capacity.setter
    def capacity(self, value):
        self._capacity = value
        if self._capacity < 0:
            raise ValueError('capacity must be >= 0')
        while len(self.cache) > self._capacity:
            self.cache.popitem(last=False)
