"""histogram -- tools for constructing and manipulating histograms

Language: Python 3
Patrick J. Fasano
University of Notre Dame

+ 09/26/19 (pjf): Created, adding HistogramMapping.
"""

import bisect
import collections.abc
import math

################################################################
# histogram binning
################################################################

class BinMapping(collections.abc.MutableMapping):
    """Custom mapping for where all keys within bins are mapped to the same value.

    This can be used to construct a histogram by automatically collecting all
    values with "similar" keys together:

        >>> binmap = BinMapping(keys=[0., 1., 2.])
        >>> binmap[1.2] += 0.2
        >>> binmap[0.8] += 0.25
        >>> binmap[1.75] += 0.8
        >>> print(binmap)
        {0.0: 0, 1.0: 0.45, 2.0: 0.8}
    """
    @classmethod
    def create_bisection_bins(cls, keys):
        """Create default bins for histogram (by bisection between keys).

        Arguments:
            keys (list of numeric): values for centers of bins
        Returns:
            (list of numeric): bin partitions
        """
        bins = [-math.inf]
        for i in range(len(keys)-1):
            bins.append((keys[i]+keys[i+1])/2)
        bins += [math.inf]

        return bins

    def __init__(self, keys, bins=None, values=None):
        self.__keys = list(keys)
        self.__len = len(keys)
        if values is None:
            self.clear()
        else:
            if len(values) != len(self):
                raise ValueError("values have wrong shape: {}".format(values))
            self.__data = list(values[:])
        if bins is None:
            self.__bins = BinMapping.create_bisection_bins(keys)
        else:
            if len(bins) > len(keys)+1:
                raise ValueError("too many bins for keys")
            elif len(bins) == len(keys)+1:
                self.__bins = bins[:]
            elif len(bins) == len(keys):
                if bins[0] > keys[0]:
                    # add dummy bin to beginning of list
                    self.__bins = [-math.inf] + bins
                elif bins[-1] < keys[-1]:
                    # add dummy bin to end of list
                    self.__bins = bins + [math.inf]
                else:
                    raise ValueError("invalid bins: {}".format(bins))
            elif len(bins) == len(keys)-1:
                self.__bins = [-math.inf] + bins + [math.inf]
            else:  # len(bins) < len(keys)-1
                raise ValueError("too few bins for keys")

        # validate bins
        for i in range(len(keys)):
            if (self.__keys[i] < self.__bins[i]) or (self.__keys[i] > self.__bins[i+1]):
                raise ValueError("invalid bin: ({}, {}, {})".format(
                    self.__bins[i], self.__keys[i], self.__bins[i+1]
                ))
            if sorted(self.__bins) != self.__bins:
                raise ValueError("bins not sorted: {}".format(self.__bins))

    def __repr__(self):
        return "{class_name}(keys={keys}, bins={bins}, values={values})".format(
            class_name=self.__class__.__name__,
            keys=self.__keys, bins=self.__bins, values=self.__data
        )

    def __str__(self):
        return str(dict(self))

    def clear(self):
        """Reset all bins to zero."""
        self.__data = [0] * len(self)

    def key_index(self, key):
        """Get index associated with key.

        Arguments:
            key (number): key to lookup
        Returns:
            (int): index for bin
        Raises:
            KeyError: if key falls outside bins
        """
        index = bisect.bisect(self.__bins, key) - 1
        if (index < 0) or (index >= self.__len):
            raise KeyError(key)

        return index

    def key(self, key):
        """Get canonical key associated with key.

        Arguments:
            key (number): key to canonicalize
        Returns:
            (number): canonical key
        Raises:
            KeyError: if key falls outside bins
        """
        return self.__keys[self.key_index(key)]

    def __len__(self):
        return self.__len

    def __getitem__(self, key):
        return self.__data[self.key_index(key)]

    def __setitem__(self, key, value):
        self.__data[self.key_index(key)] = value

    def __delitem__(self, key):
        """Not implemented."""
        raise NotImplementedError("removing bins not supported")

    def __iter__(self):
        return self.__keys.__iter__()

    def __contains__(self, key):
        """Check if key in bin range."""
        try:
            self.key_index(key)
        except KeyError:
            return False
        return True

if __name__ == "__main__":
    import timeit
    import numpy as np

    print(min(timeit.repeat(
        setup='binmap=BinMapping(keys=np.linspace(-5,5,101))',
        stmt='binmap[np.random.normal()] += 1',
        globals=globals(),
        number=100000,
        repeat=5
    )))