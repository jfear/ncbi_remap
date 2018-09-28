"""Some basic code that is useful all over.

Functions
---------
grouper
    Take an iterable and a number of chunks and split the iterable into groups
    of that many chunks
"""

from itertools import zip_longest


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks

    This useful function is from the itertools website.

    Parameter
    ---------
    iterable : iterable
        any kind of iterable object
    n : int
        number of items per group
    fillvalue : any
        fills the last group with this if it ran out of items (default is None)

    Example
    -------
    >>> groups = grouper('ABCDEFG', 3)
    >>> groups
    ... [('A', 'B', 'C'), ('D', 'E', 'F'), ('G', None, None)]

    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)
