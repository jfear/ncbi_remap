import os

from .io import remove_id


def check_alignment2(store, pattern, **kwargs):
    """Checks for ALIGNMENT_BAD file.

    If there is an ALIGNMENT_BAD file then remove from the queue, add to
    complete, and add to 'prealn/alignment_bad'.

    Parameters
    ----------
    store : pd.io.pytables.HDFStore
        The data store to save to.
    patter : str
        File naming pattern for the ALIGNEMNT_BAD file.
    **kwargs
        Keywords needed to fill the pattern.

    """
    ab = pattern.format(**kwargs)
    if os.path.exists(ab):
        remove_id(store, 'aln/queue', **kwargs)
        flags = store['aln/alignment_bad'].copy()
        flags[(kwargs['srx'], kwargs['srr'])] = True
        store['aln/alignment_bad'] = flags
        return True


