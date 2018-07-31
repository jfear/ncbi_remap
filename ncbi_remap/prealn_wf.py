import os
from pathlib import Path

from .io import remove_id


def check_indicator_file(store, key, pattern, srx, srr, **kwargs):
    """Check indicator file for presense.

    Parameters
    ----------
    store : pd.io.pytables.HDFStore
        The data store to save to.
    key : str
        Where in the store to save flag.
    pattern : str
        File naming pattern for indicator file.
    srx : str
        SRX accession
    srr : str
        SRR accession
    **kwargs
        Keywords needed to fill the pattern.

    """
    fill = {
        'srx': srx,
        'srr': srr
    }
    fill.update(kwargs)

    ab = pattern.format(**fill)
    if os.path.exists(ab):
        store.pop('prealn/queue', (srx, srr))
        flags = store[key].copy()
        flags[(srx, srr)] = True
        store[key] = flags
        return True


def check_flag_file(store, key, pattern, srx, srr, **kwargs):
    """Reads a file for its stored flag.

    This is currently sued for the layout and strand files.

    Parameters
    ----------
    store : pd.io.pytables.HDFStore
        The data store to save to.
    key : str
        Where in the store to save flag.
    pattern : str
        File naming pattern for the ALIGNEMNT_BAD file.
    srx : str
        SRX accession
    srr : str
        SRR accession
    **kwargs
        Keywords needed to fill the pattern.

    """
    fill = {
        'srx': srx,
        'srr': srr
    }
    fill.update(kwargs)

    with open(pattern.format(**fill)) as fh:
        flag = fh.read().strip()
        dat = store[key]
        dat[(srx, srr)] = flag
        store[key] = dat


def check_download(store, pattern, srx, srr, **kwargs):
    """Checks FASTQ logs determine if bad DB entry.

    There are some SRRs that are no longer accessible in the SRA database. I
    get an error in the fastq-dump log that I want to translate into the
    DOWNLOAD_BAD flag.

    Parameters
    ----------
    store : .io.HDFStore
        The data store to save to.
    pattern : str
        File naming pattern for the DOWNLOAD_BAD file.
    srx : str
        SRX accession
    srr : str
        SRR accession
    **kwargs
        Keywords needed to fill the pattern.

    """
    fill = {
        'srx': srx,
        'srr': srr
    }
    fill.update(kwargs)

    logName = Path(pattern.format(**fill))
    flagBad = Path(logName.parent, 'DOWNLOAD_BAD')
    if logName.exists():
        dat = logName.read_text()
        if 'failed to resolve accession' in dat:
            flagBad.touch()

    return check_indicator_file(store, pattern, 'prealn/download_bad',
                                srx, srr, **kwargs)


