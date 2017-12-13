import os
import hashlib


def check_fastq(fn):
    """Checks if a FASTQ file actually has data in it."""
    try:
        assert os.stat(fn).st_size > 10000
        return True
    except (FileNotFoundError, AssertionError) as err:
        return False
    except:
        raise


def md5sum(fn):
    """Calculates the md5sum of the raw fastq file."""
    md5 = hashlib.md5()
    with open(fn, 'rb') as f:
        for chunk in iter(lambda: f.read(128 * md5.block_size), b''):
            md5.update(chunk)
    return md5.hexdigest()


def fastq_stats(fn):
    """Calculates the number of reads (libsize) and the average read length."""
    libsize, lens = 0, 0

    with open(fn, 'r') as fh:
        for i, row in enumerate(fh):
            if (i % 4 == 1):
                libsize += 1
                lens += len(row.strip())

    avgLen = lens / libsize
    return libsize, avgLen

