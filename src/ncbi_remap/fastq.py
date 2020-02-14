import os
import re
from typing import TextIO, Tuple
import hashlib


def fastq_is_empty(file_object: TextIO) -> bool:
    """Checks if a FASTQ file actually has data in it."""
    file_object.seek(1000, 0)
    if file_object.tell() == 0:
        return True


def fastq_read_stats(file_object: TextIO) -> Tuple[int, float]:
    """Calculates the number of reads (libsize) and the average read length."""
    libsize, lens = 0, 0
    for i, row in enumerate(file_object):
        if i % 4 == 1:
            libsize += 1
            lens += len(row.strip())

    avgLen = lens / libsize
    return libsize, avgLen


def fastq_is_abi_solid(file_object: TextIO) -> bool:
    """Look at read and determine if using abi colorspace.

    An Abi Solid read has the format:
        @SRR######.1 solid0527_####### length=50
        T120202210232000020002.00301000012.100...00........
        +SRR######.1 solid0527_####### length=50
        !/<%2/:%*)-%%0'--'')/.!%('1'%),+/%!&',!!!'+!!!!!!!!

    """
    _, seq = file_object.readline(), file_object.readline()
    if re.match(r"^T[\d\.]+$", seq):
        return True
