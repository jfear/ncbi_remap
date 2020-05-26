import os
import re
from collections import namedtuple
from gzip import open as gopen
from io import BytesIO
from itertools import zip_longest
from typing import Optional, TextIO, Union

from more_itertools import grouper


class UnequalNumberReadsException(Exception):
    """Basic exception when PE data don't have the same number of reads"""

class MixedUpReadsException(Exception):
    """Basic exception when PE reads are not in the same order"""


Read = namedtuple("Read", "h1,seq,h2,qual")


class Fastq:
    """A simple FASTQ Parser"""

    def __init__(self, R1: str, R2: Optional[str] = None):
        self.R1 = R1
        self.R2 = R2
        self.libsize = None
        self.avgReadLen = None
        self.unequal_len = None
        self.bad_ecoding = None
        self.incomplete_read = None
        self.flags = set()

    def process(self):
        """Process FASTQ File

        Iterate over the FASTQ file returning each read as a string. If the
        FASTQ is from an ABI Solid sequencer it will raise an `AbiException`.
        The script also removes reads that have:
            - Are not encoded as ASCII
            - Have unequal length Sequence and Quality
            - Incomplete reads

        The object accumulates the Library Size (libsize) and the Average
        Read Length (avgReadLen). It also provides counts for the above
        problems.

        Yields
        ======
        String representation of each read.

        Raises
        ======
        NoReadsException: If FASTQ is empty.
        AbiException: If FASTQ is from ABI Solid.

        """
        if self._is_abi():
            self.flags.add("abi_solid")
            return

        if self._is_empty(self.R1) & self._is_empty(self.R2):  # No data
            self.flags.add("download_bad")
            return

        if (~self._is_empty(self.R1) & self._is_empty(self.R2)) | ("keep_R1" in self.flags):
            self.flags.add("keep_R1")
            self.R2 = None
            yield from self._process_single_end()
            return

        if (self._is_empty(self.R1) & ~self._is_empty(self.R2)) | ("keep_R2" in self.flags):
            self.flags.add("keep_R2")
            self.R1 = self.R2
            self.R2 = None
            yield from self._process_single_end()
            return

        if ~self._is_empty(self.R1) & (self.R2 is None):
            self.flags.add("SE")
            yield from self._process_single_end()
            return

        if ~self._is_empty(self.R1) & ~self._is_empty(self.R2):
            self.flags.add("PE")
            yield from self._process_pair_end()
            return

    def open_fastq(self, fastq: Union[str, bytes] = None):
        if fastq is None:
            fastq = self.R1

        if isinstance(fastq, bytes):
            return BytesIO(fastq)
        elif fastq.endswith(".gz"):
            return gopen(fastq, "rb")
        elif fastq.endswith(".fastq") or fastq.endswith(".fq"):
            return open(fastq, "rb")
        elif isinstance(fastq, str):
            return BytesIO(fastq.encode("ascii"))

    @staticmethod
    def iter_reads(file_handle: TextIO) -> Read:
        for h1, seq, h2, qual in grouper(file_handle, 4):
            if seq is None and qual is None:
                continue
            yield Read(h1, seq, h2, qual)

    def _process_single_end(self):
        self.libsize = 0
        self.unequal_len = 0
        self.bad_ecoding = 0
        self.incomplete_read = 0
        total_len = 0

        with self.open_fastq(self.R1) as fh:
            for read in self.iter_reads(fh):
                if self._is_incomplete(read):
                    # This should only happen if file was truncated.
                    self.incomplete_read += 1
                    continue

                if self._is_wrong_encoding(read):
                    self.bad_ecoding += 1
                    continue

                if self._is_unequal_seq_qual(read):
                    self.unequal_len += 1
                    continue

                self.libsize += 1
                total_len += len(read.seq.decode("ascii").strip())
                yield self._read_to_string(read)

        self.avgReadLen = total_len / self.libsize

    def _process_pair_end(self):
        self.libsize = 0
        self.unequal_len = 0
        self.bad_ecoding = 0
        self.incomplete_read = 0
        total_len = [0, 0]

        with self.open_fastq(self.R1) as fh1, self.open_fastq(self.R2) as fh2:
            for read1, read2 in zip_longest(self.iter_reads(fh1), self.iter_reads(fh2)):
                if read1 is None:
                    self.flags.remove("PE")
                    self.flags.add("keep_R2")
                    raise UnequalNumberReadsException

                if read2 is None:
                    self.flags.remove("PE")
                    self.flags.add("keep_R1")
                    raise UnequalNumberReadsException

                if self._is_incomplete(read1) | self._is_incomplete(read2):
                    # This should only happen if file was truncated.
                    self.incomplete_read += 1
                    continue

                if self._is_wrong_encoding(read1) | self._is_wrong_encoding(read2):
                    self.bad_ecoding += 1
                    continue

                if self._is_unequal_seq_qual(read1) | self._is_unequal_seq_qual(read2):
                    self.unequal_len += 1
                    continue

                if self._is_different_header(read1, read2):
                    self.flags.remove("PE")
                    if os.stat(self.R1).st_size < os.stat(self.R2).st_size:
                        self.flags.add("keep_R2")
                    else:
                        self.flags.add("keep_R1")
                    raise MixedUpReadsException

                self.libsize += 1
                total_len[0] += len(read1.seq.decode("ascii").strip())
                total_len[1] += len(read2.seq.decode("ascii").strip())
                yield self._read_to_string(read1), self._read_to_string(read2)

        self.avgReadLen = [total_len[0] / self.libsize, total_len[1] / self.libsize]

    @staticmethod
    def _is_empty(data: Union[str, bytes]) -> bool:
        if data is None:
            return True

        if isinstance(data, bytes) & len(data) > 0:
            return False

        if isinstance(data, str):
            if data.endswith(".gz") | data.endswith(".fq") | data.endswith(".fastq"):
                stat = os.stat(data)
                if stat.st_size > 1000:
                    return False
            elif len(data) > 0:
                return False

        return True

    def _is_abi(self):
        with self.open_fastq(self.R1) as fh:
            for read in self.iter_reads(fh):
                if self._is_abi_read(read):
                    return True
                break

        if self.R2:
            with self.open_fastq(self.R2) as fh:
                for read in self.iter_reads(fh):
                    if self._is_abi_read(read):
                        return True
                    break
        return False

    @staticmethod
    def _is_abi_read(read: Read) -> bool:
        """Look at read and determine if using abi colorspace.

        An Abi Solid read has the format:
            @SRR######.1 solid0527_####### length=50
            T120202210232000020002.00301000012.100...00........
            +SRR######.1 solid0527_####### length=50
            !/<%2/:%*)-%%0'--'')/.!%('1'%),+/%!&',!!!'+!!!!!!!!
        """
        if re.match(r"^T[\d\.]+$", read.seq.decode("ascii")):
            return True
        return False

    @staticmethod
    def _is_incomplete(read: Read) -> bool:
        """Determine if read has qual score"""
        if read.qual is None:
            return True
        return False

    @staticmethod
    def _is_wrong_encoding(read: Read) -> bool:
        """Determine if read is correctly encoded"""
        try:
            read.h1.decode("ascii")
            read.seq.decode("ascii")
            read.h2.decode("ascii")
            read.qual.decode("ascii")
        except UnicodeDecodeError:
            return True
        return False

    @staticmethod
    def _is_unequal_seq_qual(read: Read) -> bool:
        """Determine if read has equal length seq and qual"""
        seq = read.seq.decode("ascii").strip()
        qual = read.qual.decode("ascii").strip()
        if len(seq) != len(qual):
            return True
        return False

    @staticmethod
    def _is_different_header(r1: Read, r2: Read) -> bool:
        """Determine if PE has same header"""
        # Ignore read length because they can be different
        regex = re.compile(r"length=\d+") 
        r1_h1 = re.sub(regex, "", r1.h1.decode("ascii").strip())
        r2_h1 = re.sub(regex, "", r2.h1.decode("ascii").strip())
        if r1_h1 != r2_h1:
            return True
        return False

    @staticmethod
    def _read_to_string(read: Read) -> str:
        return read.h1 + read.seq + read.h2 + read.qual
