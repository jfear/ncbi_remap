from textwrap import dedent
from pathlib import Path

import pytest

from ncbi_remap.fastq import Fastq, UnequalNumberReadsException, MixedUpReadsException

HERE = Path(__file__).parents[1].absolute()

FASTQ_STRING = dedent(
    """\
    @SRR823324.sra.1 HWI-ST1024_0265:8:1101:1417:2091 length=29
    TTTCTTAAAAATAGGTAATCTAATTGAAA
    +SRR823324.sra.1 HWI-ST1024_0265:8:1101:1417:2091 length=29
    DDDDDDDDDEIEEECF4<CEFF>EF?CFF
    @SRR823324.sra.2 HWI-ST1024_0265:8:1101:3024:2199 length=26
    TGTAATTGTTTAAGCCAACTATTTGT
    +SRR823324.sra.2 HWI-ST1024_0265:8:1101:3024:2199 length=26
    FFFDHHHGHJJEHFIHJJJEGIJJGI
"""
)

POORLY_ENCODED_FASTQ_STRING = (
    dedent(
        """\
    @SRR823324.sra.1 HWI-ST1024_0265:8:1101:1417:2091 length=29
    TTTCTTAAAAATAGGTAATCTAATTGAAA
    +SRR823324.sra.1 HWI-ST1024_0265:8:1101:1417:2091 length=29
    DDDDDDDDDEIEEECF4<CEFF>EF?CFF
    @SRR823324.sra.2 HWI-ST1024_0265:8:1101:3024:2199 length=26
    TGTAATTGTTTAAGCCAACTATTTGT
    +SRR823324.sra.2 HWI-ST1024_0265:8:1101:3024:2199 length=26
"""
    ).encode("ascii")
    + "FFFDHHHGHJJEHFIHJJJEGIJJGI\n".encode("utf-16")
)


INCOMPLETE_FASTQ_STRING = dedent(
    """\
    @SRR823324.sra.1 HWI-ST1024_0265:8:1101:1417:2091 length=29
    TTTCTTAAAAATAGGTAATCTAATTGAAA
    +SRR823324.sra.1 HWI-ST1024_0265:8:1101:1417:2091 length=29
    DDDDDDDDDEIEEECF4<CEFF>EF?CFF
    @SRR823324.sra.2 HWI-ST1024_0265:8:1101:3024:2199 length=26
    TGTAATTGTTTAAGCCAACTATTTGT
    +SRR823324.sra.2 HWI-ST1024_0265:8:1101:3024:2199 length=26
"""
)

LONGER_SEQ_STRING = dedent(
    """\
    @SRR823324.sra.1 HWI-ST1024_0265:8:1101:1417:2091 length=29
    TTTCTTAAAAATAGGTAATCTAATTGAAA
    +SRR823324.sra.1 HWI-ST1024_0265:8:1101:1417:2091 length=29
    DDDDDDDDDEIEEECF4<CEFF>EF?CFF
    @SRR823324.sra.2 HWI-ST1024_0265:8:1101:3024:2199 length=26
    TGTAATTGTTTAAGCCAACTATTTGT
    +SRR823324.sra.2 HWI-ST1024_0265:8:1101:3024:2199 length=26
    FFFDHHHGHJJEHFIHJJJEGIJJ
"""
)

LONGER_QUAL_STRING = dedent(
    """\
    @SRR823324.sra.1 HWI-ST1024_0265:8:1101:1417:2091 length=29
    TTTCTTAAAAATAGGTAATCTAATTGAAA
    +SRR823324.sra.1 HWI-ST1024_0265:8:1101:1417:2091 length=29
    DDDDDDDDDEIEEECF4<CEFF>EF?CFF
    @SRR823324.sra.2 HWI-ST1024_0265:8:1101:3024:2199 length=26
    TGTAATTGTTTAAGCCAACTATTTGT
    +SRR823324.sra.2 HWI-ST1024_0265:8:1101:3024:2199 length=26
    FFFDHHHGHJJEHFIHJJJEGIJJGI!!!!
"""
)

EMPTY_FASTQ = ""

ABI_FASTQ = dedent(
    """\
    @SRR######.1 solid0527_####### length=50
    T120202210232000020002.00301000012.100...00........
    +SRR######.1 solid0527_####### length=50
    !/<%2/:%*)-%%0'--'')/.!%('1'%),+/%!&',!!!'+!!!!!!!!
"""
)


@pytest.mark.parametrize(
    "fastq,num_reads",
    [
        ((HERE / "test_data/test.fastq").as_posix(), 2),
        ((HERE / "test_data/test.fastq.gz").as_posix(), 2),
        (FASTQ_STRING, 2),
    ],
)
def test_read_fastq_file(fastq, num_reads):
    fq = Fastq(fastq)
    with fq.open_fastq() as fh:
        reads = list(fq.iter_reads(fh))
        assert len(reads) == num_reads
        assert len(reads[0]) == 4
        assert len(reads[1]) == 4


def test_non_ascii_encoding():
    fq = Fastq(POORLY_ENCODED_FASTQ_STRING)
    reads = list(fq.process())
    assert len(reads) == 1
    assert fq.bad_ecoding == 1
    assert fq.incomplete_read == 0


def test_incomplete_read():
    fq = Fastq(INCOMPLETE_FASTQ_STRING)
    reads = list(fq.process())
    assert len(reads) == 1
    assert fq.incomplete_read == 1
    assert fq.bad_ecoding == 0


@pytest.mark.parametrize("fastq", [LONGER_QUAL_STRING, LONGER_SEQ_STRING])
def test_remove_reads_with_unequal_seq_and_qual(fastq):
    fq = Fastq(fastq)
    list(fq.process())
    assert fq.unequal_len == 1
    assert fq.libsize == 1


def test_single_end():
    r1 = (HERE / "test_data/sample_1.fastq").as_posix()
    fq = Fastq(r1)
    reads = list(fq.process())
    assert len(reads) == 10
    assert fq.libsize == 10
    assert fq.bad_ecoding == 0
    assert fq.incomplete_read == 0
    assert "SE" in fq.flags


def test_pair_end():
    r1 = (HERE / "test_data/sample_1.fastq").as_posix()
    r2 = (HERE / "test_data/sample_2.fastq").as_posix()
    fq = Fastq(r1, r2)
    reads = list(fq.process())
    assert len(reads) == 10
    assert fq.libsize == 10
    assert fq.bad_ecoding == 0
    assert fq.incomplete_read == 0
    assert "PE" in fq.flags


def test_read2_too_short():
    r1 = (HERE / "test_data/sample_1.fastq").as_posix()
    r2 = (HERE / "test_data/sample_short_2.fastq").as_posix()
    fq = Fastq(r1, r2)
    with pytest.raises(UnequalNumberReadsException):
        reads = list(fq.process())

    assert "keep_R1" in fq.flags


def test_rerun_keep_r1():
    r1 = (HERE / "test_data/sample_1.fastq").as_posix()
    r2 = (HERE / "test_data/sample_short_2.fastq").as_posix()
    fq = Fastq(r1, r2)
    fq.flags.add("keep_R1")
    reads = list(fq.process())
    assert "keep_R1" in fq.flags
    assert len(reads) == 10
    assert fq.libsize == 10
    assert fq.bad_ecoding == 0
    assert fq.incomplete_read == 0


def test_read1_too_short():
    r1 = (HERE / "test_data/sample_short_1.fastq").as_posix()
    r2 = (HERE / "test_data/sample_2.fastq").as_posix()
    fq = Fastq(r1, r2)
    with pytest.raises(UnequalNumberReadsException):
        reads = list(fq.process())

    assert "keep_R2" in fq.flags

def test_rerun_keep_r2():
    r1 = (HERE / "test_data/sample_short_1.fastq").as_posix()
    r2 = (HERE / "test_data/sample_2.fastq").as_posix()
    fq = Fastq(r1, r2)
    fq.flags.add("keep_R2")
    reads = list(fq.process())
    assert "keep_R2" in fq.flags
    assert len(reads) == 10
    assert fq.libsize == 10
    assert fq.bad_ecoding == 0
    assert fq.incomplete_read == 0

def test_read1_mixed():
    with pytest.raises(MixedUpReadsException):
        r1 = (HERE / "test_data/sample_mixed_order_1.fastq").as_posix()
        r2 = (HERE / "test_data/sample_2.fastq").as_posix()
        fq = Fastq(r1, r2)
        list(fq.process())
    assert "keep_R1" in fq.flags

def test_read2_mixed():
    with pytest.raises(MixedUpReadsException):
        r1 = (HERE / "test_data/sample_1.fastq").as_posix()
        r2 = (HERE / "test_data/sample_mixed_order_2.fastq").as_posix()
        fq = Fastq(r1, r2)
        list(fq.process())
    assert "keep_R1" in fq.flags

