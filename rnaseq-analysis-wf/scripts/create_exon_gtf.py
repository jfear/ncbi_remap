"""Builds segment GTF

This script segments exons based on overlap and strand. It does the
following:

1. Sorts exons by chrom, strand, start, end
2. Checks for overlapping exons
3. Segments those exons
4. Saves a GTF formatted file.

"""
import os
from collections import defaultdict
from typing import Generator

from more_itertools import flatten, windowed
import gffutils


def main():
    db = gffutils.FeatureDB(snakemake.input[0])
    exons = db.features_of_type("exon", order_by=["seqid", "strand", "start", "end"])
    with open(snakemake.output[0], "w") as file_out:
        for i, exon in enumerate(exons):
            if i == 0:
                overlap = FeatureAccumulator(exon)
            elif overlap.add_interval(exon):
                continue
            else:
                file_out.write("".join(overlap.export_segments()))
                overlap.reset(exon)
        file_out.write("".join(overlap.export_segments()))


def feature_overlap(feature_a: gffutils.Feature, feature_b: gffutils.Feature) -> bool:
    """Tests if two features overlap."""
    if (
        (feature_a.chrom == feature_b.chrom)
        & (feature_a.strand == feature_b.strand)
        & (
            (feature_a.start <= feature_b.start <= feature_a.end)
            | (feature_b.start <= feature_a.start <= feature_b.end)
        )
    ):
        return True
    return False


class Segment:
    field_separator = "; "
    keyval_separator = " "
    multival_separator = ","
    line_terminator = ";\n"
    val_quote = True

    def __init__(self, segment_number, chrom, start, end, strand):
        self.name = f"SEG{segment_number:07d}"
        self.chrom = chrom
        self.start = start
        self.end = end or start  # If end is None then set it to start
        self.strand = strand
        self.attributes = defaultdict(set)
        self.attributes["segment_id"].add(self.name)

    def add_annotation(self, attributes: dict):
        for k, v in attributes.items():
            self.attributes[k] |= set(v)

    def _attribute_str(self):
        attrs = []
        for k, v in self.attributes.items():
            if self.val_quote:
                attrs.append(f'{k}{self.keyval_separator}"{self.multival_separator.join(v)}"')
            else:
                attrs.append(f"{k}{self.keyval_separator}{self.multival_separator.join(v)}")
        return self.field_separator.join(attrs)

    def __str__(self):
        attrs = self._attribute_str()
        return (
            "\t".join(
                [
                    f"{self.chrom}",
                    "Fear",
                    "segment",
                    f"{self.start}",
                    f"{self.end}",
                    f".",
                    f"{self.strand}",
                    f".",
                    f"{attrs}",
                ]
            )
            + self.line_terminator
        )


class FeatureAccumulator:
    def __init__(self, feature: gffutils.Feature):
        self.chrom = feature.chrom
        self.strand = feature.strand
        self.start = feature.start
        self.end = feature.end
        self.features = [feature]
        self.ranges = [(feature.start, feature.end)]
        self.segments = []
        self.counter = 1

    def add_interval(self, feature: gffutils.Feature) -> bool:
        """Adds an interval if it is overlapping."""
        if feature_overlap(self, feature):
            self.end = max(self.end, feature.end)
            self.features.append(feature)
            self.ranges.append((feature.start, feature.end))
            return True
        return False

    def reset(self, feature: gffutils.Feature):
        """Resets object to look at next interval group"""
        self.chrom = feature.chrom
        self.strand = feature.strand
        self.start = feature.start
        self.end = feature.end
        self.features = [feature]
        self.segments = []
        self.ranges = [(feature.start, feature.end)]

    def _build_segments(self):
        """Takes overlapping intervals and creates non-overlapping segments."""
        windows = list(windowed(sorted(set(flatten(self.ranges))), 2))
        for i, (start, end) in enumerate(windows):
            if end is None:
                end = start
            elif i < len(windows) - 1:
                # Tweak end so segments don't overlap in GTF format
                end -= 1
            self.segments.append(Segment(self.counter, self.chrom, start, end, self.strand))
            self.counter += 1

    def _annotate_segments(self):
        """Combines annotations from overlapping exons."""
        for segment in self.segments:
            for feature in self.features:
                if feature_overlap(segment, feature):
                    segment.add_annotation(feature.attributes)

    def export_segments(self) -> Generator[str, None, None]:
        self._build_segments()
        self._annotate_segments()
        for segment in self.segments:
            yield str(segment)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="rnaseq-analysis-wf",
            input="../../lcdb-references/dmel/r6-11/gtf/dmel_r6-11.gtf.db",
        )

    main()
