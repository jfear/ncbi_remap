import re
from collections import namedtuple
from typing import Union, List, Set
from itertools import combinations

from more_itertools import flatten, windowed

Coordinate = namedtuple("Coordinate", "start end")


class GtfFeature:
    field_separator = "; "
    keyval_separator = " "
    multival_separator = ","
    line_terminator = ";\n"
    val_quote = True

    def __init__(
        self,
        seqid: str = ".",
        source: str = ".",
        featuretype: str = ".",
        start: int = ".",
        stop: int = ".",
        score: int = ".",
        strand: str = ".",
        frame: int = ".",
        attributes: dict = None,
    ):
        self.seqid = seqid
        self.source = source
        self.featuretype = featuretype
        self.start = start
        self.stop = stop
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attributes = attributes or {}

        if self.start != ".":
            self.start = int(self.start)

        if self.stop != ".":
            self.stop = int(self.stop)

        if self.score != ".":
            self.score = int(self.score)

        if self.frame != ".":
            self.frame = int(self.frame)

    def add_annotation(self, attributes: dict):
        for k, v in attributes.items():
            if isinstance(v, str):
                v = [v]
            if k in self.attributes:
                self.attributes[k] |= set(v)
            else:
                self.attributes[k] = set(v)

    @classmethod
    def from_string(cls, gtf_row: str):
        parts = gtf_row.strip().split("\t")  # type: list
        attrs = cls._parse_attributes(parts.pop(-1))
        return cls(*parts, attributes=attrs)

    @staticmethod
    def _parse_attributes(attrs):
        pattern = re.compile(r'(?P<key>[^";]*) "(?P<value>.*?)"')
        matches = re.findall(pattern, attrs)
        return {k.strip(): v.strip() for k, v in matches}

    @property
    def chrom(self):
        return self.seqid

    @chrom.setter
    def chrom(self, value):
        self.seqid = value

    @property
    def end(self):
        return self.stop

    @end.setter
    def end(self, value):
        self.stop = value

    def _attribute_str(self):
        attrs = []
        for k, v in self.attributes.items():
            if isinstance(v, str):
                v = [v]

            if self.val_quote:
                attrs.append(f'{k}{self.keyval_separator}"{self.multival_separator.join(v)}"')
            else:
                attrs.append(f"{k}{self.keyval_separator}{self.multival_separator.join(v)}")
        return self.field_separator.join(attrs)

    def __repr__(self):
        return f"({self.__class__.__name__}: {self.chrom}({self.strand}):{self.start}-{self.stop})"

    def __str__(self):
        return (
            "\t".join(
                [
                    f"{self.chrom}",
                    f"{self.source}",
                    f"{self.featuretype}",
                    f"{self.start}",
                    f"{self.end}",
                    f".",
                    f"{self.strand}",
                    f".",
                    f"{self._attribute_str()}",
                ]
            )
            + self.line_terminator
        )


class Segmentor:
    def __init__(self):
        self.counter = 1
        self.chrom = None
        self.source = "Fear"
        self.featuretype = "segment"
        self.segments = None
        self.features = None

    def run(self, features: List[GtfFeature]) -> List[GtfFeature]:
        self.features = features
        self.chrom = features[0].chrom
        self._build()
        self._annotate()
        return self.segments

    def _build(self):
        """Takes overlapping features and creates non-overlapping segments."""
        self.segments = []
        for segment in self._segment(set([Coordinate(x.start, x.end) for x in self.features])):
            attributes = dict(ID=f"SEG{self.counter:06d}")
            self.segments.append(
                GtfFeature(
                    seqid=self.chrom,
                    source=self.source,
                    featuretype=self.featuretype,
                    start=segment.start,
                    stop=segment.end,
                    attributes=attributes,
                )
            )
            self.counter += 1

    def _segment(self, coordinates: Set[Coordinate]) -> List[Coordinate]:
        if len(coordinates) == 1:
            return coordinates

        coordinate_sorted = coordinate_sort(coordinates)  # type: list
        segments = set()
        for coord1, coord2 in windowed(coordinate_sorted, 2):
            if not coordinate_overlap(coord1, coord2):
                segments.add(coord1)
                segments.add(coord2)
                continue
            elif coord1.start == coord2.start:
                segments.add(Coordinate(coord1.start, coord1.end))
                segments.add(Coordinate(coord1.end + 1, coord2.end))
                continue
            elif coord1.end == coord2.end:
                segments.add(Coordinate(coord1.start, coord2.start - 1))
                segments.add(Coordinate(coord2.start, coord1.end))
                continue
            elif coord1.start <= coord2.start <= coord2.end <= coord1.end:
                segments.add(Coordinate(coord1.start, coord2.start - 1))
                segments.add(Coordinate(coord2.start, coord2.end))
                segments.add(Coordinate(coord2.end + 1, coord1.end))
                continue
            elif coord2.start <= coord1.start <= coord1.end <= coord2.end:
                segments.add(Coordinate(coord2.start, coord1.start - 1))
                segments.add(Coordinate(coord1.start, coord1.end))
                segments.add(Coordinate(coord1.end + 1, coord2.end))
                continue
            else:
                segments.add(Coordinate(coord1.start, coord2.start - 1))
                segments.add(Coordinate(coord2.start, coord1.end))
                segments.add(Coordinate(coord1.end + 1, coord2.end))
                continue

        segments_sorted = coordinate_sort(segments)  # type: list
        if coordinate_sorted != segments_sorted:
            return self._segment(segments)

        return segments_sorted

    def _annotate(self):
        """Combines annotations from overlapping exons."""
        for segment in self.segments:  # type: GtfFeature
            strands = set()
            for feature in self.features:  # type: GtfFeature
                if feature_overlap(segment, feature, stranded=False):
                    strands.add(feature.strand)
                    segment.add_annotation(feature.attributes)

            if len(strands) == 1:
                segment.strand = strands.pop()


class Fuser:
    def __init__(self):
        self.counter = 1
        self.chrom = None
        self.source = "Fear"
        self.featuretype = "fusion"
        self.features = None
        self.fusion = None

    def run(self, features: List[GtfFeature]) -> List[GtfFeature]:
        self.features = features
        self.chrom = features[0].chrom
        self._build()
        self._annotate()
        return [self.fusion]

    def _build(self):
        """Takes overlapping features and fuses them."""
        coords = list(flatten(((x.start, x.end) for x in self.features)))
        start = min(coords)
        end = max(coords)

        attributes = dict(ID=f"FUS{self.counter:06d}")
        self.fusion = GtfFeature(
            seqid=self.chrom,
            source=self.source,
            featuretype=self.featuretype,
            start=start,
            stop=end,
            attributes=attributes,
        )
        self.counter += 1

    def _annotate(self):
        """Combines annotations from overlapping exons."""
        strands = set()

        for feature in self.features:
            strands.add(feature.strand)
            self.fusion.add_annotation(feature.attributes)

        if len(strands) == 1:
            self.fusion.strand = strands.pop()


class FeatureAccumulator:
    """Aggregates features together using two methods.

    Parameters
    ----------
    agg_method : str ["segment", "fusion"]
        If you want features to be fused (merged) or split into
        non-overlapping segments.
    stranded : bool
        If True then strand will be accounted for when running. Note that
        this will assume that features are sorted by chromosome, strand,
        start, end. For example using `gffutils`.

    Example
    -------
    >>> db = gffutils.FeatureDB("your_reference_file.gtf.db")
    >>> exons = db.features_of_type("exon", order_by=["seqid", "strand", "start", "end"])
    >>> accumulator = FeatureAccumulator("segment", stranded=True)
    >>> for exon in exons:
    >>>     if accumulator.add_feature(exon):
    >>>         continue
    >>>     ...

    """

    def __init__(self, agg_method: str = "segment", stranded: bool = True):
        self.agg_method = agg_method
        self.stranded = stranded
        self.chrom = None
        self.strand = None
        self.start = None
        self.end = None
        self.features = None

        if self.agg_method == "segment":
            self.merger = Segmentor()
        if self.agg_method == "fusion":
            self.merger = Fuser()

    def add_feature(self, feature: GtfFeature) -> bool:
        """Adds an interval if it is overlapping."""
        if (self.start is None) & (self.end is None):
            self.reset(feature)
            return True
        elif feature_overlap(self, feature, stranded=self.stranded):
            self.start = min(self.start, feature.start)
            self.end = max(self.end, feature.end)
            self.features.append(feature)
            return True
        return False

    def reset(self, feature: GtfFeature):
        """Resets object to look at next interval group"""
        self.chrom = feature.chrom
        self.strand = feature.strand
        self.start = feature.start
        self.end = feature.end
        self.features = [feature]

    def merge(self):
        return "".join([str(x) for x in self.merger.run(self.features)])


def feature_overlap(feature_a: GtfFeature, feature_b: GtfFeature, stranded: bool = True) -> bool:
    """Tests if two features overlap."""
    if feature_a.chrom != feature_b.chrom:
        return False

    if stranded & (feature_a.strand != feature_b.strand):
        return False

    if (feature_a.start <= feature_b.start < feature_a.end) | (
        feature_b.start <= feature_a.start < feature_b.end
    ):
        return True

    return False


def coordinate_overlap(coord1: Coordinate, coord2: Coordinate) -> bool:
    if (coord1.start <= coord2.start < coord1.end) | (coord2.start <= coord1.start < coord2.end):
        return True
    return False


def coordinate_sort(coordinates: List[Coordinate]) -> List[Coordinate]:
    return sorted(coordinates, key=lambda x: (x.start, x.end))
