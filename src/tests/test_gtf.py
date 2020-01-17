import pytest

from ncbi_remap import gtf


def test_feature_from_string():
    row = 'chr3R	FlyBase	exon	7966463	7967408	3	-	.	gene_id "FBgn0000504"; gene_symbol "dsx"; transcript_id "FBtr0081759"; transcript_symbol "dsx-RA"; '
    feature = gtf.GtfFeature.from_string(row)
    assert feature.chrom == "chr3R"
    assert feature.start == 7966463
    assert feature.attributes["gene_symbol"] == "dsx"
    assert feature.attributes["transcript_symbol"] == "dsx-RA"


def test_segmentation_exact_overlap():
    rows = [
        'chr3R	Test	exon	1	10	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
        'chr3R	Test	exon	1	10	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
        'chr3R	Test	exon	1	10	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
        'chr3R	Test	exon	1	10	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
    ]
    features = [gtf.GtfFeature.from_string(x) for x in rows]
    segmentor = gtf.Segmentor()
    segments = segmentor.run(features)
    assert len(segments) == 1
    assert (segments[0].start, segments[0].stop) == (1, 10)


def test_segmentation_simple_overlap():
    rows = [
        'chr3R	Test	exon	1	10	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
        'chr3R	Test	exon	5	20	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
    ]
    features = [gtf.GtfFeature.from_string(x) for x in rows]
    segmentor = gtf.Segmentor()
    segments = segmentor.run(features)
    assert len(segments) == 3
    for segment in segments:
        assert (segment.start, segment.stop) in [(1, 4), (5, 10), (11, 20)]


def test_segmentation_simple_inset():
    rows = [
        'chr3R	Test	exon	1	10	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
        'chr3R	Test	exon	3	8	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
    ]
    features = [gtf.GtfFeature.from_string(x) for x in rows]
    segmentor = gtf.Segmentor()
    segments = segmentor.run(features)
    assert len(segments) == 3
    for segment in segments:
        assert (segment.start, segment.stop) in [(1, 2), (3, 8), (9, 10)]


def test_segmentation_simple_extension_same_start():
    rows = [
        'chr3R	Test	exon	1	10	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
        'chr3R	Test	exon	1	7	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
    ]

    features = [gtf.GtfFeature.from_string(x) for x in rows]
    segmentor = gtf.Segmentor()
    segments = segmentor.run(features)
    assert len(segments) == 2
    for segment in segments:
        assert (segment.start, segment.stop) in [(1, 7), (8, 10)]


def test_segmentation_simple_extension_same_end():
    rows = [
        'chr3R	Test	exon	1	10	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
        'chr3R	Test	exon	5	10	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
    ]

    features = [gtf.GtfFeature.from_string(x) for x in rows]
    segmentor = gtf.Segmentor()
    segments = segmentor.run(features)
    assert len(segments) == 2
    for segment in segments:
        assert (segment.start, segment.stop) in [(1, 4), (5, 10)]


def test_segmentation_three():
    rows = [
        'chr3R	Test	exon	1	10	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
        'chr3R	Test	exon	5	20	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
        'chr3R	Test	exon	3	8	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
    ]

    features = [gtf.GtfFeature.from_string(x) for x in rows]
    segmentor = gtf.Segmentor()
    segments = segmentor.run(features)
    assert len(segments) == 5
    for segment in segments:
        assert (segment.start, segment.stop) in [(1, 2), (3, 4), (5, 8), (9, 10), (11, 20)]


def test_segment_attribute_aggregation():
    rows = [
        'chr3R	FlyBase	exon	7924323	7925360	3	-	.	gene_id "FBgn0000504"; gene_symbol "dsx"; transcript_id "FBtr0081759"; transcript_symbol "dsx-RA";',
        'chr3R	FlyBase	exon	7924323	7925360	3	-	.	gene_id "FBgn0000504"; gene_symbol "dsx"; transcript_id "FBtr0330073"; transcript_symbol "dsx-RD";',
        'chr3R	FlyBase	exon	7924323	7925360	3	-	.	gene_id "FBgn0000504"; gene_symbol "dsx"; transcript_id "FBtr0330074"; transcript_symbol "dsx-RE";',
        'chr3R	FlyBase	exon	7923468	7924491	15	+	.	gene_id "FBgn0002542"; gene_symbol "lds"; transcript_id "FBtr0081758"; transcript_symbol "lds-RA";',
    ]

    features = [gtf.GtfFeature.from_string(x) for x in rows]
    segmentor = gtf.Segmentor()
    segments = segmentor.run(features)
    for segment in segments:
        if segment.stop < 7924323:
            assert segment.attributes["gene_symbol"] == set(["lds"])

        elif (7924323 <= segment.start) & (segment.end <= 7924491):
            assert segment.attributes["gene_symbol"] == set(["dsx", "lds"])
        else:
            assert segment.attributes["gene_symbol"] == set(["dsx"])


def test_fusion_exact_overlap():
    rows = [
        'chr3R	Test	exon	1	10	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
        'chr3R	Test	exon	1	10	.	-	.	gene_id "Test2"; gene_symbol "test"; ',
        'chr3R	Test	exon	1	10	.	-	.	gene_id "Test3"; gene_symbol "test"; ',
        'chr3R	Test	exon	1	10	.	-	.	gene_id "Test4"; gene_symbol "test"; ',
    ]

    features = [gtf.GtfFeature.from_string(x) for x in rows]
    fuser = gtf.Fuser()
    fusion = fuser.run(features)
    assert (fusion.start, fusion.stop) == (1, 10)
    assert len(fusion.attributes["gene_id"]) == 4


def test_fusion_simple_overlap():
    rows = [
        'chr3R	Test	exon	1	10	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
        'chr3R	Test	exon	5	20	.	-	.	gene_id "Test2"; gene_symbol "test"; ',
    ]

    features = [gtf.GtfFeature.from_string(x) for x in rows]
    fuser = gtf.Fuser()
    fusion = fuser.run(features)
    assert (fusion.start, fusion.stop) == (1, 20)


def test_fusion_simple_inset():
    rows = [
        'chr3R	Test	exon	1	20	.	-	.	gene_id "Test1"; gene_symbol "test"; ',
        'chr3R	Test	exon	5	10	.	-	.	gene_id "Test2"; gene_symbol "test"; ',
    ]

    features = [gtf.GtfFeature.from_string(x) for x in rows]
    fuser = gtf.Fuser()
    fusion = fuser.run(features)
    assert (fusion.start, fusion.stop) == (1, 20)


def test_fusion_attribute_aggregation():
    rows = [
        'chr3R	FlyBase	exon	7924323	7925360	3	-	.	gene_id "FBgn0000504"; gene_symbol "dsx"; transcript_id "FBtr0081759"; transcript_symbol "dsx-RA";',
        'chr3R	FlyBase	exon	7924323	7925360	3	-	.	gene_id "FBgn0000504"; gene_symbol "dsx"; transcript_id "FBtr0330073"; transcript_symbol "dsx-RD";',
        'chr3R	FlyBase	exon	7924323	7925360	3	-	.	gene_id "FBgn0000504"; gene_symbol "dsx"; transcript_id "FBtr0330074"; transcript_symbol "dsx-RE";',
        'chr3R	FlyBase	exon	7923468	7924491	15	+	.	gene_id "FBgn0002542"; gene_symbol "lds"; transcript_id "FBtr0081758"; transcript_symbol "lds-RA";',
    ]

    features = [gtf.GtfFeature.from_string(x) for x in rows]
    fuser = gtf.Fuser()
    fusion = fuser.run(features)
    assert (fusion.start, fusion.stop) == (7923468, 7925360)
    assert fusion.attributes["gene_symbol"] == set(["dsx", "lds"])

