"""Unit tests for GeneFilter."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import pytest

from vartriage.filter.gene_filter import GeneFilter
from vartriage.models.config import GeneFilterConfig
from vartriage.models.variant import (
    AnnotatedVariant,
    FunctionalConsequence,
    Variant,
)


def _make_variant(
    gene_name: Optional[str] = "BRCA1",
    chrom: str = "chr1",
    pos: int = 100,
) -> AnnotatedVariant:
    """Build a minimal AnnotatedVariant for testing."""
    return AnnotatedVariant(
        variant=Variant(
            chrom=chrom,
            pos=pos,
            id=None,
            ref="A",
            alt="T",
            qual=30.0,
            filter_status="PASS",
        ),
        consequence=FunctionalConsequence.MISSENSE,
        gene_name=gene_name,
    )


def _write_gene_file(tmp_path: Path, content: str) -> Path:
    """Write content to a gene list file and return the path."""
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text(content)
    return gene_file


class TestGeneListLoading:
    """Gene list file parsing and validation."""

    def test_loads_mixed_case_symbols_as_uppercase(
        self, tmp_path: Path
    ) -> None:
        gene_file = _write_gene_file(
            tmp_path, "brca1\nTP53\nMlh1\n"
        )
        config = GeneFilterConfig(gene_list_path=gene_file)
        gf = GeneFilter(config)

        assert gf.genes == frozenset({"BRCA1", "TP53", "MLH1"})

    def test_skips_blank_lines(self, tmp_path: Path) -> None:
        gene_file = _write_gene_file(
            tmp_path, "BRCA1\n\n\nTP53\n   \nMLH1\n"
        )
        config = GeneFilterConfig(gene_list_path=gene_file)
        gf = GeneFilter(config)

        assert gf.genes == frozenset({"BRCA1", "TP53", "MLH1"})

    def test_skips_comment_lines(self, tmp_path: Path) -> None:
        gene_file = _write_gene_file(
            tmp_path,
            "# This is a comment\nBRCA1\n# Another comment\nTP53\n",
        )
        config = GeneFilterConfig(gene_list_path=gene_file)
        gf = GeneFilter(config)

        assert gf.genes == frozenset({"BRCA1", "TP53"})

    def test_strips_leading_trailing_whitespace(
        self, tmp_path: Path
    ) -> None:
        gene_file = _write_gene_file(
            tmp_path, "  BRCA1  \n\tTP53\t\n"
        )
        config = GeneFilterConfig(gene_list_path=gene_file)
        gf = GeneFilter(config)

        assert gf.genes == frozenset({"BRCA1", "TP53"})

    def test_raises_file_not_found_on_missing_file(
        self, tmp_path: Path
    ) -> None:
        missing = tmp_path / "nonexistent.txt"
        config = GeneFilterConfig(gene_list_path=missing)

        with pytest.raises(FileNotFoundError, match="not found"):
            GeneFilter(config)

    def test_raises_value_error_on_only_comments_and_blanks(
        self, tmp_path: Path
    ) -> None:
        gene_file = _write_gene_file(
            tmp_path, "# comment\n\n   \n# another\n"
        )
        config = GeneFilterConfig(gene_list_path=gene_file)

        with pytest.raises(ValueError, match="no valid gene"):
            GeneFilter(config)


class TestGeneFilterInclusion:
    """Filter correctly includes matching variants."""

    def test_variant_with_matching_gene_passes(
        self, tmp_path: Path
    ) -> None:
        gene_file = _write_gene_file(tmp_path, "BRCA1\nTP53\n")
        config = GeneFilterConfig(gene_list_path=gene_file)
        gf = GeneFilter(config)

        variant = _make_variant(gene_name="BRCA1")
        result = list(gf.apply(iter([variant])))

        assert result == [variant]

    def test_case_insensitive_matching(
        self, tmp_path: Path
    ) -> None:
        gene_file = _write_gene_file(tmp_path, "BRCA1\n")
        config = GeneFilterConfig(gene_list_path=gene_file)
        gf = GeneFilter(config)

        variant = _make_variant(gene_name="brca1")
        result = list(gf.apply(iter([variant])))

        assert result == [variant]

    def test_mixed_case_variant_matches_uppercase_gene_set(
        self, tmp_path: Path
    ) -> None:
        gene_file = _write_gene_file(tmp_path, "tp53\n")
        config = GeneFilterConfig(gene_list_path=gene_file)
        gf = GeneFilter(config)

        variant = _make_variant(gene_name="Tp53")
        result = list(gf.apply(iter([variant])))

        assert result == [variant]


class TestGeneFilterExclusion:
    """Filter correctly excludes non-matching variants."""

    def test_variant_with_non_matching_gene_excluded(
        self, tmp_path: Path
    ) -> None:
        gene_file = _write_gene_file(tmp_path, "BRCA1\nTP53\n")
        config = GeneFilterConfig(gene_list_path=gene_file)
        gf = GeneFilter(config)

        variant = _make_variant(gene_name="EGFR")
        result = list(gf.apply(iter([variant])))

        assert result == []

    def test_intergenic_variant_excluded(
        self, tmp_path: Path
    ) -> None:
        gene_file = _write_gene_file(tmp_path, "BRCA1\nTP53\n")
        config = GeneFilterConfig(gene_list_path=gene_file)
        gf = GeneFilter(config)

        variant = _make_variant(gene_name=None)
        result = list(gf.apply(iter([variant])))

        assert result == []


class TestUnmatchedGeneDetection:
    """Unmatched gene tracking after stream exhaustion."""

    def test_unmatched_genes_reported_after_exhaustion(
        self, tmp_path: Path
    ) -> None:
        gene_file = _write_gene_file(
            tmp_path, "BRCA1\nTP53\nMLH1\n"
        )
        config = GeneFilterConfig(gene_list_path=gene_file)
        gf = GeneFilter(config)

        variants = [
            _make_variant(gene_name="BRCA1"),
            _make_variant(gene_name="TP53"),
        ]
        list(gf.apply(iter(variants)))

        assert gf.unmatched_genes == frozenset({"MLH1"})

    def test_no_unmatched_when_all_genes_hit(
        self, tmp_path: Path
    ) -> None:
        gene_file = _write_gene_file(tmp_path, "BRCA1\nTP53\n")
        config = GeneFilterConfig(gene_list_path=gene_file)
        gf = GeneFilter(config)

        variants = [
            _make_variant(gene_name="BRCA1"),
            _make_variant(gene_name="TP53"),
        ]
        list(gf.apply(iter(variants)))

        assert gf.unmatched_genes == frozenset()

    def test_unmatched_genes_logged_as_warning(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        gene_file = _write_gene_file(
            tmp_path, "BRCA1\nTP53\nMLH1\n"
        )
        config = GeneFilterConfig(gene_list_path=gene_file)
        gf = GeneFilter(config)

        variants = [_make_variant(gene_name="BRCA1")]

        with caplog.at_level(logging.WARNING):
            list(gf.apply(iter(variants)))

        assert "MLH1" in caplog.text
        assert "TP53" in caplog.text
