"""Unit tests for the ClinGen LoF-mechanism gene-set loader."""

from __future__ import annotations

from pathlib import Path

from vartriage.knowledge.lof_mechanism import LofMechanismGeneSet


def _write(tmp_path: Path) -> Path:
    content = (
        "# comment line, ignored\n"
        "# another comment\n"
        "gene_symbol\thi_score\n"
        "BRCA2\t3\n"
        "CFTR\t30\n"
        "MSH2\t3\n"
        "PAH\t30\n"
    )
    tsv = tmp_path / "lof_mechanism_genes.tsv"
    tsv.write_text(content)
    return tsv


def test_loads_dominant_and_recessive_lof_genes(tmp_path: Path) -> None:
    gene_set = LofMechanismGeneSet(_write(tmp_path))
    assert "BRCA2" in gene_set  # HI=3 dominant haploinsufficiency
    assert "CFTR" in gene_set  # HI=30 recessive LoF phenotype
    assert "MSH2" in gene_set
    assert "PAH" in gene_set
    assert gene_set.gene_count == 4


def test_absent_gene_is_not_a_member(tmp_path: Path) -> None:
    gene_set = LofMechanismGeneSet(_write(tmp_path))
    # GAA is not ClinGen dosage-curated: it must not appear, so PVS1 fails closed.
    assert "GAA" not in gene_set
    assert "NOVELGENE" not in gene_set


def test_missing_file_loads_empty(tmp_path: Path) -> None:
    gene_set = LofMechanismGeneSet(tmp_path / "missing.tsv")
    assert gene_set.gene_count == 0
    assert "BRCA2" not in gene_set
