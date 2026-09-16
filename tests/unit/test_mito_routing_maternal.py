"""Mitochondrial routing and maternal-inheritance no-call handling.

Routing decides whether a variant goes to the mitochondrial classifier or
the nuclear one, so a mis-route sends mtDNA through the wrong biology. And a
maternal genotype that is a no-call must not be read as "absent in mother",
which would manufacture a de novo call.
"""

from __future__ import annotations

from vartriage.mito.genetic_code import is_mitochondrial
from vartriage.mito.maternal import check_maternal_inheritance
from vartriage.models.variant import Variant


class TestIsMitochondrial:
    def test_recognizes_common_names(self) -> None:
        assert is_mitochondrial("chrM")
        assert is_mitochondrial("MT")
        assert is_mitochondrial("chrMT")
        assert is_mitochondrial("M")

    def test_recognizes_refseq_accession(self) -> None:
        assert is_mitochondrial("NC_012920.1")

    def test_rejects_nuclear_contigs(self) -> None:
        assert not is_mitochondrial("chr1")
        assert not is_mitochondrial("1")
        assert not is_mitochondrial("chrX")
        assert not is_mitochondrial("chr3")

    def test_rejects_a_contig_that_reduces_to_m_under_char_stripping(self) -> None:
        # A defensive case: the old str.lstrip("CHR") logic could reduce an
        # unrelated name to "M"/"MT"; a real routing check must not.
        assert not is_mitochondrial("CHACM")


def _sample(gt: tuple) -> dict:  # noqa: ANN001
    return {"GT": list(gt)}


def _variant() -> Variant:
    return Variant(
        chrom="chrM",
        pos=3243,
        id=None,
        ref="A",
        alt="G",
        qual=None,
        filter_status="PASS",
        info={
            "_pysam_samples": {
                "proband": _sample((1, 1)),
                "mother": _sample((None, None)),  # maternal no-call
                "father": _sample((0, 0)),
            }
        },
    )


class TestMaternalNoCall:
    def test_maternal_no_call_is_unknown_not_de_novo(self) -> None:
        result = check_maternal_inheritance(
            _variant(),
            proband_name="proband",
            mother_name="mother",
            father_name="father",
        )
        assert result is not None
        assert result.status == "unknown"
        assert result.status != "de_novo"

    def test_confirmed_absent_in_both_is_de_novo(self) -> None:
        v = _variant()
        v.info["_pysam_samples"]["mother"] = _sample((0, 0))
        result = check_maternal_inheritance(
            v,
            proband_name="proband",
            mother_name="mother",
            father_name="father",
        )
        assert result is not None
        assert result.status == "de_novo"
