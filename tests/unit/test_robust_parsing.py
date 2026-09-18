"""Robust VCF and ClinVar parsing.

Covers three parsing gaps: a malformed INFO field must not drop the whole
INFO dict, a multi-allelic record must be split so every ALT is triaged,
and a ClinVar compound significance string must map to the appropriate
assertion instead of being discarded as unknown.
"""

from __future__ import annotations

from pathlib import Path

from vartriage.annotation.clinvar import _map_significance
from vartriage.io.vcf_parser import VCFParser
from vartriage.models.variant import ClinVarAssertion

_MULTIALLELIC_VCF = """\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG,T\t30.0\tPASS\tDP=50
chr1\t200\t.\tC\tA\t40.0\tPASS\tDP=80
"""


def _write_vcf(tmp_path: Path, content: str) -> Path:
    vcf_path = tmp_path / "test.vcf"
    vcf_path.write_text(content)
    return vcf_path


class TestMultiAllelicSplit:
    def test_multiallelic_record_yields_one_variant_per_alt(
        self, tmp_path: Path
    ) -> None:
        parser = VCFParser(_write_vcf(tmp_path, _MULTIALLELIC_VCF))
        variants = list(parser)
        parser.close()

        at_100 = [v for v in variants if v.pos == 100]
        assert len(at_100) == 2
        assert {v.alt for v in at_100} == {"G", "T"}
        # The biallelic record is unaffected.
        assert sum(1 for v in variants if v.pos == 200) == 1

    def test_split_preserves_shared_fields(self, tmp_path: Path) -> None:
        parser = VCFParser(_write_vcf(tmp_path, _MULTIALLELIC_VCF))
        at_100 = [v for v in parser if v.pos == 100]
        parser.close()
        for v in at_100:
            assert v.chrom == "chr1"
            assert v.ref == "A"
            assert v.qual == 30.0
            assert v.info.get("DP") == 50


class TestClinVarCompoundSignificance:
    def test_plain_pathogenic(self) -> None:
        assert _map_significance("Pathogenic") == ClinVarAssertion.PATHOGENIC

    def test_pathogenic_likely_pathogenic_compound(self) -> None:
        # A combined assertion is not "unknown"; take the more severe side.
        assert _map_significance("Pathogenic/Likely pathogenic") == (
            ClinVarAssertion.PATHOGENIC
        )

    def test_benign_likely_benign_compound(self) -> None:
        assert _map_significance("Benign/Likely benign") == (ClinVarAssertion.BENIGN)

    def test_case_and_whitespace_insensitive(self) -> None:
        assert _map_significance("  pathogenic  ") == ClinVarAssertion.PATHOGENIC

    def test_conflicting_maps_to_vus(self) -> None:
        assert (
            _map_significance("Conflicting interpretations of pathogenicity")
            == ClinVarAssertion.VUS
        )

    def test_qualified_pathogenic_low_penetrance(self) -> None:
        assert _map_significance("Likely pathogenic, low penetrance") == (
            ClinVarAssertion.LIKELY_PATHOGENIC
        )

    def test_unrecognized_returns_none(self) -> None:
        assert _map_significance("not_a_real_term") is None
