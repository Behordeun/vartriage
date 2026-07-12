"""Trio-based Mendelian inheritance pattern classification.

Replaces SampleExtractor when trio mode is active. Extracts genotype
data for proband-mother-father, classifies each variant into one or
more inheritance patterns, and emits enriched Variant objects.
"""

from __future__ import annotations

import re
from typing import Iterator

from vartriage.models.config import InheritanceConfig
from vartriage.models.variant import Variant


class InheritanceFilter:
    """Streaming filter for trio inheritance pattern classification.

    Extracts trio genotypes from multi-sample VCF variants, classifies
    each into Mendelian inheritance patterns (de_novo, dominant,
    recessive, compound_het, x_linked), and yields enriched Variant
    objects compatible with downstream pipeline stages.

    Parameters
    ----------
    config : InheritanceConfig
        Trio sample names and patterns to evaluate.
    sample_names : list[str]
        Sample names from the VCF header.

    Raises
    ------
    ValueError
        If any trio sample name is not found in sample_names.
    """

    def __init__(
        self,
        config: InheritanceConfig,
        sample_names: list[str],
    ) -> None:
        for name, role in [
            (config.proband, "proband"),
            (config.mother, "mother"),
            (config.father, "father"),
        ]:
            if name not in sample_names:
                raise ValueError(
                    f"{role.capitalize()} sample '{name}' not found "
                    f"in VCF. Available samples: {sample_names}"
                )

        self._proband = config.proband
        self._mother = config.mother
        self._father = config.father
        self._patterns = config.patterns

    def apply(self, variants: Iterator[Variant]) -> Iterator[Variant]:
        """Classify variants by inheritance pattern and yield results.

        Evaluates non-compound patterns per-variant in a streaming
        pass. When compound_het is active, buffers variants by gene
        annotation and flushes on gene boundary change or stream end.

        Parameters
        ----------
        variants : Iterator[Variant]
            Input stream with _pysam_samples in info dict.

        Yields
        ------
        Variant
            Enriched variants with inheritance_pattern, sample_gt,
            sample_name, and optionally sample_gq in info.
        """
        compound_het_active = "compound_het" in self._patterns
        gene_buffer: list[tuple[Variant, str, str, str, list[str]]] = []
        current_gene: str | None = None

        for variant in variants:
            proband_gt = self._extract_genotype(variant, self._proband)
            if proband_gt is None or proband_gt == "./.":
                continue

            if not self._has_alt_allele(proband_gt):
                continue

            mother_gt = self._extract_genotype(variant, self._mother)
            father_gt = self._extract_genotype(variant, self._father)
            if mother_gt is None:
                mother_gt = "./."
            if father_gt is None:
                father_gt = "./."

            patterns: list[str] = []
            if "de_novo" in self._patterns:
                if self._classify_de_novo(proband_gt, mother_gt, father_gt):
                    patterns.append("de_novo")
            if "dominant" in self._patterns:
                if self._classify_dominant(proband_gt, mother_gt, father_gt):
                    patterns.append("dominant")
            if "recessive" in self._patterns:
                if self._classify_recessive(proband_gt, mother_gt, father_gt):
                    patterns.append("recessive")
            if "x_linked" in self._patterns:
                if self._classify_x_linked(proband_gt, mother_gt, variant.chrom):
                    patterns.append("x_linked")

            if compound_het_active:
                gene = variant.info.get("gene")
                if gene is None:
                    yield self._build_output(variant, proband_gt, patterns)
                    continue

                # Gene boundary flush: assumes variants arrive grouped
                # by gene (standard for coordinate-sorted VCFs where
                # same-gene variants are contiguous). Using a single
                # buffer rather than per-gene dicts keeps memory bounded
                # for large inputs.
                if current_gene is not None and gene != current_gene:
                    yield from self._flush_gene_buffer(gene_buffer)
                    gene_buffer = []

                current_gene = gene
                gene_buffer.append(
                    (
                        variant,
                        proband_gt,
                        mother_gt,
                        father_gt,
                        patterns,
                    )
                )
            else:
                yield self._build_output(variant, proband_gt, patterns)

        if compound_het_active and gene_buffer:
            yield from self._flush_gene_buffer(gene_buffer)

    def _extract_genotype(self, variant: Variant, sample_name: str) -> str | None:
        """Pull GT from _pysam_samples and format as VCF string.

        Parameters
        ----------
        variant : Variant
            Variant with _pysam_samples in info.
        sample_name : str
            Sample to extract genotype for.

        Returns
        -------
        str or None
            VCF-style GT string (e.g. "0/1") or None if unavailable.
        """
        sample_data = variant.info.get("_pysam_samples", {})
        sample_entry = sample_data.get(sample_name, {})
        gt = sample_entry.get("GT")

        if gt is None:
            return None

        return self._format_gt(gt)

    @staticmethod
    def _format_gt(gt: tuple[int | None, ...]) -> str:
        """Convert pysam GT tuple to VCF-style string.

        Parameters
        ----------
        gt : tuple
            Allele tuple from pysam, e.g. (0, 1), (None, None).

        Returns
        -------
        str
            VCF-style genotype, e.g. "0/1", "./.".
        """
        alleles = []
        for a in gt:
            alleles.append(str(a) if a is not None else ".")
        return "/".join(alleles)

    @staticmethod
    def _has_alt_allele(gt_str: str) -> bool:
        """Check if genotype contains at least one alt allele (>= 1).

        Parameters
        ----------
        gt_str : str
            VCF-style genotype string.

        Returns
        -------
        bool
        """
        alleles = re.split(r"[/|]", gt_str)
        for allele in alleles:
            if allele == ".":
                continue
            try:
                if int(allele) >= 1:
                    return True
            except ValueError:
                continue
        return False

    @staticmethod
    def _is_het(gt_str: str) -> bool:
        """True if genotype has exactly one ref and one alt allele.

        Parameters
        ----------
        gt_str : str
            VCF-style genotype string.

        Returns
        -------
        bool
        """
        alleles = re.split(r"[/|]", gt_str)
        numeric = []
        for a in alleles:
            if a == ".":
                return False
            try:
                numeric.append(int(a))
            except ValueError:
                return False

        if len(numeric) != 2:
            return False

        has_ref = any(a == 0 for a in numeric)
        has_alt = any(a >= 1 for a in numeric)
        return has_ref and has_alt

    @staticmethod
    def _is_hom_alt(gt_str: str) -> bool:
        """True if all alleles are alternate (>= 1).

        Parameters
        ----------
        gt_str : str
            VCF-style genotype string.

        Returns
        -------
        bool
        """
        alleles = re.split(r"[/|]", gt_str)
        for a in alleles:
            if a == ".":
                return False
            try:
                if int(a) < 1:
                    return False
            except ValueError:
                return False
        return len(alleles) > 0

    @staticmethod
    def _is_hom_ref(gt_str: str) -> bool:
        """True if all alleles are reference (== 0).

        Parameters
        ----------
        gt_str : str
            VCF-style genotype string.

        Returns
        -------
        bool
        """
        alleles = re.split(r"[/|]", gt_str)
        for a in alleles:
            if a == ".":
                return False
            try:
                if int(a) != 0:
                    return False
            except ValueError:
                return False
        return len(alleles) > 0

    def _classify_de_novo(
        self,
        proband_gt: str,
        mother_gt: str,
        father_gt: str,
    ) -> bool:
        """Classify de novo: proband has alt, both parents hom-ref.

        Returns True iff proband carries at least one alt allele AND
        both mother and father are homozygous-reference.
        """
        if not self._has_alt_allele(proband_gt):
            return False
        if not self._is_hom_ref(mother_gt):
            return False
        if not self._is_hom_ref(father_gt):
            return False
        return True

    def _classify_dominant(
        self,
        proband_gt: str,
        mother_gt: str,
        father_gt: str,
    ) -> bool:
        """Classify dominant: proband het, exactly one parent het.

        Returns True iff proband is heterozygous AND exactly one of
        mother or father is heterozygous.
        """
        if not self._is_het(proband_gt):
            return False

        mother_het = self._is_het(mother_gt)
        father_het = self._is_het(father_gt)

        return (mother_het and not father_het) or (father_het and not mother_het)

    def _classify_recessive(
        self,
        proband_gt: str,
        mother_gt: str,
        father_gt: str,
    ) -> bool:
        """Classify recessive: proband hom-alt, both parents het.

        Returns True iff proband is homozygous-alternate AND both
        mother and father are heterozygous.
        """
        if not self._is_hom_alt(proband_gt):
            return False
        if not self._is_het(mother_gt):
            return False
        if not self._is_het(father_gt):
            return False
        return True

    def _classify_x_linked(
        self,
        proband_gt: str,
        mother_gt: str,
        chrom: str,
    ) -> bool:
        """Classify X-linked: chrX, proband has alt, mother is het.

        Returns True iff chromosome is X AND proband carries at least
        one alt allele AND the mother is heterozygous.
        """
        if chrom not in ("chrX", "X"):
            return False
        if not self._has_alt_allele(proband_gt):
            return False
        if not self._is_het(mother_gt):
            return False
        return True

    def _evaluate_compound_het(
        self,
        buffer: list[tuple[Variant, str, str, str, list[str]]],
    ) -> set[int]:
        """Find indices of variants forming compound het trans pairs.

        A trans pair requires at least one variant with the alt allele
        from the mother (mother het, father hom-ref) AND at least one
        with the alt from the father (father het, mother hom-ref).

        Parameters
        ----------
        buffer : list
            Tuples of (variant, proband_gt, mother_gt, father_gt,
            patterns) for a single gene group.

        Returns
        -------
        set[int]
            Indices of variants in trans pairs.
        """
        het_indices = [
            i for i, (_, gt, _, _, _) in enumerate(buffer) if self._is_het(gt)
        ]

        if len(het_indices) < 2:
            return set()

        maternal_indices: list[int] = []
        paternal_indices: list[int] = []

        for i in het_indices:
            _, _, mother_gt, father_gt, _ = buffer[i]
            if self._has_alt_allele(mother_gt) and self._is_hom_ref(father_gt):
                maternal_indices.append(i)
            if self._has_alt_allele(father_gt) and self._is_hom_ref(mother_gt):
                paternal_indices.append(i)

        if maternal_indices and paternal_indices:
            return set(maternal_indices + paternal_indices)

        return set()

    def _flush_gene_buffer(
        self,
        buffer: list[tuple[Variant, str, str, str, list[str]]],
    ) -> Iterator[Variant]:
        """Evaluate compound het and yield all buffered variants.

        Adds "compound_het" to the patterns list for variants that
        form trans pairs, then builds output for each variant.

        Parameters
        ----------
        buffer : list
            Gene-grouped variant tuples.

        Yields
        ------
        Variant
            Enriched variants with inheritance_pattern assigned.
        """
        compound_het_indices = self._evaluate_compound_het(buffer)

        for i, (
            variant,
            proband_gt,
            _,
            _,
            patterns,
        ) in enumerate(buffer):
            if i in compound_het_indices:
                patterns = patterns + ["compound_het"]
            yield self._build_output(variant, proband_gt, patterns)

    def _build_output(
        self,
        variant: Variant,
        proband_gt: str,
        patterns: list[str],
    ) -> Variant:
        """Build output Variant with enriched info dict.

        Strips _pysam_samples and attaches sample_gt, sample_name,
        sample_gq, and inheritance_pattern.

        Parameters
        ----------
        variant : Variant
            Original variant.
        proband_gt : str
            Formatted proband genotype string.
        patterns : list[str]
            Matched inheritance patterns.

        Returns
        -------
        Variant
            New Variant with enriched info dict.
        """
        sample_data = variant.info.get("_pysam_samples", {})
        proband_entry = sample_data.get(self._proband, {})
        gq = proband_entry.get("GQ")

        new_info = dict(variant.info)
        new_info.pop("_pysam_samples", None)
        new_info["sample_gt"] = proband_gt
        new_info["sample_name"] = self._proband
        if gq is not None:
            new_info["sample_gq"] = gq
        new_info["inheritance_pattern"] = patterns

        return Variant(
            chrom=variant.chrom,
            pos=variant.pos,
            id=variant.id,
            ref=variant.ref,
            alt=variant.alt,
            qual=variant.qual,
            filter_status=variant.filter_status,
            info=new_info,
        )
