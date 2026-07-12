"""Sample-based genotype extraction and filtering for multi-sample VCFs."""

from __future__ import annotations

import warnings
from typing import Iterator

from vartriage.models.config import SampleConfig
from vartriage.models.variant import Variant
from vartriage.models.warnings import MissingDataWarning


class SampleExtractor:
    """Extract and filter variants by sample genotype.

    Keeps variants where the selected sample carries at least one
    alt allele. Attaches GT, GQ, and sample_name to Variant.info.

    Parameters
    ----------
    config : SampleConfig
        Sample name and optional min_gq threshold.
    sample_names : list[str]
        Sample names from the VCF header.

    Raises
    ------
    ValueError
        If config.sample_name is not in sample_names.
    """

    def __init__(
        self, config: SampleConfig, sample_names: list[str]
    ) -> None:
        if config.sample_name not in sample_names:
            raise ValueError(
                f"Sample '{config.sample_name}' not found in VCF. "
                f"Available samples: {sample_names}"
            )
        self._sample_name = config.sample_name
        self._min_gq = config.min_gq

    def apply(
        self, variants: Iterator[Variant]
    ) -> Iterator[Variant]:
        """Filter variants by alt allele presence and GQ threshold.

        Strips _pysam_samples from info and attaches sample_gt,
        sample_gq, sample_name instead.

        Parameters
        ----------
        variants : Iterator[Variant]
            Must have "_pysam_samples" in info (from VCFParser).

        Yields
        ------
        Variant
            Variants passing genotype and GQ filters.
        """
        for variant in variants:
            sample_data = variant.info.get("_pysam_samples", {})
            sample_entry = sample_data.get(self._sample_name, {})
            gt = sample_entry.get("GT")

            if gt is None:
                continue

            gt_str = self._format_gt(gt)

            if not self._has_alt_allele(gt_str):
                continue

            gq = sample_entry.get("GQ")

            if self._min_gq is not None:
                if gq is None:
                    reason = (
                        f"Missing GQ for sample "
                        f"'{self._sample_name}' at "
                        f"{variant.chrom}:{variant.pos}"
                    )
                    warning_data = MissingDataWarning(
                        chrom=variant.chrom,
                        pos=variant.pos,
                        ref=variant.ref,
                        alt=variant.alt,
                        source="GQ",
                        reason=reason,
                    )
                    warnings.warn(
                        UserWarning(warning_data),
                        stacklevel=2,
                    )
                    continue
                if gq < self._min_gq:
                    continue

            new_info = dict(variant.info)
            new_info.pop("_pysam_samples", None)
            new_info["sample_gt"] = gt_str
            new_info["sample_name"] = self._sample_name
            if gq is not None:
                new_info["sample_gq"] = gq

            yield Variant(
                chrom=variant.chrom,
                pos=variant.pos,
                id=variant.id,
                ref=variant.ref,
                alt=variant.alt,
                qual=variant.qual,
                filter_status=variant.filter_status,
                info=new_info,
            )

    @staticmethod
    def _format_gt(gt: tuple[int | None, ...]) -> str:
        """Convert pysam GT tuple to VCF-style string.

        Parameters
        ----------
        gt : tuple
            e.g. (0, 1), (None, None).

        Returns
        -------
        str
            e.g. "0/1", "./.".
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
            e.g. "0/1", "1/1", "./.", "0|1".

        Returns
        -------
        bool
        """
        import re

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
