"""Annotation engine with backend auto-detection.

Composes consequence assignment, gnomAD frequency lookups, and ClinVar
lookups into a batch-processing pipeline. Picks the fastest available
backend (pyranges, polars) at init time, falling back to pure-Python
when the optional deps aren't installed.
"""

from __future__ import annotations

import logging
import re as _re
from collections.abc import Iterator
from itertools import islice
from pathlib import Path
from typing import Any

from vartriage._internal.path_safety import resolve_path
from vartriage.models.config import AnnotationConfig
from vartriage.models.variant import (
    AnnotatedVariant,
    ClinVarAssertion,
    ProteinChange,
    Variant,
)
from vartriage.models.warnings import MissingDataWarning
from vartriage.protocols import ClinVarDatabase, FrequencyDatabase, IntervalIndex

logger = logging.getLogger(__name__)


def _pyranges_available() -> bool:
    """Check whether pyranges is importable."""
    try:
        import pyranges  # noqa: F401

        return True
    except ImportError:
        return False


def _polars_available() -> bool:
    """Check whether polars is importable."""
    try:
        import polars  # noqa: F401

        return True
    except ImportError:
        return False


_GTF_ATTR_CACHE: dict[str, _re.Pattern[str]] = {}


def _extract_gtf_attr(attrs: str, key: str) -> str | None:
    """Extract a single attribute value from a GTF attributes string."""
    pattern = _GTF_ATTR_CACHE.get(key)
    if pattern is None:
        pattern = _re.compile(rf'{key}\s+"([^"]+)"')
        _GTF_ATTR_CACHE[key] = pattern
    match = pattern.search(attrs)
    return match.group(1) if match else None


class AnnotationEngine:
    """Annotates variants with consequence, frequency, and ClinVar data.

    Processes in configurable batch sizes. Picks the fastest backend
    at construction time:

    - Consequence: pyranges if available, else pure-Python sorted intervals
    - Frequency/ClinVar: polars if available, else dict-based lookups

    Parameters
    ----------
    config : AnnotationConfig
        Reference file paths and batch_size.

    Raises
    ------
    FileNotFoundError
        If a required reference file is missing.
    ValueError
        If batch_size is outside [1_000, 100_000].
    """

    def __init__(self, config: AnnotationConfig) -> None:
        self._config = config
        self._warnings: list[MissingDataWarning] = []

        # Validate file paths upfront (fail-fast)
        self._validate_paths(config)

        # Initialize consequence annotator
        self._consequence_annotator: IntervalIndex = self._build_consequence_annotator(
            config.gene_annotation_path
        )

        # Attach CodonResolver for proper amino acid-level consequence calling
        if config.reference_fasta_path is not None:
            self._attach_codon_resolver(config.reference_fasta_path)

        # Initialize normalizer for consistent database lookups
        self._normalizer: object = None
        if config.reference_fasta_path is not None:
            from vartriage._internal.normalizer import VariantNormalizer

            self._normalizer = VariantNormalizer(config.reference_fasta_path)
            logger.info(
                "VariantNormalizer active: indels will be left-aligned before lookups"
            )

        # Initialize frequency database
        self._frequency_db: FrequencyDatabase | None = None
        if config.gnomad_path is not None:
            self._frequency_db = self._build_frequency_db(config.gnomad_path)

        # Initialize ClinVar database (optional)
        self._clinvar_db: ClinVarDatabase | None = self._build_clinvar_db(
            config.clinvar_path
        )

        # Capability flags (determined once at init, not per-call hasattr)
        self._supports_batch_gene_names: bool = hasattr(
            self._consequence_annotator, "gene_names_batch"
        )
        self._supports_batch_cds: bool = hasattr(
            self._consequence_annotator, "cds_overlaps_batch"
        )

        # Build CodonResolver for pyranges backend protein change resolution.
        # For the pure-Python backend, the resolver is attached via
        # _attach_codon_resolver above; for pyranges, we build it here.
        self._pyranges_codon_resolver: Any = None
        if self._supports_batch_cds and config.reference_fasta_path is not None:
            self._pyranges_codon_resolver = self._build_pyranges_codon_resolver()

    @property
    def warnings(self) -> list[MissingDataWarning]:
        """Warnings accumulated for variants missing from references."""
        return self._warnings

    def set_frequency_db(self, db: FrequencyDatabase) -> None:
        """Replace the frequency database with an external implementation.

        Used by the pipeline to inject a remote tabix gnomAD backend
        when no local gnomAD file is available.

        Parameters
        ----------
        db : FrequencyDatabase
            Replacement frequency database (must satisfy the protocol).
        """
        self._frequency_db = db

    @property
    def has_frequency_db(self) -> bool:
        """Whether a frequency database is available."""
        return self._frequency_db is not None

    def annotate(self, variants: Iterator[Variant]) -> Iterator[AnnotatedVariant]:
        """Annotate variants with consequence, frequency, and ClinVar data.

        Processes variants in batches of ``config.batch_size`` (default
        10,000). Each batch undergoes consequence assignment, frequency
        lookup, and ClinVar lookup via the Protocol interfaces selected
        at construction time.

        Parameters
        ----------
        variants : Iterator[Variant]
            Input stream of raw variants.

        Yields
        ------
        AnnotatedVariant
            Variants enriched with functional consequence, allele
            frequency, and ClinVar clinical significance.
        """
        batch_size = self._config.batch_size

        while True:
            batch = list(islice(variants, batch_size))
            if not batch:
                break

            yield from self._annotate_batch(batch)

    def _annotate_batch(self, batch: list[Variant]) -> list[AnnotatedVariant]:
        """Run consequence + frequency + ClinVar on a single batch."""
        # Consequence assignment (uses original coordinates for overlap)
        consequences = self._consequence_annotator.assign_batch(batch)

        # Gene name + protein change extraction via overlap queries
        gene_names, protein_changes = self._extract_gene_and_protein(batch)

        # Normalize coordinates for database lookups (gnomAD, ClinVar)
        # Consequence calling uses original coords (GTF overlap is position-based)
        # but frequency/ClinVar lookups need normalized coords for matching
        variant_keys = [(v.chrom, v.pos, v.ref, v.alt) for v in batch]
        if self._normalizer is not None:
            variant_keys = [
                self._normalizer.normalize(c, p, r, a)  # type: ignore[attr-defined]
                for c, p, r, a in variant_keys
            ]

        if self._frequency_db is not None:
            frequencies = self._frequency_db.lookup_batch(variant_keys)
        else:
            frequencies = [None] * len(batch)

        # ClinVar lookup
        clinvar_assertions: list[ClinVarAssertion | None] = []
        if self._clinvar_db is not None:
            clinvar_assertions = self._clinvar_db.lookup_batch(variant_keys)
        else:
            clinvar_assertions = [None] * len(batch)

        # Compose results
        results: list[AnnotatedVariant] = []
        for i, variant in enumerate(batch):
            freq = frequencies[i]
            clinvar = clinvar_assertions[i]

            frequency_unknown = freq is None
            clinvar_unknown = clinvar is None

            # Emit warnings for missing data
            if frequency_unknown:
                self._warnings.append(
                    MissingDataWarning(
                        chrom=variant.chrom,
                        pos=variant.pos,
                        ref=variant.ref,
                        alt=variant.alt,
                        source="gnomAD",
                        reason="not_found",
                    )
                )

            if clinvar_unknown and self._clinvar_db is not None:
                self._warnings.append(
                    MissingDataWarning(
                        chrom=variant.chrom,
                        pos=variant.pos,
                        ref=variant.ref,
                        alt=variant.alt,
                        source="ClinVar",
                        reason="not_found",
                    )
                )

            results.append(
                AnnotatedVariant(
                    variant=variant,
                    consequence=consequences[i],
                    allele_frequency=freq,
                    clinvar_assertion=clinvar,
                    frequency_unknown=frequency_unknown,
                    clinvar_unknown=clinvar_unknown,
                    gene_name=gene_names[i],
                    protein_change=protein_changes[i],
                )
            )

        return results

    def _extract_gene_and_protein(
        self, batch: list[Variant]
    ) -> tuple[list[str | None], list[ProteinChange | None]]:
        """Extract gene names and protein changes for a batch of variants.

        Uses the vectorized gene_names_batch() when available (pyranges
        backend) for O(1) bulk join instead of per-variant overlap queries.
        For protein changes, uses batched CDS overlap detection + CodonResolver
        on just the CDS-overlapping variants.

        Falls back to per-variant overlap() only for backends without
        vectorized methods (pure-Python SortedArrayIntervalIndex).

        Parameters
        ----------
        batch : list[Variant]
            Variants to look up.

        Returns
        -------
        tuple[list[Optional[str]], list[Optional[ProteinChange]]]
            Gene names and protein changes positionally matched to the batch.
        """
        # Fast path: use vectorized gene_names_batch when available
        # (pyranges backend). This does ONE interval join for the whole
        # batch instead of N individual joins.
        if self._supports_batch_gene_names:
            fast_gene_names = self._consequence_annotator.gene_names_batch(batch)  # type: ignore[attr-defined]

            # Protein changes: use batched CDS overlap + CodonResolver
            fast_protein_changes = self._resolve_protein_changes_batched(batch)

            # If a protein change was found, override gene_name with the
            # gene from the transcript that produced it (consistency rule)
            for i, pc in enumerate(fast_protein_changes):
                if pc is not None:
                    fast_gene_names[i] = pc.gene_name

            return fast_gene_names, fast_protein_changes

        # Slow path: per-variant overlap (pure-Python backend with codon resolver)
        gene_names: list[str | None] = []
        protein_changes: list[ProteinChange | None] = []

        for variant in batch:
            overlaps = self._consequence_annotator.overlap(
                chrom=variant.chrom,
                pos=variant.pos,
                ref=variant.ref,
                alt=variant.alt,
            )
            if not overlaps:
                gene_names.append(None)
                protein_changes.append(None)
                continue

            protein_change = self._first_nonsynonymous_protein_change(overlaps)
            if protein_change is not None:
                gene_names.append(protein_change.gene_name)
            else:
                gene_names.append(overlaps[0].get("gene_name"))
            protein_changes.append(protein_change)

        return gene_names, protein_changes

    def _resolve_protein_changes_batched(
        self, batch: list[Variant]
    ) -> list[ProteinChange | None]:
        """Resolve protein changes using batched CDS detection + CodonResolver.

        Steps:
        1. Use cds_overlaps_batch() for ONE bulk interval join to find
           which variants overlap CDS regions and their transcript IDs.
        2. For only those variants (typically <5% of a chr22 VCF), call
           CodonResolver.resolve() to get the amino acid change.

        If no reference FASTA is configured (no CodonResolver available),
        returns all None (protein changes cannot be determined without
        a reference genome).
        """
        protein_changes: list[ProteinChange | None] = [None] * len(batch)

        if not self._supports_batch_cds or self._pyranges_codon_resolver is None:
            return protein_changes

        # Step 1: batched CDS overlap detection (one bulk join)
        cds_overlaps = self._consequence_annotator.cds_overlaps_batch(batch)  # type: ignore[attr-defined]

        # Step 2: resolve protein changes only for CDS-overlapping variants
        for i, transcript_ids in enumerate(cds_overlaps):
            if not transcript_ids:
                continue

            variant = batch[i]
            # Only resolve SNVs (CodonResolver handles single-base substitutions)
            if len(variant.ref) != 1 or len(variant.alt) != 1:
                continue

            # Try each overlapping transcript until we get a non-synonymous change
            for tid in transcript_ids:
                ctx = self._pyranges_codon_resolver.resolve(
                    chrom=variant.chrom,
                    pos=variant.pos,
                    ref=variant.ref,
                    alt=variant.alt,
                    transcript_id=tid,
                )
                if ctx is not None and not ctx.is_synonymous:
                    protein_changes[i] = ProteinChange(
                        gene_name=ctx.gene_name,
                        position=ctx.codon_index + 1,
                        reference_aa=ctx.reference_aa,
                        altered_aa=ctx.altered_aa,
                    )
                    break

        return protein_changes

    def _build_pyranges_codon_resolver(self) -> Any:
        """Build a CodonResolver for the pyranges backend path.

        Called once at init time. Returns None if construction fails.
        """
        try:
            if self._config.gene_annotation_path is None:
                logger.info("No gene annotation path, skipping CodonResolver")
                return None

            from vartriage.annotation.codon_resolver import CodonResolver
            from vartriage.annotation.transcript_index import TranscriptCDSIndex

            # Build transcript index by parsing CDS features from GTF
            transcript_index = TranscriptCDSIndex()
            gtf_path = self._config.gene_annotation_path

            with open(gtf_path, encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 9 or parts[2] != "CDS":
                        continue

                    chrom = parts[0]
                    start = int(parts[3]) - 1  # GTF is 1-based
                    end = int(parts[4])  # GTF end is inclusive, convert to exclusive
                    strand = parts[6]
                    try:
                        frame = int(parts[7]) if parts[7] != "." else 0
                    except ValueError:
                        frame = 0

                    # Parse attributes for gene_name and transcript_id
                    attrs = parts[8]
                    gene_name = _extract_gtf_attr(attrs, "gene_name")
                    transcript_id = _extract_gtf_attr(attrs, "transcript_id")
                    if not transcript_id:
                        continue

                    transcript_index.add_cds_exon(
                        transcript_id=transcript_id,
                        gene_name=gene_name or "unknown",
                        chrom=chrom,
                        start=start,
                        end=end,
                        strand=strand,
                        frame=frame,
                    )

            transcript_index.finalize()

            resolver = CodonResolver(
                self._config.reference_fasta_path,  # type: ignore[arg-type]
                transcript_index,
            )
            logger.info(
                "CodonResolver created for pyranges backend (FASTA: %s)",
                self._config.reference_fasta_path,
            )
            return resolver
        except Exception as exc:
            logger.warning(
                "Could not create CodonResolver for pyranges backend: %s", exc
            )
            return None

    @staticmethod
    def _first_nonsynonymous_protein_change(
        overlaps: list[dict[str, Any]],
    ) -> ProteinChange | None:
        """Extract ProteinChange from the first non-synonymous codon context."""
        for overlap in overlaps:
            ctx = overlap.get("codon_context")
            if ctx is None or ctx.is_synonymous:
                continue
            return ProteinChange(
                gene_name=ctx.gene_name,
                position=ctx.codon_index + 1,
                reference_aa=ctx.reference_aa,
                altered_aa=ctx.altered_aa,
            )
        return None

    def _extract_gene_names(self, batch: list[Variant]) -> list[str | None]:
        """Extract gene names for a batch of variants.

        Delegates to _extract_gene_and_protein for consistency.
        """
        gene_names, _ = self._extract_gene_and_protein(batch)
        return gene_names

    def _validate_paths(self, config: AnnotationConfig) -> None:
        """Fail fast if any required reference file is missing."""
        if not config.gene_annotation_path.exists():
            raise FileNotFoundError(
                f"Gene annotation file not found: {config.gene_annotation_path}"
            )

        if config.gnomad_path is not None and not config.gnomad_path.exists():
            raise FileNotFoundError(
                f"gnomAD reference file not found: {config.gnomad_path}"
            )

        if config.clinvar_path is not None and not config.clinvar_path.exists():
            raise FileNotFoundError(
                f"ClinVar reference file not found: {config.clinvar_path}"
            )

        if config.reference_fasta_path is not None:
            if not config.reference_fasta_path.exists():
                raise FileNotFoundError(
                    f"Reference FASTA not found: {config.reference_fasta_path}"
                )
            fai_path = Path(str(config.reference_fasta_path) + ".fai")
            if not fai_path.exists():
                raise FileNotFoundError(
                    f"FASTA index (.fai) not found: {fai_path}. "
                    f"Run 'samtools faidx {config.reference_fasta_path}' to create it."
                )

    def _build_consequence_annotator(self, annotation_path: Path) -> IntervalIndex:
        """Pick the best consequence annotator available.

        Tries pyranges first, falls back to the pure-Python interval tree.

        Parameters
        ----------
        annotation_path : Path
            GTF/GFF gene annotation file.

        Returns
        -------
        IntervalIndex
            Loaded consequence annotator.
        """
        annotation_path = resolve_path(annotation_path)
        if _pyranges_available():
            try:
                from vartriage.annotation.consequence_pyranges import (
                    PyRangesConsequenceAnnotator,
                )

                logger.info("Using pyranges backend for consequence annotation")
                return PyRangesConsequenceAnnotator(annotation_path)
            except Exception as exc:
                logger.warning(
                    "pyranges backend failed, falling back to pure-Python: %s",
                    exc,
                )

        from vartriage.annotation.consequence import ConsequenceAnnotator

        logger.info("Using pure-Python backend for consequence annotation")
        return ConsequenceAnnotator(annotation_path)

    def _attach_codon_resolver(self, fasta_path: Path) -> None:
        """Create and attach a CodonResolver for proper consequence calling.

        Only works with the pure-Python SortedArrayIntervalIndex backend
        (which exposes set_codon_resolver). The pyranges backend handles
        consequence differently and doesn't need this.
        """
        annotator = self._consequence_annotator
        if not hasattr(annotator, "set_codon_resolver"):
            logger.info(
                "Codon resolver not supported by %s backend, using positional heuristic",
                type(annotator).__name__,
            )
            return

        transcript_index = getattr(annotator, "transcript_index", None)
        if transcript_index is None:
            logger.warning("No TranscriptCDSIndex available, skipping codon resolution")
            return

        from vartriage.annotation.codon_resolver import CodonResolver

        resolver = CodonResolver(fasta_path, transcript_index)
        annotator.set_codon_resolver(resolver)
        logger.info("CodonResolver attached: using FASTA-backed consequence calling")

    def _build_frequency_db(self, gnomad_path: Path) -> FrequencyDatabase:
        """Pick the best frequency database available.

        Selects backend based on file extension:
        - .vcf.bgz / .vcf.gz → tabix VCF backend (on-the-fly queries)
        - .tsv / .tsv.gz → polars if available, otherwise pure-Python dict

        Parameters
        ----------
        gnomad_path : Path
            gnomAD reference file (VCF or TSV).

        Returns
        -------
        FrequencyDatabase
            Loaded frequency database.
        """
        gnomad_name = gnomad_path.name
        if gnomad_name.endswith((".vcf.bgz", ".vcf.gz")):
            from vartriage.annotation.frequency_tabix import TabixFrequencyDatabase

            logger.info("Using tabix VCF backend for frequency lookup")
            freq_db: FrequencyDatabase = TabixFrequencyDatabase()
            freq_db.load(gnomad_path)
            return freq_db

        if _polars_available():
            try:
                from vartriage.annotation.frequency_polars import (
                    PolarsFrequencyDatabase,
                )

                logger.info("Using polars backend for frequency lookup")
                polars_db: FrequencyDatabase = PolarsFrequencyDatabase()
                polars_db.load(gnomad_path)
                return polars_db
            except Exception as exc:
                logger.warning(
                    "polars frequency backend failed, falling back to pure-Python: %s",
                    exc,
                )

        from vartriage.annotation.frequency import DictFrequencyDatabase

        logger.info("Using pure-Python backend for frequency lookup")
        freq_db = DictFrequencyDatabase()
        freq_db.load(gnomad_path)
        return freq_db

    def _build_clinvar_db(self, clinvar_path: Path | None) -> ClinVarDatabase | None:
        """Pick the best ClinVar backend, or None if no path given.

        Parameters
        ----------
        clinvar_path : Optional[Path]
            ClinVar reference file, or None to skip.

        Returns
        -------
        Optional[ClinVarDatabase]
            Loaded ClinVar database, or None.
        """
        if clinvar_path is None:
            return None

        if _polars_available():
            try:
                from vartriage.annotation.clinvar_polars import PolarsClinVarDatabase

                logger.info("Using polars backend for ClinVar lookup")
                clinvar_db: ClinVarDatabase = PolarsClinVarDatabase()
                clinvar_db.load(clinvar_path)
                return clinvar_db
            except Exception as exc:
                logger.warning(
                    "polars ClinVar backend failed, falling back to pure-Python: %s",
                    exc,
                )

        from vartriage.annotation.clinvar import DictClinVarDatabase

        logger.info("Using pure-Python backend for ClinVar lookup")
        clinvar_db = DictClinVarDatabase()
        clinvar_db.load(clinvar_path)
        return clinvar_db
