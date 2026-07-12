"""Evidence narrative builder for clinical variant reports.

Transforms raw ClassifiedVariant data into human-readable clinical
evidence narratives using hardcoded string templates with data
interpolation. No LLM or generative AI is invoked.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

from vartriage.models.variant import (
    ClassifiedVariant,
    EvidenceTag,
    FunctionalConsequence,
)

if TYPE_CHECKING:
    pass


# Template constants for narrative construction.
_GENE_CONSEQUENCE_TEMPLATE = "{gene}: {consequence} at {chrom}:{pos} ({ref}>{alt})."
_INTERGENIC_CONSEQUENCE_TEMPLATE = (
    "Intergenic variant at {chrom}:{pos} ({ref}>{alt})."
)
_AF_TEMPLATE = "Population frequency: {af_formatted} in gnomAD."
_AF_MISSING_NOTE = "Population frequency data not available from gnomAD."
_CLINVAR_TEMPLATE = "ClinVar: {assertion}."
_CLINVAR_MISSING_NOTE = "ClinVar data not available for this variant."
_INHERITANCE_TEMPLATE = "Inheritance pattern: {pattern}."
_CLASSIFICATION_TEMPLATE = "Classification: {classification}."
_ACMG_CRITERIA_HEADER = "ACMG criteria met:"
_DATA_GAP_TEMPLATE = "{source} score not available."


# Predictor score scale context for known score types.
_SCORE_CONTEXT: dict[str, str] = {
    "REVEL": "scale 0-1, threshold 0.7",
    "CADD": "Phred-scaled",
    "SpliceAI": "scale 0-1, threshold 0.2",
}

# Plain-language explanations for each ACMG evidence tag.
_TAG_EXPLANATIONS: dict[EvidenceTag, str] = {
    EvidenceTag.PVS1: "null variant, {consequence_detail} truncates protein",
    EvidenceTag.PM2: "absent/rare in population databases",
    EvidenceTag.PP3: "computational predictors support a deleterious effect",
    EvidenceTag.PP5: "ClinVar pathogenic assertion",
}

# Consequence descriptions for PVS1 tag explanations.
_CONSEQUENCE_DESCRIPTIONS: dict[FunctionalConsequence, str] = {
    FunctionalConsequence.FRAMESHIFT: "frameshift",
    FunctionalConsequence.NONSENSE: "nonsense",
    FunctionalConsequence.SPLICE_SITE: "splice site disruption",
    FunctionalConsequence.MISSENSE: "missense",
    FunctionalConsequence.IN_FRAME_INSERTION: "in-frame insertion",
    FunctionalConsequence.IN_FRAME_DELETION: "in-frame deletion",
    FunctionalConsequence.SYNONYMOUS: "synonymous",
    FunctionalConsequence.INTERGENIC: "intergenic",
}


class EvidenceNarrativeBuilder:
    """Produces human-readable evidence narratives from variant data.

    Uses hardcoded string templates with data interpolation.
    No LLM or generative AI is invoked.
    """

    BANNED_VOCABULARY: frozenset[str] = frozenset({
        "leverage",
        "comprehensive",
        "robust",
        "seamlessly",
        "furthermore",
        "utilize",
        "facilitate",
        "delve",
        "foster",
        "crucial",
        "it is important to note",
        "notably",
        "moreover",
        "additionally",
    })

    EM_DASH = "\u2014"

    def build_narrative(self, variant: ClassifiedVariant) -> str:
        """Generate a complete evidence narrative for a single variant.

        Parameters
        ----------
        variant : ClassifiedVariant
            The classified variant to narrate.

        Returns
        -------
        str
            A multi-sentence narrative describing the variant evidence.
        """
        parts: list[str] = []

        annotated = variant.scored.annotated
        raw_variant = annotated.variant

        # Gene and consequence line.
        if annotated.gene_name is not None:
            parts.append(
                _GENE_CONSEQUENCE_TEMPLATE.format(
                    gene=annotated.gene_name,
                    consequence=annotated.consequence.value,
                    chrom=raw_variant.chrom,
                    pos=raw_variant.pos,
                    ref=raw_variant.ref,
                    alt=raw_variant.alt,
                )
            )
        else:
            parts.append(
                _INTERGENIC_CONSEQUENCE_TEMPLATE.format(
                    chrom=raw_variant.chrom,
                    pos=raw_variant.pos,
                    ref=raw_variant.ref,
                    alt=raw_variant.alt,
                )
            )

        # Allele frequency.
        if annotated.allele_frequency is not None:
            af_str = self.format_allele_frequency(annotated.allele_frequency)
            parts.append(_AF_TEMPLATE.format(af_formatted=af_str))
        else:
            parts.append(_AF_MISSING_NOTE)

        # Computational predictor scores.
        self._append_score_parts(parts, variant.scored)

        # ClinVar assertion.
        if annotated.clinvar_assertion is not None:
            parts.append(
                _CLINVAR_TEMPLATE.format(
                    assertion=annotated.clinvar_assertion.value
                )
            )
        else:
            parts.append(_CLINVAR_MISSING_NOTE)

        # Inheritance pattern from variant info dict.
        inheritance = raw_variant.info.get("inheritance_pattern")
        if inheritance is not None:
            parts.append(
                _INHERITANCE_TEMPLATE.format(pattern=inheritance)
            )

        # Evidence tags with explanations.
        if variant.evidence_tags:
            tag_explanations: list[str] = []
            for tag in sorted(variant.evidence_tags, key=lambda t: t.value):
                tag_explanations.append(
                    self.format_evidence_tag(tag, variant)
                )
            parts.append(
                _ACMG_CRITERIA_HEADER + " " + ", ".join(tag_explanations) + "."
            )

        # Final classification.
        parts.append(
            _CLASSIFICATION_TEMPLATE.format(
                classification=variant.classification.value
            )
        )

        narrative = " ".join(parts)

        self._validate_output(narrative)
        return narrative

    def _append_score_parts(
        self, parts: list[str], scored: "ScoredVariant",
    ) -> None:
        """Append predictor score sentences to the narrative parts."""
        score_entries: list[str] = []
        has_score = False

        for name, value in [
            ("REVEL", scored.revel_score),
            ("CADD", scored.cadd_phred),
            ("SpliceAI", scored.spliceai_score),
        ]:
            if value is not None:
                score_entries.append(
                    self.format_predictor_score(name, value)
                )
                has_score = True
            elif name != "SpliceAI":
                score_entries.append(
                    _DATA_GAP_TEMPLATE.format(source=name)
                )

        if has_score:
            parts.append(
                "Computational evidence: " + ", ".join(
                    e for e in score_entries if "not available" not in e
                ) + "."
            )
        for entry in score_entries:
            if "not available" in entry:
                parts.append(entry)

    def format_allele_frequency(self, af: float) -> str:
        """Format AF with decimal representation and denominator context.

        Parameters
        ----------
        af : float
            Allele frequency value in range (0, 1].

        Returns
        -------
        str
            Formatted string like "0.000008 (1 in 125,000)".
        """
        if af <= 0:
            return "0 (absent)"

        if af >= 1.0:
            return "1.0 (1 in 1)"

        denominator = round(1.0 / af)
        denominator_str = f"{denominator:,}"

        # Format the AF value: use enough decimal places to show
        # at least 2 significant figures.
        sig_digits = max(2, -int(math.floor(math.log10(af))) + 1)
        af_str = f"{af:.{sig_digits}f}"

        return f"{af_str} (1 in {denominator_str})"

    def format_predictor_score(
        self, score_name: str, value: float
    ) -> str:
        """Format a predictor score with its standard scale context.

        Parameters
        ----------
        score_name : str
            Name of the score (e.g., "REVEL", "CADD").
        value : float
            Numeric score value.

        Returns
        -------
        str
            Formatted string like "REVEL 0.95 (scale 0-1, threshold 0.7)".
        """
        context = _SCORE_CONTEXT.get(score_name, "")
        if score_name == "CADD":
            value_str = f"{value:.1f}"
        else:
            value_str = f"{value:.2f}"

        if context:
            return f"{score_name} {value_str} ({context})"
        return f"{score_name} {value_str}"

    def format_evidence_tag(
        self, tag: EvidenceTag, variant: ClassifiedVariant
    ) -> str:
        """Format one evidence tag with a plain-language explanation.

        Parameters
        ----------
        tag : EvidenceTag
            The ACMG evidence tag to explain.
        variant : ClassifiedVariant
            The variant providing context for the explanation.

        Returns
        -------
        str
            Formatted tag like "PVS1 (null variant, frameshift truncates protein)".
        """
        explanation_template = _TAG_EXPLANATIONS.get(tag)
        if explanation_template is None:
            return f"{tag.value} (evidence criterion met)"

        consequence = variant.scored.annotated.consequence
        consequence_detail = _CONSEQUENCE_DESCRIPTIONS.get(
            consequence, consequence.value.lower()
        )

        explanation = explanation_template.format(
            consequence_detail=consequence_detail
        )
        return f"{tag.value} ({explanation})"

    def _validate_output(self, text: str) -> None:
        """Assert the narrative contains no banned words or em dashes.

        Parameters
        ----------
        text : str
            The generated narrative text to validate.

        Raises
        ------
        AssertionError
            If banned content is found in the narrative.
        """
        assert self.EM_DASH not in text, (
            f"Narrative contains em dash (U+2014): {text!r}"
        )

        text_lower = text.lower()
        for word in self.BANNED_VOCABULARY:
            assert word not in text_lower, (
                f"Narrative contains banned vocabulary '{word}': {text!r}"
            )
