"""Data models for clinical report sections.

All dataclasses are frozen (immutable after creation) and represent
the structured content for each report section.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class HeaderData:
    """Report header section data.

    Parameters
    ----------
    patient_id : str
        Patient identifier.
    panel_name : str
        Gene panel name used for the analysis.
    analysis_date : str
        ISO 8601 formatted date of the analysis.
    pipeline_version : str
        Version string of the vartriage pipeline.
    """

    patient_id: str
    panel_name: str
    analysis_date: str
    pipeline_version: str


@dataclass(frozen=True)
class ExecutiveSummaryData:
    """Executive summary section data with variant counts.

    Parameters
    ----------
    total_variants_analyzed : int
        Total number of variants in the input.
    variants_passed_filters : int
        Number of variants that passed quality filters.
    pathogenic_count : int
        Count of variants classified as Pathogenic.
    likely_pathogenic_count : int
        Count of variants classified as Likely Pathogenic.
    vus_count : int
        Count of variants classified as VUS.
    """

    total_variants_analyzed: int
    variants_passed_filters: int
    pathogenic_count: int
    likely_pathogenic_count: int
    vus_count: int


@dataclass(frozen=True)
class FindingsRow:
    """Single row in the findings table.

    Parameters
    ----------
    gene_name : str | None
        Gene symbol, or None for intergenic variants.
    consequence : str
        Functional consequence (e.g., Missense, Frameshift).
    classification : str
        ACMG/AMP classification tier.
    composite_rank : float | None
        Composite pathogenicity rank, or None if unavailable.
    chromosome : str
        Chromosome name.
    position : int
        1-based genomic position.
    """

    gene_name: str | None
    consequence: str
    classification: str
    composite_rank: float | None
    chromosome: str
    position: int


@dataclass(frozen=True)
class EvidenceCardData:
    """Per-variant evidence card data.

    Parameters
    ----------
    gene_name : str | None
        Gene symbol, or None for intergenic variants.
    consequence : str
        Functional consequence string.
    allele_frequency_formatted : str | None
        Formatted AF string with denominator context, or None.
    predictor_scores_formatted : list[str]
        List of formatted predictor score strings.
    clinvar_assertion : str | None
        ClinVar clinical significance, or None.
    inheritance_pattern : str | None
        Inheritance mode if available, or None.
    evidence_tags_with_explanations : list[str]
        ACMG tags each followed by plain-language explanation.
    narrative : str
        Full evidence narrative text for this variant.
    """

    gene_name: str | None
    consequence: str
    allele_frequency_formatted: str | None
    predictor_scores_formatted: list[str]
    clinvar_assertion: str | None
    inheritance_pattern: str | None
    evidence_tags_with_explanations: list[str]
    narrative: str


@dataclass(frozen=True)
class MethodologyData:
    """Methodology section data.

    Parameters
    ----------
    pipeline_version : str
        Version string of the pipeline.
    reference_files : dict[str, str]
        Mapping of reference file paths to SHA-256 checksums.
    classification_parameters : dict[str, str]
        Key classification parameters and their values.
    analysis_timestamp : str
        ISO 8601 timestamp of the analysis run.
    """

    pipeline_version: str
    reference_files: dict[str, str]
    classification_parameters: dict[str, str]
    analysis_timestamp: str


@dataclass(frozen=True)
class SignOffData:
    """Sign-off section with placeholder fields.

    Parameters
    ----------
    reviewer_name_placeholder : str
        Placeholder for the reviewer name.
    review_date_placeholder : str
        Placeholder for the review date.
    digital_signature_placeholder : str
        Placeholder for the digital signature.
    """

    reviewer_name_placeholder: str = "[Reviewer Name]"
    review_date_placeholder: str = "[Review Date]"
    digital_signature_placeholder: str = "[Digital Signature]"


@dataclass(frozen=True)
class ReportSections:
    """Container for all assembled report section data.

    Parameters
    ----------
    header : HeaderData
        Report header content.
    executive_summary : ExecutiveSummaryData
        Summary statistics.
    findings_table : list[FindingsRow]
        Ranked variant findings.
    evidence_cards : list[EvidenceCardData]
        Per-variant evidence cards.
    limitations : list[str]
        Data source limitations encountered.
    methodology : MethodologyData
        Analysis methodology details.
    sign_off : SignOffData
        Sign-off placeholder section.
    """

    header: HeaderData
    executive_summary: ExecutiveSummaryData
    findings_table: list[FindingsRow]
    evidence_cards: list[EvidenceCardData]
    limitations: list[str]
    methodology: MethodologyData
    sign_off: SignOffData = field(default_factory=SignOffData)
