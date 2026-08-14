"""Named URL presets for common remote tabix-indexed databases.

Users pass a preset name (e.g., "cadd-v1.7-grch38") instead of memorizing
full download URLs. The resolver checks the registry and returns the
canonical URL, or passes through a raw URL unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True)
class PresetEntry:
    """A single preset definition.

    Parameters
    ----------
    name : str
        Short identifier used on the CLI (e.g., "cadd-v1.7-grch38").
    url : str
        Full download URL. May contain a {chrom} placeholder for
        per-chromosome resources.
    source : str
        Which scoring system this serves ("cadd" or "gnomad").
    genome_build : str
        Reference build ("grch38" or "grch37").
    description : str
        Human-readable one-liner shown by list-presets.
    """

    name: str
    url: str
    source: Literal["cadd", "gnomad"]
    genome_build: str
    description: str


_PRESET_REGISTRY: tuple[PresetEntry, ...] = (
    PresetEntry(
        name="cadd-v1.7-grch38",
        url="https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz",
        source="cadd",
        genome_build="grch38",
        description="CADD v1.7 all possible SNVs, GRCh38",
    ),
    PresetEntry(
        name="cadd-v1.7-grch37",
        url="https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh37/whole_genome_SNVs.tsv.gz",
        source="cadd",
        genome_build="grch37",
        description="CADD v1.7 all possible SNVs, GRCh37",
    ),
    PresetEntry(
        name="cadd-v1.7-indels-grch38",
        url="https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz",
        source="cadd",
        genome_build="grch38",
        description="CADD v1.7 pre-scored gnomAD indels, GRCh38",
    ),
    PresetEntry(
        name="gnomad-exomes-v4-grch38",
        url="https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1.1/vcf/exomes/gnomad.exomes.v4.1.1.sites.{chrom}.vcf.bgz",
        source="gnomad",
        genome_build="grch38",
        description="gnomAD exomes v4.1.1 per-chromosome VCF, GRCh38",
    ),
)

_PRESET_BY_NAME: dict[str, PresetEntry] = {p.name: p for p in _PRESET_REGISTRY}


def resolve_preset(name_or_url: str) -> str:
    """Resolve a preset name to its URL, or pass through a raw URL.

    Parameters
    ----------
    name_or_url : str
        Either a registered preset name or a full URL starting with
        http:// or https://.

    Returns
    -------
    str
        The resolved URL.

    Raises
    ------
    ValueError
        If the string is not a valid URL and not a registered preset.
    """
    if name_or_url.startswith(("http://", "https://")):
        return name_or_url

    preset = _PRESET_BY_NAME.get(name_or_url)
    if preset is None:
        available = ", ".join(sorted(_PRESET_BY_NAME.keys()))
        raise ValueError(
            f"Unknown preset '{name_or_url}'. Available presets: {available}"
        )
    return preset.url


def get_preset(name: str) -> PresetEntry | None:
    """Look up a preset by name. Returns None if not found."""
    return _PRESET_BY_NAME.get(name)


def list_presets(source: str | None = None) -> list[PresetEntry]:
    """Return all registered presets, optionally filtered by source.

    Parameters
    ----------
    source : str | None
        Filter to "cadd" or "gnomad". None returns all.

    Returns
    -------
    list[PresetEntry]
        Matching presets sorted by name.
    """
    entries = list(_PRESET_REGISTRY)
    if source is not None:
        entries = [p for p in entries if p.source == source]
    return sorted(entries, key=lambda p: p.name)


def is_preset_name(value: str) -> bool:
    """Check whether a string is a registered preset name."""
    return value in _PRESET_BY_NAME
