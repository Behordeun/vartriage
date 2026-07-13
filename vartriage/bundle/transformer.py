"""Post-download transformers for reference files.

Converts raw downloaded files (VCF, CSV, etc.) into the TSV format
expected by vartriage's annotation and scoring engines.
"""

from __future__ import annotations

import csv
import gzip
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Protocol, Sequence


@dataclass
class TransformResult:
    """Result of a file transformation.

    Attributes
    ----------
    output_path : Path
        Path to the transformed file.
    rows_written : int
        Number of data rows in the output.
    source_path : Path
        Path to the raw input file.
    """

    output_path: Path
    rows_written: int
    source_path: Path


class TransformStrategy(Protocol):
    """Protocol for file transformation strategies."""

    def transform(self, source: Path, dest: Path, build: str) -> TransformResult:
        """Transform source file into vartriage-compatible TSV.

        Parameters
        ----------
        source : Path
            Raw downloaded file.
        dest : Path
            Output path for the transformed TSV.
        build : str
            Genome build (for chr prefix decisions).

        Returns
        -------
        TransformResult
            Metadata about the transformation.
        """
        ...


class VcfToTsvTransformer:
    """Extract columns from a VCF into vartriage TSV format.

    Attempts bcftools first (fast, handles bgzipped VCFs natively).
    Falls back to a pure-Python parser using pysam if bcftools is
    unavailable.
    """

    def __init__(
        self,
        columns: str = "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\n",
        header: str = "chrom\tpos\tref\talt\taf",
        add_chr_prefix: bool = False,
    ) -> None:
        """Configure the VCF-to-TSV transformer.

        Parameters
        ----------
        columns : str
            bcftools query format string.
        header : str
            TSV header line (tab-separated column names).
        add_chr_prefix : bool
            If True, prepend 'chr' to chromosome names that lack it.
        """
        self._columns = columns
        self._header = header
        self._add_chr_prefix = add_chr_prefix

    def transform(self, source: Path, dest: Path, _build: str) -> TransformResult:
        """Transform VCF to TSV using bcftools or pysam fallback."""
        if self._bcftools_available():
            return self._transform_bcftools(source, dest)
        return self._transform_pysam(source, dest)

    def _bcftools_available(self) -> bool:
        return shutil.which("bcftools") is not None

    def _transform_bcftools(self, source: Path, dest: Path) -> TransformResult:
        """Use bcftools query for fast extraction."""
        cmd = [
            "bcftools",
            "query",
            "-f",
            self._columns,
            str(source),
        ]

        dest.parent.mkdir(parents=True, exist_ok=True)

        with open(dest, "w", encoding="utf-8") as out:
            out.write(self._header + "\n")
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
                text=True,
            )
            rows = 0
            for line in result.stdout.splitlines():
                if not line.strip():
                    continue
                if self._add_chr_prefix and not line.startswith("chr"):
                    line = "chr" + line
                out.write(line + "\n")
                rows += 1

        return TransformResult(output_path=dest, rows_written=rows, source_path=source)

    def _transform_pysam(self, source: Path, dest: Path) -> TransformResult:
        """Pure-Python fallback using pysam."""
        try:
            import pysam
        except ImportError as exc:
            raise ImportError(
                "Neither bcftools nor pysam is available for VCF transformation. "
                "Install bcftools or run: pip install pysam"
            ) from exc

        dest.parent.mkdir(parents=True, exist_ok=True)
        rows = 0

        vcf = pysam.VariantFile(str(source))
        with open(dest, "w", encoding="utf-8") as out:
            out.write(self._header + "\n")
            for record in vcf:
                chrom = record.chrom
                if self._add_chr_prefix and not chrom.startswith("chr"):
                    chrom = "chr" + chrom
                pos = record.pos
                ref = record.ref
                for alt in record.alts or []:
                    af = record.info.get("AF", [None])[0] if "AF" in record.info else ""
                    af_str = str(af) if af is not None else "."
                    out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{af_str}\n")
                    rows += 1
        vcf.close()

        return TransformResult(output_path=dest, rows_written=rows, source_path=source)


class ClinvarVcfTransformer(VcfToTsvTransformer):
    """ClinVar-specific VCF transformer.

    Extracts CLNSIG field and normalizes clinical significance values.
    Always adds chr prefix since ClinVar GRCh38 VCF uses bare names.
    """

    def __init__(self) -> None:
        super().__init__(
            columns="%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/CLNSIG\\n",
            header="chrom\tpos\tref\talt\tclinical_significance",
            add_chr_prefix=True,
        )

    def transform(self, source: Path, dest: Path, build: str) -> TransformResult:
        """Transform ClinVar VCF with significance normalization."""
        if self._bcftools_available():
            return self._transform_clinvar_bcftools(source, dest)
        return self._transform_clinvar_pysam(source, dest)

    def _transform_clinvar_bcftools(self, source: Path, dest: Path) -> TransformResult:
        """Extract ClinVar with bcftools, normalizing CLNSIG values."""
        cmd = [
            "bcftools",
            "query",
            "-f",
            "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/CLNSIG\\n",
            str(source),
        ]

        sig_map = {
            "Pathogenic": "Pathogenic",
            "Likely_pathogenic": "Likely pathogenic",
            "Uncertain_significance": "Uncertain significance",
            "Likely_benign": "Likely benign",
            "Benign": "Benign",
        }

        dest.parent.mkdir(parents=True, exist_ok=True)
        rows = 0

        with open(dest, "w", encoding="utf-8") as out:
            out.write("chrom\tpos\tref\talt\tclinical_significance\n")
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
                text=True,
            )
            for line in result.stdout.splitlines():
                if not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    continue

                chrom = parts[0]
                if not chrom.startswith("chr"):
                    chrom = "chr" + chrom

                # Normalize significance
                raw_sig = parts[4].split("/")[0].split(",")[0]
                normalized = None
                for key, val in sig_map.items():
                    if raw_sig.startswith(key) or raw_sig.lower().startswith(
                        key.lower()
                    ):
                        normalized = val
                        break

                if normalized is None:
                    continue

                out.write(
                    f"{chrom}\t{parts[1]}\t{parts[2]}\t{parts[3]}\t{normalized}\n"
                )
                rows += 1

        return TransformResult(output_path=dest, rows_written=rows, source_path=source)

    _SIG_MAP = {
        "Pathogenic": "Pathogenic",
        "Likely_pathogenic": "Likely pathogenic",
        "Uncertain_significance": "Uncertain significance",
        "Likely_benign": "Likely benign",
        "Benign": "Benign",
    }

    def _transform_clinvar_pysam(self, source: Path, dest: Path) -> TransformResult:
        """ClinVar pysam fallback."""
        try:
            import pysam
        except ImportError as exc:
            raise ImportError(
                "Neither bcftools nor pysam available for ClinVar transformation."
            ) from exc

        dest.parent.mkdir(parents=True, exist_ok=True)
        rows = 0

        vcf = pysam.VariantFile(str(source))
        with open(dest, "w", encoding="utf-8") as out:
            out.write("chrom\tpos\tref\talt\tclinical_significance\n")
            for record in vcf:
                chrom = self._normalize_chrom(record.chrom)
                normalized = self._extract_clnsig(record)
                if normalized is None:
                    continue

                for alt in record.alts or []:
                    out.write(
                        f"{chrom}\t{record.pos}\t{record.ref}\t{alt}\t{normalized}\n"
                    )
                    rows += 1
        vcf.close()

        return TransformResult(output_path=dest, rows_written=rows, source_path=source)

    @staticmethod
    def _normalize_chrom(chrom: str) -> str:
        """Ensure chromosome has 'chr' prefix."""
        return chrom if chrom.startswith("chr") else "chr" + chrom

    def _extract_clnsig(self, record: object) -> Optional[str]:
        """Extract and normalize CLNSIG from a pysam record."""
        clnsig = record.info.get("CLNSIG", [None])  # type: ignore[attr-defined]
        if isinstance(clnsig, tuple):
            clnsig = clnsig[0] if clnsig else None
        if clnsig is None:
            return None

        raw = str(clnsig).split("/")[0].split(",")[0]
        for key, val in self._SIG_MAP.items():
            if raw.startswith(key) or raw.lower().startswith(key.lower()):
                return val
        return None


class CsvToTsvTransformer:
    """Convert CSV files (REVEL) to vartriage TSV format.

    Handles column renaming and chr prefix normalization.
    """

    def __init__(
        self,
        column_map: Optional[dict[str, str]] = None,
        add_chr_prefix: bool = True,
    ) -> None:
        """Configure CSV-to-TSV transformer.

        Parameters
        ----------
        column_map : dict, optional
            Mapping from source column names to output names.
            If None, uses REVEL defaults.
        add_chr_prefix : bool
            Prepend 'chr' to bare chromosome names.
        """
        self._column_map = column_map or {
            "chr": "chrom",
            "grch38_pos": "pos",
            "ref": "ref",
            "alt": "alt",
            "REVEL": "score",
        }
        self._add_chr_prefix = add_chr_prefix

    def transform(self, source: Path, dest: Path, _build: str) -> TransformResult:
        """Transform CSV to TSV with column renaming."""
        dest.parent.mkdir(parents=True, exist_ok=True)
        rows = 0

        opener = gzip.open if source.suffix == ".gz" else open
        open_kwargs = {"mode": "rt", "encoding": "utf-8"}

        with opener(source, **open_kwargs) as infile:  # type: ignore[call-overload]
            reader = csv.DictReader(infile)
            if reader.fieldnames is None:
                raise ValueError(f"CSV file has no header: {source}")

            source_cols, output_cols = self._map_columns(
                list(reader.fieldnames), source
            )

            with open(dest, "w", encoding="utf-8") as out:
                out.write("\t".join(output_cols) + "\n")
                for row in reader:
                    values = [self._transform_value(row, src) for src in source_cols]
                    out.write("\t".join(values) + "\n")
                    rows += 1

        return TransformResult(output_path=dest, rows_written=rows, source_path=source)

    def _map_columns(
        self, fieldnames: Sequence[str], source: Path
    ) -> tuple[list[str], list[str]]:
        """Map source columns to output columns, raising if none match."""
        source_cols = []
        output_cols = []
        for src_name, out_name in self._column_map.items():
            if src_name in fieldnames:
                source_cols.append(src_name)
                output_cols.append(out_name)

        if not output_cols:
            raise ValueError(
                f"No matching columns found in {source}. "
                f"Expected: {list(self._column_map.keys())}, "
                f"Found: {fieldnames}"
            )
        return source_cols, output_cols

    def _transform_value(self, row: dict[str, str], src: str) -> str:
        """Get value for a source column, applying chr prefix if needed."""
        val = row.get(src, "")
        if (
            self._add_chr_prefix
            and src in ("chr", "chrom")
            and val
            and not val.startswith("chr")
        ):
            return "chr" + val
        return val


class SpliceAIExtractor:
    """Extract max delta score from SpliceAI VCF INFO field.

    SpliceAI VCF contains: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|...
    We extract max(DS_AG, DS_AL, DS_DG, DS_DL) as the score.
    """

    def transform(self, source: Path, dest: Path, _build: str) -> TransformResult:
        """Extract SpliceAI max delta scores to TSV."""
        if shutil.which("bcftools"):
            return self._transform_bcftools(source, dest)
        return self._transform_pysam(source, dest)

    def _transform_bcftools(self, source: Path, dest: Path) -> TransformResult:
        """Use bcftools to extract SpliceAI INFO field."""
        cmd = [
            "bcftools",
            "query",
            "-f",
            "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/SpliceAI\\n",
            str(source),
        ]

        dest.parent.mkdir(parents=True, exist_ok=True)
        rows = 0

        with open(dest, "w", encoding="utf-8") as out:
            out.write("chrom\tpos\tref\talt\tscore\n")
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
                text=True,
            )
            for line in result.stdout.splitlines():
                if not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    continue

                chrom, pos, ref, alt, info = (
                    parts[0],
                    parts[1],
                    parts[2],
                    parts[3],
                    parts[4],
                )
                score = self._parse_max_delta(info)
                if score is not None:
                    out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{score:.4f}\n")
                    rows += 1

        return TransformResult(output_path=dest, rows_written=rows, source_path=source)

    def _transform_pysam(self, source: Path, dest: Path) -> TransformResult:
        """Pure-Python fallback using pysam."""
        try:
            import pysam
        except ImportError as exc:
            raise ImportError(
                "Neither bcftools nor pysam available for SpliceAI extraction."
            ) from exc

        dest.parent.mkdir(parents=True, exist_ok=True)
        rows = 0

        vcf = pysam.VariantFile(str(source))
        with open(dest, "w", encoding="utf-8") as out:
            out.write("chrom\tpos\tref\talt\tscore\n")
            for record in vcf:
                info_val = record.info.get("SpliceAI", None)
                if info_val is None:
                    continue
                raw = str(info_val[0]) if isinstance(info_val, tuple) else str(info_val)
                score = self._parse_max_delta(raw)
                if score is not None:
                    for alt in record.alts or []:
                        out.write(
                            f"{record.chrom}\t{record.pos}\t{record.ref}\t{alt}\t{score:.4f}\n"
                        )
                        rows += 1
        vcf.close()

        return TransformResult(output_path=dest, rows_written=rows, source_path=source)

    @staticmethod
    def _parse_max_delta(info_str: str) -> Optional[float]:
        """Parse SpliceAI INFO field and return max delta score.

        Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
        """
        try:
            # May contain multiple annotations separated by comma
            max_score = 0.0
            for annotation in info_str.split(","):
                fields = annotation.split("|")
                if len(fields) < 6:
                    continue
                scores = [float(fields[i]) for i in range(2, 6)]
                max_score = max(max_score, max(scores))
            return max_score if max_score > 0 else None
        except (ValueError, IndexError):
            return None


class PassthroughTransformer:
    """No transformation needed — copy or decompress only.

    Used for GENCODE GTF files that are used as-is by the pipeline.
    Handles .gz decompression if needed.
    """

    def transform(self, source: Path, dest: Path, _build: str) -> TransformResult:
        """Copy or decompress source to dest."""
        dest.parent.mkdir(parents=True, exist_ok=True)

        if source.suffix == ".gz":
            with gzip.open(source, "rb") as f_in, open(dest, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        else:
            shutil.copy2(source, dest)

        # Count lines (approximate row count)
        rows = sum(1 for _ in open(dest, "rb")) - 1

        return TransformResult(
            output_path=dest, rows_written=max(0, rows), source_path=source
        )


# Registry mapping transform_type strings to transformer classes
TRANSFORM_REGISTRY: dict[str, type] = {
    "vcf_to_tsv": VcfToTsvTransformer,
    "csv_to_tsv": CsvToTsvTransformer,
    "spliceai_vcf": SpliceAIExtractor,
    "none": PassthroughTransformer,
}


def get_transformer(transform_type: str, bundle_name: str = "") -> "TransformStrategy":
    """Get the appropriate transformer for a bundle.

    Parameters
    ----------
    transform_type : str
        Transform type from registry entry.
    bundle_name : str
        Bundle name for special-case handling (e.g., "clinvar").

    Returns
    -------
    TransformStrategy
        Transformer instance with a transform() method.
    """
    # Special case: ClinVar needs its own transformer
    if bundle_name == "clinvar" and transform_type == "vcf_to_tsv":
        return ClinvarVcfTransformer()  # type: ignore[return-value,unused-ignore]

    cls = TRANSFORM_REGISTRY.get(transform_type)
    if cls is None:
        raise ValueError(f"Unknown transform type: {transform_type}")
    instance = cls()
    return instance  # type: ignore[no-any-return]
