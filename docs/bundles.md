# Score Bundles

vartriage requires reference files (gnomAD, ClinVar, CADD, REVEL, SpliceAI, GENCODE) for annotation and classification. The bundle system automates downloading, transforming, and managing these files.

## Quick start

```bash
# List available bundles
vartriage bundle list

# Download ClinVar (smallest, ~80 MB)
vartriage bundle download --bundle clinvar

# Run pipeline using installed bundles
vartriage --vcf input.vcf.gz --output results.json --use-bundles
```

## Available bundles

| Bundle | Source | Size | Builds |
| -------- | -------- | ------ | -------- |
| `clinvar` | NCBI ClinVar | ~80 MB | GRCh37, GRCh38 |
| `gnomad-exomes-chr22` | gnomAD v4.1.1 | ~4.7 GB | GRCh38 |
| `revel` | REVEL v1.3 | ~2 GB | GRCh38 |
| `gencode` | GENCODE v46 | ~50 MB | GRCh38 |
| `spliceai` | SpliceAI (Illumina) | ~30 GB | GRCh37, GRCh38 |

## Commands

### `vartriage bundle download`

Download and transform a reference bundle.

```bash
vartriage bundle download --bundle clinvar --build grch38
vartriage bundle download --bundle gnomad-exomes-chr22
vartriage bundle download --bundle revel --no-progress
```

Options:

- `--bundle NAME` — bundle to download (required)
- `--build BUILD` — genome build, default from config (grch38)
- `--dest DIR` — custom storage directory
- `--no-transform` — download raw file only, skip TSV conversion
- `--no-progress` — suppress progress bar

### `vartriage bundle list`

Show available bundles and their installation status.

```bash
vartriage bundle list
vartriage bundle list --build grch37
vartriage bundle list --json
```

### `vartriage bundle verify`

Check integrity of installed bundles.

```bash
vartriage bundle verify
vartriage bundle verify --bundle clinvar
```

### `vartriage bundle status`

Show installed bundles, versions, and disk usage.

```bash
vartriage bundle status
```

### `vartriage bundle update-registry`

Check for registry updates (future: fetch from GitHub releases).

```bash
vartriage bundle update-registry
```

## Using bundles with the pipeline

Pass `--use-bundles` to auto-resolve reference paths from installed bundles:

```bash
vartriage \
  --vcf patient.vcf.gz \
  --output report.json \
  --use-bundles
```

Explicitly passed paths always take precedence over bundles:

```bash
# Uses bundle for gnomAD but a custom ClinVar file
vartriage \
  --vcf patient.vcf.gz \
  --output report.json \
  --use-bundles \
  --clinvar /path/to/my/clinvar.tsv
```

## Storage layout

Bundles are stored at `~/.vartriage/bundles/` by default:

```text
~/.vartriage/
  bundles/
    grch38/
      clinvar/
        raw/clinvar.vcf.gz    # Original download
        clinvar.tsv           # Transformed output
        manifest.json         # Version, checksums, timestamps
      gnomad-exomes-chr22/
        raw/gnomad.exomes.v4.1.1.sites.chr22.vcf.bgz
        gnomad-exomes-chr22.tsv
        manifest.json
```

Override with environment variable:

```bash
export VARTRIAGE_BUNDLE_STORAGE=/data/vartriage/bundles
```

## Configuration

Create `~/.vartriage/config.toml` for persistent settings:

```toml
[bundle]
default_build = "grch38"
download_concurrency = 2
storage_path = "~/.vartriage/bundles"
auto_verify = true

[bundle.proxy]
http_proxy = ""
https_proxy = ""
```

## Chromosome naming

The bundle transformers automatically handle chromosome name normalization:

- **ClinVar**: Adds `chr` prefix (ClinVar GRCh38 VCF uses bare names)
- **gnomAD**: Already uses `chr` prefix in GRCh38 builds
- **GENCODE**: Already uses `chr` prefix

No manual normalization needed when using bundles.

## Requirements

- `bcftools` (for VCF-to-TSV transformation) OR `pysam` (Python fallback)
- `wget` or `curl` (auto-detected for downloads)
- Sufficient disk space (check with `vartriage bundle list` before downloading)
