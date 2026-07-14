# API Annotation Mode

Query remote services (Ensembl VEP, ClinVar, CADD, SpliceAI) instead of maintaining local reference files. Useful for gene panels, exploratory triage, and environments where downloading 100+ GB of reference data is impractical.

## Installation

```bash
pip install vartriage[api]
```

This pulls in `httpx` for HTTP/2 connections with connection pooling.

## Quick Start

```bash
# Annotate a gene panel with zero local reference files
vartriage --vcf panel.vcf --output results.json --mode api

# Use an NCBI API key for higher ClinVar rate limits
vartriage --vcf panel.vcf --output results.json --mode api --api-key YOUR_KEY

# Hybrid: local gnomAD + API for everything else
vartriage --vcf panel.vcf --output results.json --mode hybrid --gnomad local_gnomad.tsv
```

## Modes

| Mode | Behavior |
|------|----------|
| `local` (default) | File-based backends only. Requires local reference files. |
| `api` | All annotation via remote APIs. No local files needed (except REVEL for PP3). |
| `hybrid` | Local files where available, API fills in any missing sources. |

## What Each API Provides

| Service | Data | Rate Limit |
|---------|------|------------|
| Ensembl VEP | Consequence, gene name, gnomAD frequency, CADD (plugin) | 15 req/sec, 55K/day |
| ClinVar E-utilities | Clinical significance, review status | 10 req/sec (with API key) |
| CADD API | Phred score (fallback when VEP plugin has no data) | 2 req/sec |
| SpliceAI Lookup | Delta scores for splice-relevant variants | 5 req/min |

## REVEL Limitation

REVEL has no public API. In pure API mode, PP3 evidence from REVEL is unavailable. To get REVEL-based PP3:

```bash
vartriage --vcf panel.vcf --output results.json --mode api --revel-scores revel_scores.tsv
```

## Response Caching

All API responses are cached locally in SQLite at `~/.vartriage/api_cache.db`. Cached variants skip the network entirely on subsequent runs.

Default TTL: 7 days for VEP/CADD/SpliceAI, 30 days for ClinVar.

### Cache Management

```bash
# Show cache statistics
vartriage api cache stats

# Clear all cached responses
vartriage api cache clear
```

### Pinned Cache (Clinical Reproducibility)

Set `cache_ttl_days = -1` in config to disable expiry. Cached responses persist indefinitely, producing bit-identical results across runs.

```toml
# ~/.vartriage/config.toml
[api]
cache_ttl_days = -1
```

## Configuration

Settings load from (highest priority wins):

1. CLI flags (`--mode`, `--api-key`)
2. Environment variables (`NCBI_API_KEY`, `HTTPS_PROXY`)
3. TOML config file (`~/.vartriage/config.toml`)
4. Built-in defaults

### TOML Configuration

```toml
# ~/.vartriage/config.toml
[api]
mode = "api"
genome_build = "grch38"
ncbi_api_key = ""
cache_ttl_days = 7
vep_batch_size = 200

[api.rate_limits]
vep_requests_per_second = 15
clinvar_requests_per_second = 10
cadd_requests_per_second = 2
spliceai_requests_per_minute = 5

[api.timeouts]
connect_seconds = 10
read_seconds = 30

[api.proxy]
url = "http://proxy.hospital.internal:8080"
```

### Environment Variables

| Variable | Purpose |
|----------|---------|
| `NCBI_API_KEY` | ClinVar rate limit upgrade (3 req/sec to 10 req/sec) |
| `HTTPS_PROXY` | HTTP proxy for outbound connections |
| `HTTP_PROXY` | Fallback proxy (HTTPS_PROXY takes priority) |

## Performance Expectations

| Workload | Approximate Time |
|----------|-----------------|
| Gene panel (50 variants) | < 30 seconds |
| Targeted sequencing (500 variants) | 2-3 minutes |
| Exome subset (5,000 variants) | ~15 minutes |
| Cached re-run (any size) | < 5 seconds |

API mode is 20-130x slower than local mode. For whole-exome or whole-genome scale, use local mode with the bundle downloader:

```bash
vartriage bundle download --bundle all
vartriage --vcf wgs.vcf --output results.json --use-bundles
```

## Error Handling

- **Rate limiting**: the client respects each service's limits automatically and waits when throttled.
- **Circuit breaker**: after 5 consecutive failures to a service, that source is skipped for 60 seconds. Remaining sources continue.
- **Retry**: transient errors (5xx, timeouts) are retried 3 times with exponential backoff (1s, 2s, 4s).
- **Graceful degradation**: if a service is unreachable, variants are annotated with available sources. Missing sources appear in the clinical report limitations section.

## Output Differences vs Local Mode

API mode uses Ensembl VEP for consequence calling, which is more granular than the local positional heuristic. Expect:

- ~10% of variants may get a different consequence (VEP uses full codon-level analysis)
- Gene names may differ slightly (VEP picks the canonical transcript differently)
- Allele frequencies match >99% (same gnomAD source, possible version skew)
- ClinVar assertions match >99% (live query vs monthly snapshot)

These differences mean API mode is generally **more accurate** for consequence calling. For validation studies comparing modes, use `--strict-concordance` in hybrid mode to log per-variant differences.

## Genome Build Support

Both GRCh38 (default) and GRCh37 are supported:

```bash
vartriage --vcf hg19_data.vcf --output results.json --mode api --genome-build grch37
```

GRCh37 queries use `grch37.rest.ensembl.org` instead of the default endpoint.
