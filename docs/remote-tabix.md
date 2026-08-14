# Remote Tabix Score Backend

Query CADD and gnomAD databases directly from their public HTTP servers using tabix byte-range requests. No local download, no disk quota issues, no rate limits.

## Quick Start

```bash
# CADD scores via remote tabix (named preset)
vartriage --vcf panel.vcf --output results.json \
  --gene-annotation gencode.gtf \
  --cadd-remote cadd-v1.7-grch38

# gnomAD frequencies via remote tabix
vartriage --vcf panel.vcf --output results.json \
  --gene-annotation gencode.gtf \
  --gnomad-remote gnomad-exomes-v4-grch38

# Both together
vartriage --vcf patient.vcf.gz --output report.json \
  --gene-annotation gencode.gtf --clinvar clinvar.tsv \
  --cadd-remote cadd-v1.7-grch38 \
  --gnomad-remote gnomad-exomes-v4-grch38
```

## How It Works

Pathogenicity databases like CADD (80+ GB) and gnomAD (15+ GB per chromosome) publish their data as bgzipped, tabix-indexed files on public HTTP servers. The tabix index maps genomic coordinates to byte offsets in the compressed file.

`pysam.TabixFile` natively supports opening remote URLs. It downloads only the index (`.tbi`) and the specific compressed blocks needed for each query region via HTTP byte-range requests. A gene panel query fetching 500 scores might transfer a few hundred KB total instead of 80 GB.

The pipeline:

1. Resolves the preset name (or uses the raw URL)
2. Checks the local SQLite cache for each variant
3. Groups uncached variants by chromosome, then batches nearby variants (within 10 kb) into range queries
4. Queries the remote tabix file for each batch
5. Parses the returned records and matches to requested variants
6. Caches results locally for future runs

## Named Presets

List available presets:

```bash
vartriage remote list-presets
```

```
Name                           Source   Build    Description
--------------------------------------------------------------------------------
cadd-v1.7-grch37               cadd     grch37   CADD v1.7 all possible SNVs, GRCh37
cadd-v1.7-grch38               cadd     grch38   CADD v1.7 all possible SNVs, GRCh38
cadd-v1.7-indels-grch38        cadd     grch38   CADD v1.7 pre-scored gnomAD indels, GRCh38
gnomad-exomes-v4-grch38        gnomad   grch38   gnomAD exomes v4.1.1 per-chromosome VCF, GRCh38
```

Filter by source:

```bash
vartriage remote list-presets --source cadd
vartriage remote list-presets --source gnomad
```

You can also pass a full URL directly:

```bash
--cadd-remote "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"
```

## CLI Reference

| Flag | Description |
|------|-------------|
| `--cadd-remote <preset-or-url>` | Remote CADD score source. Ignored when `--cadd-scores` is set (local takes priority). |
| `--gnomad-remote <preset-or-url>` | Remote gnomAD frequency source. Ignored when `--gnomad` is set. |
| `--remote-cache-ttl <days>` | Cache TTL in days (default: 30). Use `-1` for pinned mode (never expire). |

### Priority Rules

When multiple score sources are available:

1. **Local file** (`--cadd-scores /path/to/file.tsv`) — always preferred
2. **Remote tabix** (`--cadd-remote cadd-v1.7-grch38`) — used when no local file
3. **REST API** (`--mode api`) — lowest priority for CADD

If both `--cadd-scores` and `--cadd-remote` are specified, the local file wins and a log message notes the override.

## Caching

Scores are cached in a SQLite database at `~/.vartriage/remote_cache.db`.

**Default TTL:** 30 days. CADD scores don't change between version releases, so a month of caching is safe.

**Pinned mode:** `--remote-cache-ttl -1` disables expiry entirely. Use this in clinical settings where bit-identical reproducibility across runs is required.

**Cache schema:** `(source, chrom, pos, ref, alt)` compound primary key with `score` and `fetched_at` columns.

**Warm cache performance:** A gene panel (500 variants) with warm cache completes in under 5 seconds. A chr22 WGS (42K variants) with warm cache completes in under 10 seconds.

## Circuit Breaker

Network failures are handled gracefully:

- After 5 consecutive failures within 60 seconds, the circuit breaker opens
- All remaining variants for that run are treated as "score unavailable"
- The pipeline continues (scoring is supplementary, not gating)
- A WARNING is logged when the breaker opens

The breaker resets on the next successful connection or after a 30-second recovery period.

## Python API

```python
from vartriage.remote.config import RemoteTabixConfig
from vartriage.remote.cadd import RemoteTabixCADD
from vartriage.remote.gnomad import RemoteTabixGnomAD

# Configure
config = RemoteTabixConfig(
    cadd_remote_url="cadd-v1.7-grch38",
    gnomad_remote_url="gnomad-exomes-v4-grch38",
    cache_ttl_days=30,
)

# CADD lookups
cadd = RemoteTabixCADD(config)
scores = cadd.lookup_batch([
    ("chr22", 28695868, "A", "G"),
    ("chr22", 28695869, "C", "T"),
])
# Returns: {("chr22", 28695868, "A", "G"): 23.5, ...}
cadd.close()

# gnomAD lookups (satisfies FrequencyDatabase protocol)
gnomad = RemoteTabixGnomAD(config)
frequencies = gnomad.lookup_batch([
    ("chr22", 28695868, "A", "G"),
])
# Returns: [0.00032]  (positional list, None for not found)
gnomad.close()
```

### Presets API

```python
from vartriage.remote.presets import resolve_preset, list_presets

# Resolve a name to URL
url = resolve_preset("cadd-v1.7-grch38")

# List all presets
for preset in list_presets():
    print(f"{preset.name}: {preset.description}")

# Filter by source
cadd_presets = list_presets(source="cadd")
```

## Configuration File

Remote tabix settings can be persisted in `~/.vartriage/config.toml`:

```toml
[remote]
cadd_url = "cadd-v1.7-grch38"
gnomad_url = "gnomad-exomes-v4-grch38"
cache_ttl_days = 30
cache_path = "~/.vartriage/remote_cache.db"
connect_timeout = 10.0
read_timeout = 30.0
max_retries = 3
```

CLI flags override config file values.

## Performance

| Workload | Variants | Cold Cache | Warm Cache |
|----------|----------|------------|------------|
| Gene panel | 500 | < 30 seconds | < 5 seconds |
| chr22 WGS | 42,000 | < 10 minutes | < 10 seconds |

Cold cache times depend on network latency and server response time. The batch window optimization (grouping variants within 10 kb) reduces round-trips significantly for sorted VCF input.

**Memory:** Peak memory increase from remote tabix is under 10 MB regardless of variant count. Scores are fetched per-region (streaming) and cached on-disk.

## Limitations

- **Remote REVEL/SpliceAI:** No public tabix-indexed files exist for these databases. Use local files or API mode.
- **Offline use:** Remote tabix requires network access. Use `--remote-cache-ttl -1` to pin the cache, then subsequent offline runs use cached scores.
- **Concurrency:** Queries are serial with batched range optimization. Parallel queries add complexity without proportional gain since HTTP pipelining is more effective.
- **Custom self-hosted files:** Any URL that pysam can open works automatically. Pass the full URL to `--cadd-remote` or `--gnomad-remote`.

## Troubleshooting

**"Unknown preset" error:**
Check available presets with `vartriage remote list-presets`. Preset names are case-sensitive.

**Slow first run:**
Cold cache queries are network-bound. Subsequent runs with warm cache are fast. Consider running a priming pass on a representative gene panel to populate the cache.

**Circuit breaker activating:**
Check network connectivity to the remote server. The CADD server (`krishna.gs.washington.edu`) and gnomAD S3 bucket must be reachable. Firewalls blocking HTTPS outbound on port 443 will trigger the breaker.

**Cache growing large:**
Run `vartriage remote list-presets` and check `~/.vartriage/remote_cache.db` size. Each cached score is ~100 bytes. A full WGS (5M variants) uses approximately 500 MB of cache.
