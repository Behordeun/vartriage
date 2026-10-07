"""Confirmed-absent caching for the global-AF remote tabix path.

A variant gnomAD genuinely does not carry has no allele frequency. The global
path (lookup_batch -> remote_scores) must record that confirmed absence so a
repeat run is a cache hit returning None, rather than a fresh remote query every
time. The population path (lookup_batch_populations) already does this; these
tests hold the global-float path to the same contract.

The safety boundary mirrors the population path: only a query that actually
completed may be cached as absent. A fetch that exhausted its retries returns no
records, which is indistinguishable from a genuine absence by result alone;
caching it under clinical pinning (cache_ttl_days=-1) would permanently suppress
frequency evidence for a transiently-unreachable variant.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

from vartriage.remote.config import RemoteTabixConfig


def _make_config(tmp_path: Path) -> RemoteTabixConfig:
    return RemoteTabixConfig(
        gnomad_remote_url="https://example.com/{chrom}.vcf.bgz",
        cache_path=tmp_path / "absent_cache.db",
        cache_ttl_days=-1,
    )


class TestConfirmedAbsentGlobalCache:
    """lookup_batch persists a confirmed absence as a cache hit."""

    @patch("vartriage.remote.gnomad.pysam.TabixFile")
    def test_confirmed_absence_is_cached_not_requeried(
        self, mock_tabix_cls: MagicMock, tmp_path: Path
    ) -> None:
        from vartriage.remote.gnomad import RemoteTabixGnomAD

        config = _make_config(tmp_path)
        mock_tabix = MagicMock()
        # Successful query that returns no record for the variant: a genuine
        # absence, not a failure.
        mock_tabix.fetch.return_value = iter([])
        mock_tabix_cls.return_value = mock_tabix

        backend = RemoteTabixGnomAD(config)
        variant = ("chr22", 100, "A", "T")

        first = backend.lookup_batch([variant])
        assert first == [None]
        assert mock_tabix.fetch.call_count == 1

        # Second lookup must be served from the cache: the confirmed absence was
        # stored, so no second remote query happens.
        second = backend.lookup_batch([variant])
        assert second == [None]
        assert mock_tabix.fetch.call_count == 1, (
            "a confirmed absence must be served from cache, not re-queried"
        )
        backend.close()

    @patch("vartriage.remote.gnomad.pysam.TabixFile")
    def test_present_variant_still_caches_as_hit(
        self, mock_tabix_cls: MagicMock, tmp_path: Path
    ) -> None:
        from vartriage.remote.gnomad import RemoteTabixGnomAD

        config = _make_config(tmp_path)
        mock_tabix = MagicMock()
        mock_tabix.fetch.return_value = iter(
            ["chr22\t100\t.\tA\tT\t.\tPASS\tAF=0.0015;AN=100000"]
        )
        mock_tabix_cls.return_value = mock_tabix

        backend = RemoteTabixGnomAD(config)
        variant = ("chr22", 100, "A", "T")

        first = backend.lookup_batch([variant])
        assert first == [0.0015]

        # Present variants must remain cache hits (regression guard for the
        # existing score path).
        mock_tabix.fetch.return_value = iter([])
        second = backend.lookup_batch([variant])
        assert second == [0.0015]
        assert backend.cache_hits >= 1
        backend.close()

    @patch("vartriage.remote.gnomad.pysam.TabixFile")
    def test_mixed_present_and_absent_both_cache(
        self, mock_tabix_cls: MagicMock, tmp_path: Path
    ) -> None:
        from vartriage.remote.gnomad import RemoteTabixGnomAD

        config = _make_config(tmp_path)
        mock_tabix = MagicMock()
        # One variant present, one absent in the same window.
        mock_tabix.fetch.return_value = iter(
            ["chr22\t100\t.\tA\tT\t.\tPASS\tAF=0.0015;AN=100000"]
        )
        mock_tabix_cls.return_value = mock_tabix

        backend = RemoteTabixGnomAD(config)
        present = ("chr22", 100, "A", "T")
        absent = ("chr22", 200, "G", "C")

        first = backend.lookup_batch([present, absent])
        assert first == [0.0015, None]
        assert mock_tabix.fetch.call_count == 1

        # Both resolve from cache on the second pass: no further remote query.
        second = backend.lookup_batch([present, absent])
        assert second == [0.0015, None]
        assert mock_tabix.fetch.call_count == 1, (
            "present and absent variants must both be cached after one query"
        )
        backend.close()


class TestTransientFailureNotCachedAsAbsent:
    """A failed fetch must never be recorded as a confirmed absence."""

    @patch("vartriage.remote.gnomad.pysam.TabixFile")
    def test_failed_fetch_does_not_cache_absence(
        self, mock_tabix_cls: MagicMock, tmp_path: Path
    ) -> None:
        from vartriage.remote.gnomad import RemoteTabixGnomAD

        config = RemoteTabixConfig(
            gnomad_remote_url="https://example.com/{chrom}.vcf.bgz",
            cache_path=tmp_path / "transient_cache.db",
            cache_ttl_days=-1,
            max_retries=0,
        )
        mock_tabix = MagicMock()
        mock_tabix.fetch.side_effect = OSError("Connection refused")
        mock_tabix_cls.return_value = mock_tabix

        backend = RemoteTabixGnomAD(config)
        variant = ("chr22", 100, "A", "T")

        first = backend.lookup_batch([variant])
        assert first == [None]

        # The failure must NOT have been cached as absent. A later successful
        # query must re-fetch and return the real frequency, not a cached None.
        mock_tabix.fetch.side_effect = None
        mock_tabix.fetch.return_value = iter(
            ["chr22\t100\t.\tA\tT\t.\tPASS\tAF=0.004;AN=100000"]
        )
        second = backend.lookup_batch([variant])
        assert second == [0.004], (
            "a transiently-failed fetch must be re-queried, never cached as absent"
        )
        backend.close()
