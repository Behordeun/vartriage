"""CLI handlers for the 'vartriage bundle' subcommand group.

Implements: download, list, verify, status, update-registry.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from vartriage.bundle._checksums import compute_sha256, verify_checksum
from vartriage.bundle._disk import format_bytes
from vartriage.bundle.config import BundleConfig
from vartriage.bundle.downloader import BundleDownloader, DownloadError
from vartriage.bundle.manifest import BundleManifest
from vartriage.bundle.registry import BundleRegistry
from vartriage.bundle.storage import BundleStorage
from vartriage.bundle.transformer import get_transformer


def add_bundle_subcommands(subparsers: Any) -> None:
    """Add bundle subcommands to an argparse subparsers group.

    This is called from the main CLI module during parser setup.
    """
    dl = subparsers.add_parser("download", help="Download a reference bundle")
    dl.add_argument("--bundle", required=True, help="Bundle name to download")
    dl.add_argument("--build", default=None, help="Genome build (default: from config)")
    dl.add_argument("--dest", default=None, help="Custom storage directory")
    dl.add_argument("--no-transform", action="store_true", help="Skip transformation")
    dl.add_argument("--no-progress", action="store_true", help="Suppress progress bar")

    ls = subparsers.add_parser("list", help="List available and installed bundles")
    ls.add_argument("--build", default=None, help="Filter by genome build")
    ls.add_argument(
        "--json", action="store_true", dest="json_output", help="JSON output"
    )

    vf = subparsers.add_parser("verify", help="Verify installed bundle checksums")
    vf.add_argument("--bundle", default=None, help="Specific bundle to verify")
    vf.add_argument("--build", default=None, help="Genome build")

    subparsers.add_parser("status", help="Show installed bundles and disk usage")
    subparsers.add_parser("update-registry", help="Fetch latest bundle registry")


def run_bundle_command(args: argparse.Namespace) -> int:
    """Dispatch to the appropriate bundle subcommand handler.

    Returns
    -------
    int
        Exit code (0 for success, non-zero for errors).
    """
    config = BundleConfig.load()
    build = getattr(args, "build", None) or config.default_build

    command = getattr(args, "bundle_command", None)

    if command == "download":
        return _cmd_download(args, config, build)
    elif command == "list":
        return _cmd_list(args, config, build)
    elif command == "verify":
        return _cmd_verify(args, config, build)
    elif command == "status":
        _cmd_status(config, build)
        return 0
    elif command == "update-registry":
        _cmd_update_registry()
        return 0
    else:
        print("Usage: vartriage bundle <command>", file=sys.stderr)
        print(
            "Commands: download, list, verify, status, update-registry", file=sys.stderr
        )
        return 1


def _cmd_download(args: argparse.Namespace, config: BundleConfig, build: str) -> int:
    """Handle 'vartriage bundle download'."""
    registry = BundleRegistry.load()
    bundle_name = args.bundle

    entry = registry.get(bundle_name, build)
    if entry is None:
        print(
            f"Error: Bundle '{bundle_name}' not found for build '{build}'.",
            file=sys.stderr,
        )
        print(
            f"Available bundles: {', '.join(registry.bundle_names())}", file=sys.stderr
        )
        return 1

    # Resolve URL for this build
    url = entry.source_urls.get(build.lower())
    if not url:
        print(
            f"Error: No URL for build '{build}' in bundle '{bundle_name}'.",
            file=sys.stderr,
        )
        return 1

    # Set up storage
    dest_path = Path(args.dest) if args.dest else config.storage_path
    storage = BundleStorage(dest_path)
    storage.ensure_dirs(build, bundle_name)

    # Determine raw file destination
    raw_dir = storage.raw_dir(build, bundle_name)
    raw_filename = url.split("/")[-1]
    raw_dest = raw_dir / raw_filename

    print(f"Downloading {entry.display_name} v{entry.version} ({build})...")
    print(f"  Source: {url}")
    print(f"  Destination: {raw_dest}")

    # Download
    show_progress = not getattr(args, "no_progress", False)
    downloader = BundleDownloader(show_progress=show_progress)

    try:
        result = downloader.download(
            url=url,
            dest=raw_dest,
            expected_size=entry.expected_size_bytes,
            expected_checksum=entry.checksum_sha256,
        )
    except DownloadError as exc:
        print(f"\nError: {exc}", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"\nError: {exc}", file=sys.stderr)
        return 1

    print(
        f"\n  Downloaded: {format_bytes(result.bytes_downloaded)} "
        f"in {result.duration_seconds:.1f}s"
    )

    # Transform (unless --no-transform)
    if getattr(args, "no_transform", False):
        transformed_filename = raw_filename
    else:
        print(f"  Transforming ({entry.transform_type})...")
        transformer = get_transformer(entry.transform_type, bundle_name)
        bundle_dir = storage.bundle_dir(build, bundle_name)
        output_name = f"{bundle_name}.tsv"
        if entry.transform_type == "none":
            # For passthrough, keep original extension
            output_name = raw_filename.replace(".gz", "")

        transform_dest = bundle_dir / output_name
        try:
            t_result = transformer.transform(raw_dest, transform_dest, build)
            print(f"  Transformed: {t_result.rows_written} records -> {output_name}")
            transformed_filename = output_name
        except (ImportError, OSError, ValueError) as exc:
            print(f"  Transform failed: {exc}", file=sys.stderr)
            print(f"  Raw file retained at: {raw_dest}", file=sys.stderr)
            transformed_filename = raw_filename

    # Write manifest
    timestamp = datetime.now(timezone.utc).isoformat()
    manifest = BundleManifest(
        bundle_name=bundle_name,
        version=entry.version,
        genome_build=build,
        download_timestamp=timestamp,
        source_url=url,
        raw_checksum=compute_sha256(raw_dest) if raw_dest.exists() else "",
        transformed_checksum="",
        transformed_filename=transformed_filename,
        raw_filename=raw_filename,
        transform_type=entry.transform_type,
        file_size_bytes=raw_dest.stat().st_size if raw_dest.exists() else 0,
    )
    manifest.save(storage.manifest_path(build, bundle_name))
    print("  Manifest written.")
    print(f"\nDone. Bundle '{bundle_name}' installed for {build}.")
    return 0


def _cmd_list(args: argparse.Namespace, config: BundleConfig, build: str) -> int:
    """Handle 'vartriage bundle list'."""
    registry = BundleRegistry.load()
    storage = BundleStorage(config.storage_path)

    available = registry.available_for_build(build)

    if getattr(args, "json_output", False):
        _print_list_json(available, storage, build)
    else:
        _print_list_table(available, storage, build)

    return 0


def _print_list_json(available: list[Any], storage: BundleStorage, build: str) -> None:
    """Print bundle list as JSON."""
    output = []
    for entry in available:
        installed_ver = storage.installed_version(build, entry.name)
        output.append(
            {
                "name": entry.name,
                "display_name": entry.display_name,
                "version": entry.version,
                "installed_version": installed_ver,
                "status": "installed" if installed_ver else "available",
                "genome_builds": entry.genome_builds,
            }
        )
    print(json.dumps(output, indent=2))


def _print_list_table(available: list[Any], storage: BundleStorage, build: str) -> None:
    """Print bundle list as formatted table."""
    print(f"Score bundles for {build}:")
    print(f"{'Name':<22} {'Version':<10} {'Status':<12} {'Description'}")
    print("-" * 78)

    for entry in sorted(available, key=lambda e: e.name):
        status = _install_status(storage, build, entry)
        desc = (
            entry.description[:35] + "..."
            if len(entry.description) > 38
            else entry.description
        )
        print(f"  {entry.name:<20} {entry.version:<10} {status:<12} {desc}")

    print(f"\n{len(available)} bundles available for {build}.")


def _install_status(storage: BundleStorage, build: str, entry: object) -> str:
    """Determine install status string for a bundle entry."""
    installed_ver = storage.installed_version(build, entry.name)  # type: ignore[attr-defined]
    if not installed_ver:
        return "available"
    if installed_ver == entry.version:  # type: ignore[attr-defined]
        return "installed"
    return f"outdated ({installed_ver})"


def _cmd_verify(args: argparse.Namespace, config: BundleConfig, build: str) -> int:
    """Handle 'vartriage bundle verify'."""
    storage = BundleStorage(config.storage_path)
    manifests = storage.list_installed(build)

    if not manifests:
        print(f"No bundles installed for {build}.")
        return 0

    bundle_filter = getattr(args, "bundle", None)
    if bundle_filter:
        manifests = [m for m in manifests if m.bundle_name == bundle_filter]
        if not manifests:
            print(
                f"Bundle '{bundle_filter}' not installed for {build}.", file=sys.stderr
            )
            return 1

    errors = 0
    for manifest in manifests:
        if not _verify_manifest(manifest, storage, build):
            errors += 1

    if errors:
        print(f"\n{errors} bundle(s) failed verification.")
        return 1

    print(f"\nAll {len(manifests)} bundle(s) verified.")
    return 0


def _verify_manifest(
    manifest: BundleManifest, storage: BundleStorage, build: str
) -> bool:
    """Verify a single bundle manifest. Returns True if OK, False on failure."""
    bundle_dir = storage.bundle_dir(build, manifest.bundle_name)
    transformed = bundle_dir / manifest.transformed_filename

    if not transformed.exists():
        print(f"  {manifest.bundle_name}: MISSING ({transformed})")
        return False

    if manifest.raw_checksum:
        raw_path = storage.raw_dir(build, manifest.bundle_name) / manifest.raw_filename
        if raw_path.exists():
            actual = compute_sha256(raw_path)
            if actual != manifest.raw_checksum:
                print(f"  {manifest.bundle_name}: CHECKSUM MISMATCH (raw)")
                print(f"    expected: {manifest.raw_checksum}")
                print(f"    actual:   {actual}")
                return False

    print(f"  {manifest.bundle_name}: OK (v{manifest.version})")
    return True


def _cmd_status(config: BundleConfig, build: str) -> None:
    """Handle 'vartriage bundle status'."""
    storage = BundleStorage(config.storage_path)
    manifests = storage.list_installed(build)
    usage = storage.disk_usage(build)

    if not manifests:
        print(f"No bundles installed for {build}.")
        print("Run 'vartriage bundle download --bundle <name>' to get started.")
        return

    total_bytes = sum(usage.values())
    print(f"Installed bundles ({build}): {len(manifests)}")
    print(f"Total disk usage: {format_bytes(total_bytes)}")
    print()
    print(f"{'Bundle':<22} {'Version':<10} {'Size':<10} {'Downloaded'}")
    print("-" * 65)

    for manifest in manifests:
        size = usage.get(manifest.bundle_name, 0)
        date = (
            manifest.download_timestamp[:10]
            if manifest.download_timestamp
            else "unknown"
        )
        print(
            f"  {manifest.bundle_name:<20} {manifest.version:<10} "
            f"{format_bytes(size):<10} {date}"
        )

    print(f"\nStorage path: {config.storage_path}")


def _cmd_update_registry() -> None:
    """Handle 'vartriage bundle update-registry'."""
    registry = BundleRegistry.load()
    print(f"Registry version: {registry.version}")
    print(f"Last updated: {registry.updated_at}")
    print(f"Bundles: {len(registry.bundles)}")
    print()
    print("Note: In v0.6.0, the registry is bundled with the package.")
    print("Remote registry updates will be available in a future release.")
