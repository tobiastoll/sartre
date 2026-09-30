#!/usr/bin/env python3
"""

getSartreTables.py

Copyright (C) 2026 Tobias Toll and Thomas Ullrich

This file is part of the Sartre event generator.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation. See <http://www.gnu.org/licenses/>.

---------

Download and install a compatible Sartre table bundle.

The public table index is intentionally kept outside the source repository:
the ROOT tables are too large for normal source-control distribution.

---------

$Date: 2026-09-18 13:46:08 -0400 (Fri, 18 Sep 2026) $
$Author: ullrich $
"""

from __future__ import annotations

import argparse
import csv
import io
import os
from pathlib import Path
import shutil
import sys
import tarfile
import tempfile
from typing import Iterable
from urllib import error, parse, request


TABLE_SERVER = "https://rhig.physics.yale.edu/~ullrich/sartre-tables/"
INDEX_URL = parse.urljoin(TABLE_SERVER, "tableIndex.csv")


class TableDownloadError(RuntimeError):
    """An expected error while locating, downloading, or installing tables."""


def normalise(value: str) -> str:
    """Make user input and CSV values comparable without case or punctuation."""
    return "".join(character for character in value.casefold() if character.isalnum())


def normalise_particle(value: str) -> str:
    """Accept common particle-name and PDG-ID spellings in the table index."""
    aliases = {
        "22": "dvcs",
        "113": "rho",
        "333": "phi",
        "443": "jpsi",
        "553": "upsilon",
        "ups": "upsilon",
        "jpsi": "jpsi",
        "inclusive": "x",
    }
    canonical = normalise(value)
    return aliases.get(canonical, canonical)


def is_true(value: str) -> bool:
    """Interpret the isUPC column of the public CSV index."""
    return normalise(value) in {"1", "true", "yes", "y", "upc"}


def require_sartre_directory() -> Path:
    """Return SARTRE_DIR, or infer it from the installed script location."""
    sartre_dir = os.environ.get("SARTRE_DIR")
    if sartre_dir:
        installation = Path(sartre_dir).expanduser().resolve()
        source = f"SARTRE_DIR={installation}"
    else:
        # The installed path is normally /usr/local/sartre/bin/getSartreTables.py.
        # This fallback also works when sudo removes SARTRE_DIR from the environment.
        installation = Path(__file__).resolve().parent.parent
        source = f"the script location ({Path(__file__).resolve()})"

    tables = installation / "tables"
    if not installation.is_dir() or not tables.is_dir():
        if sartre_dir:
            raise TableDownloadError(
                f"{source} is not a Sartre installation with a tables/ directory."
            )
        raise TableDownloadError(
            f"Could not locate a Sartre installation from {source}. "
            "Set SARTRE_DIR to the Sartre installation directory."
        )
    return installation


def fetch_index() -> list[dict[str, str]]:
    """Fetch and validate the public CSV table index."""
    try:
        http_request = request.Request(INDEX_URL, headers={"User-Agent": "Sartre-table-downloader/1.0"})
        with request.urlopen(http_request, timeout=30) as response:
            contents = response.read().decode("utf-8-sig")
    except (error.HTTPError, error.URLError, TimeoutError, UnicodeDecodeError) as exc:
        raise TableDownloadError(f"Could not read table index {INDEX_URL}: {exc}") from exc

    reader = csv.DictReader(io.StringIO(contents))
    if not reader.fieldnames:
        raise TableDownloadError("The table index is empty or has no header row.")

    headers = {normalise(header): header for header in reader.fieldnames if header}
    required_headers = {
        "system": "system",
        "satmode": "sat mode",
        "particle": "particle",
        "parameterset": "parameter set",
        "isupc": "isUPC",
    }
    missing = [display for key, display in required_headers.items() if key not in headers]
    if missing:
        raise TableDownloadError(
            "The table index is missing required column(s): " + ", ".join(missing) + "."
        )

    archive_header = reader.fieldnames[-1]
    if not archive_header or archive_header in required_headers.values():
        raise TableDownloadError("The final table-index column must contain the tar.gz filename.")

    records: list[dict[str, str]] = []
    for line_number, row in enumerate(reader, start=2):
        record = {
            "system": (row.get(headers["system"]) or "").strip(),
            "sat_mode": (row.get(headers["satmode"]) or "").strip(),
            "particle": (row.get(headers["particle"]) or "").strip(),
            "parameter_set": (row.get(headers["parameterset"]) or "").strip(),
            "is_upc": (row.get(headers["isupc"]) or "").strip(),
            "archive": (row.get(archive_header) or "").strip(),
        }
        if not all(record.values()):
            raise TableDownloadError(f"Table-index row {line_number} has an empty required value.")
        if Path(record["archive"]).name != record["archive"] or not record["archive"].endswith(".tar.gz"):
            raise TableDownloadError(
                f"Table-index row {line_number} has an invalid tar.gz filename: {record['archive']!r}."
            )
        records.append(record)

    if not records:
        raise TableDownloadError("The table index contains no table sets.")
    return records


def print_index(records: Iterable[dict[str, str]], verbose: bool = False) -> None:
    """Print the public catalog, adding archive names only in verbose mode."""
    rows = [
        (
            record["system"],
            record["sat_mode"],
            record["particle"],
            record["parameter_set"],
            "yes" if is_true(record["is_upc"]) else "no",
            record["archive"],
        )
        for record in records
    ]
    headings = ("System", "Sat mode", "Particle", "Parameter set", "UPC", "Archive")
    if not verbose:
        headings = headings[:-1]
        rows = [row[:-1] for row in rows]
    widths = [len(heading) for heading in headings]
    for row in rows:
        for index, value in enumerate(row):
            widths[index] = max(widths[index], len(value))

    def format_row(row: tuple[str, ...]) -> str:
        return "  ".join(value.ljust(widths[index]) for index, value in enumerate(row))

    print(format_row(headings))
    print(format_row(tuple("-" * width for width in widths)))
    for row in rows:
        print(format_row(row))


def matching_record(records: Iterable[dict[str, str]], arguments: argparse.Namespace) -> dict[str, str] | None:
    """Find the one catalog entry requested by the user."""
    requested_upc = bool(arguments.upc)
    requested_particle = normalise_particle(arguments.particle)
    for record in records:
        if (
            normalise(record["system"]) == normalise(arguments.system)
            and normalise(record["sat_mode"]) == normalise(arguments.sat_mode)
            and normalise_particle(record["particle"]) == requested_particle
            and normalise(record["parameter_set"]) == normalise(arguments.parameter_set)
            and is_true(record["is_upc"]) == requested_upc
        ):
            return record
    return None


def download_archive(filename: str, table_directory: Path) -> Path:
    """Download an archive to a private temporary file below tables/."""
    archive_url = parse.urljoin(TABLE_SERVER, parse.quote(filename))
    try:
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=".getSartreTables-", suffix=".tar.gz", dir=table_directory
        )
    except OSError as exc:
        raise TableDownloadError(
            f"Cannot write to {table_directory}: {exc}. "
            "Use sudo when installing tables into a system-wide Sartre installation."
        ) from exc
    temporary_path = Path(temporary_name)
    try:
        http_request = request.Request(archive_url, headers={"User-Agent": "Sartre-table-downloader/1.0"})
        with request.urlopen(http_request, timeout=120) as response, os.fdopen(descriptor, "wb") as output:
            shutil.copyfileobj(response, output)
    except (error.HTTPError, error.URLError, TimeoutError, OSError) as exc:
        temporary_path.unlink(missing_ok=True)
        raise TableDownloadError(f"Could not download {archive_url}: {exc}") from exc
    return temporary_path


def validate_archive_members(archive: tarfile.TarFile, destination: Path) -> list[tarfile.TarInfo]:
    """Reject path traversal, links, and non-file archive members before extraction."""
    destination_text = str(destination.resolve())
    members = archive.getmembers()
    if not members:
        raise TableDownloadError("The downloaded archive is empty.")

    for member in members:
        member_path = Path(member.name)
        extracted_path = (destination / member_path).resolve()
        if member_path.is_absolute() or ".." in member_path.parts:
            raise TableDownloadError(f"Archive contains an unsafe path: {member.name!r}.")
        try:
            common_path = os.path.commonpath((destination_text, str(extracted_path)))
        except ValueError as exc:
            raise TableDownloadError(f"Archive contains an unsafe path: {member.name!r}.") from exc
        if common_path != destination_text:
            raise TableDownloadError(f"Archive contains an unsafe path: {member.name!r}.")
        if not (member.isdir() or member.isfile()):
            raise TableDownloadError(
                f"Archive contains unsupported member {member.name!r}; symbolic and hard links are not allowed."
            )
    return members


def install_archive(archive_path: Path, table_directory: Path) -> int:
    """Safely unpack a new archive without overwriting existing table files."""
    try:
        staging_directory = Path(tempfile.mkdtemp(prefix=".getSartreTables-", dir=table_directory))
    except OSError as exc:
        raise TableDownloadError(
            f"Cannot write to {table_directory}: {exc}. "
            "Use sudo when installing tables into a system-wide Sartre installation."
        ) from exc
    try:
        try:
            with tarfile.open(archive_path, mode="r:gz") as archive:
                members = validate_archive_members(archive, staging_directory)
                # Python 3.12+ requires an explicit extraction policy before
                # its 3.14 default changes.  The data filter complements the
                # path and member-type checks above.
                if sys.version_info >= (3, 12):
                    archive.extractall(staging_directory, members=members, filter="data")
                else:
                    archive.extractall(staging_directory, members=members)
        except (tarfile.TarError, OSError) as exc:
            raise TableDownloadError(f"Downloaded archive is not a usable tar.gz file: {exc}") from exc

        files = [path for path in staging_directory.rglob("*") if path.is_file()]
        if not files:
            raise TableDownloadError("The archive does not contain any table files.")

        collisions = [
            table_directory / path.relative_to(staging_directory)
            for path in files
            if (table_directory / path.relative_to(staging_directory)).exists()
        ]
        if collisions:
            examples = ", ".join(str(path.relative_to(table_directory)) for path in collisions[:3])
            suffix = " ..." if len(collisions) > 3 else ""
            raise TableDownloadError(
                "Refusing to overwrite existing table file(s): " + examples + suffix +
                ". Remove or move the conflicting set before retrying."
            )

        for source in files:
            destination = table_directory / source.relative_to(staging_directory)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(source), str(destination))
        return len(files)
    finally:
        shutil.rmtree(staging_directory, ignore_errors=True)


def rewrite_bare_options(arguments: list[str]) -> list[str]:
    """Also accept the documented convenience spelling part=phi and set=KMW."""
    option_names = {
        "sys": "--sys",
        "system": "--sys",
        "sat": "--sat",
        "satmode": "--sat",
        "part": "--part",
        "particle": "--part",
        "set": "--set",
        "parameterset": "--set",
    }
    rewritten: list[str] = []
    for argument in arguments:
        if not argument.startswith("-") and "=" in argument:
            name, value = argument.split("=", 1)
            option = option_names.get(normalise(name))
            if option:
                rewritten.append(f"{option}={value}")
                continue
        rewritten.append(argument)
    return rewritten


def parse_arguments(arguments: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="List, download, and install published Sartre table sets.",
        epilog=(
            "Examples:\n"
            "  getSartreTables.py list\n"
            "  getSartreTables.py --sys=ePb --sat=bNonSat part=phi set=KMW\n"
            "  getSartreTables.py --sys=ePb --sat=bSat --part=X --set=KMW --upc"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("command", nargs="?", choices=("list",), help="print the published table catalog")
    parser.add_argument("--sys", "--system", dest="system", metavar="SYSTEM", help="collision system, e.g. ePb")
    parser.add_argument("--sat", "--sat-mode", dest="sat_mode", metavar="MODE", help="dipole model, e.g. bSat")
    parser.add_argument("--part", "--particle", dest="particle", metavar="PARTICLE", help="particle name, PDG ID, or X for inclusive tables")
    parser.add_argument("--set", "--parameter-set", dest="parameter_set", metavar="SET", help="dipole-model parameter set, e.g. KMW")
    parser.add_argument("--upc", action="store_true", help="select a UPC table set (default: non-UPC)")
    parser.add_argument("-v", "--verbose", action="store_true", help="show archive filenames when listing table sets")
    parser.add_argument("-?", action="help", help="show this help message and exit")
    return parser.parse_args(rewrite_bare_options(arguments))


def main(arguments: list[str] | None = None) -> int:
    options = parse_arguments(sys.argv[1:] if arguments is None else arguments)
    try:
        installation = require_sartre_directory()
        table_directory = installation / "tables"
        records = fetch_index()

        if options.command == "list":
            print_index(records, verbose=options.verbose)
            return 0

        missing = [
            option
            for option, value in (
                ("--sys", options.system),
                ("--sat", options.sat_mode),
                ("--part", options.particle),
                ("--set", options.parameter_set),
            )
            if not value
        ]
        if missing:
            raise TableDownloadError(
                "A table download requires " + ", ".join(missing) + ". Use --help for examples."
            )

        record = matching_record(records, options)
        if record is None:
            mode = "UPC" if options.upc else "non-UPC"
            raise TableDownloadError(
                "No published table set matches "
                f"system={options.system}, sat mode={options.sat_mode}, "
                f"particle={options.particle}, parameter set={options.parameter_set}, {mode}. "
                "Run 'getSartreTables.py list' to see the available table sets."
            )

        print(f"Downloading {record['archive']} ...")
        archive_path = download_archive(record["archive"], table_directory)
        try:
            count = install_archive(archive_path, table_directory)
        finally:
            archive_path.unlink(missing_ok=True)
        print(f"Installed {count} table file(s) in {table_directory}.")
        print("The requested Sartre table set is ready to use.")
        return 0
    except TableDownloadError as exc:
        print(f"getSartreTables.py: error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
