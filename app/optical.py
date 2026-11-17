"""Optical-duplicate detection for re-paired paired-end FASTQ results.

The checker reads the *published* R1/R2 result files of a successful re-pairing
task, verifies that they still satisfy the original pairing contract, and
clusters read pairs that are optical duplicates of each other.

Two read pairs are linked when they share the instrument, run, flowcell, lane,
tile and both read sequences, and the Euclidean distance between their pixel
coordinates is at most the requested threshold.  Clusters are the connected
components of that graph.  Every cluster keeps the pair with the highest
combined Phred+33 quality score (ties: the pair appearing earlier in R1).

The implementation never writes to the result files.  It holds only compact
per-pair metadata (no sequences or quality strings), so 50,000 pairs of up to
500 bases per end stay comfortably inside the 256 MiB process limit.  Any
missing file, malformed record, invalid read id/coordinate or out-of-range
quality byte aborts the whole check before a report is produced.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from . import config

# Read ids look like INST:RUN:FLOWCELL:LANE:TILE:X:Y.  Keep the scratch data
# small even if a malformed file presents an enormous first header field.
MAX_READ_ID_BYTES = 512
# A single FASTQ line can never legitimately exceed this for our inputs; the
# cap also keeps a corrupt, newline-free file from exhausting memory.
MAX_LINE_BYTES = 65_536


class OpticalCheckError(ValueError):
    """Raised when the result files cannot be analysed as a whole."""


@dataclass(slots=True)
class OpticalPair:
    read_id: str
    x: int
    y: int
    quality_total: int


@dataclass(slots=True)
class _Record:
    read_id: str
    instrument: str
    run: str
    flowcell: str
    lane: int
    tile: int
    x: int
    y: int
    sequence: bytes
    quality_total: int


def _read_line(fh, label: str) -> bytes:
    """Read one bounded line; a too-long line is a corrupted record."""
    line = fh.readline(MAX_LINE_BYTES + 1)
    if len(line) > MAX_LINE_BYTES:
        raise OpticalCheckError(
            f"{label}: line exceeds {MAX_LINE_BYTES} bytes; record is corrupt"
        )
    return line


def _strip_eol(line: bytes) -> bytes:
    return line.rstrip(b"\r\n")


def normalize_read_id(header_line: bytes, side_label: str) -> str:
    """Apply the service's pairing normalization to a FASTQ header line."""
    header = _strip_eol(header_line)
    if not header.startswith(b"@"):
        raise OpticalCheckError(
            f"{side_label}: header does not start with '@' (corrupt record)"
        )
    if len(header) == 1:
        raise OpticalCheckError(f"{side_label}: empty read id")
    first_field = header[1:].split(None, 1)[0]
    if not first_field:
        raise OpticalCheckError(f"{side_label}: empty read id")
    if len(first_field) >= 2 and first_field[-2:-1] == b"/" and first_field[-1:] in (
        b"1",
        b"2",
    ):
        first_field = first_field[:-2]
    if not first_field:
        raise OpticalCheckError(f"{side_label}: empty read id after /1 or /2")
    if len(first_field) > MAX_READ_ID_BYTES:
        raise OpticalCheckError(
            f"{side_label}: read id exceeds {MAX_READ_ID_BYTES} bytes"
        )
    return first_field.decode("latin-1")


def _parse_coordinate_field(value: str, name: str, read_id: str) -> int:
    if not value or not value.isascii() or not value.isdigit():
        raise OpticalCheckError(
            f"invalid {name} coordinate in read id '{read_id}'"
        )
    return int(value)


def parse_read_id(read_id: str) -> _Record:
    """Parse INST:RUN:FLOWCELL:LANE:TILE:X:Y; invalid coordinates abort."""
    parts = read_id.split(":")
    if len(parts) != 7:
        raise OpticalCheckError(
            f"read id '{read_id}' does not match INST:RUN:FLOWCELL:LANE:TILE:X:Y"
        )
    instrument, run, flowcell, lane_s, tile_s, x_s, y_s = parts
    for name, value in (("instrument", instrument), ("run", run),
                        ("flowcell", flowcell)):
        if not value:
            raise OpticalCheckError(f"read id '{read_id}' has an empty {name} field")
    lane = _parse_coordinate_field(lane_s, "lane", read_id)
    tile = _parse_coordinate_field(tile_s, "tile", read_id)
    x = _parse_coordinate_field(x_s, "x", read_id)
    y = _parse_coordinate_field(y_s, "y", read_id)
    if lane < 1 or tile < 1:
        raise OpticalCheckError(
            f"invalid lane/tile coordinate in read id '{read_id}'"
        )
    # Sequence and quality are filled in by _record_from_raw; this helper is
    # only responsible for parsing and validating the identifier.
    return _Record(read_id, instrument, run, flowcell, lane, tile, x, y, b"", 0)


def _next_raw_record(
    fh, side_label: str
) -> tuple[bytes, bytes, bytes, bytes] | None:
    """Return (header, sequence, plus, quality) or None at clean EOF."""
    header = _read_line(fh, side_label)
    if not header:
        return None
    sequence = _read_line(fh, side_label)
    plus = _read_line(fh, side_label)
    quality = _read_line(fh, side_label)
    if not sequence or not plus or not quality:
        raise OpticalCheckError(
            f"{side_label}: truncated record (FASTQ records must span four lines)"
        )
    if not _strip_eol(header).startswith(b"@"):
        raise OpticalCheckError(
            f"{side_label}: header does not start with '@' (corrupt record)"
        )
    if not _strip_eol(plus).startswith(b"+"):
        raise OpticalCheckError(
            f"{side_label}: third line does not start with '+' (corrupt record)"
        )
    return header, sequence, plus, quality


def _record_from_raw(
    raw: tuple[bytes, bytes, bytes, bytes], side_label: str
) -> _Record:
    header, sequence, _plus, quality = raw
    read_id = normalize_read_id(header, side_label)
    rec = parse_read_id(read_id)

    seq = _strip_eol(sequence)
    qual = _strip_eol(quality)
    if not seq:
        raise OpticalCheckError(f"{side_label}: empty sequence in read '{read_id}'")
    if not qual:
        raise OpticalCheckError(
            f"{side_label}: empty quality string in read '{read_id}'"
        )
    if len(seq) != len(qual):
        raise OpticalCheckError(
            f"{side_label}: sequence length {len(seq)} != quality length "
            f"{len(qual)} in read '{read_id}'"
        )
    if len(seq) > config.OPTICAL_MAX_READ_LENGTH:
        raise OpticalCheckError(
            f"{side_label}: read '{read_id}' has {len(seq)} bases, exceeding the "
            f"{config.OPTICAL_MAX_READ_LENGTH}-base limit"
        )
    low = min(qual)
    high = max(qual)
    if low < 33 or high > 126:
        raise OpticalCheckError(
            f"{side_label}: quality byte outside ASCII 33..126 in read '{read_id}'"
        )
    rec.sequence = seq
    rec.quality_total = sum(qual) - 33 * len(qual)
    return rec


def _load_pairs(
    r1_path: Path, r2_path: Path
) -> tuple[list[OpticalPair], list[list[int]]]:
    """Parse both files in lockstep and group identical-key read pairs."""
    pairs: list[OpticalPair] = []
    groups: dict[tuple, int] = {}
    group_members: list[list[int]] = []
    seen_ids: set[str] = set()

    try:
        f1 = r1_path.open("rb")
    except OSError as exc:
        raise OpticalCheckError(f"R1 result file is missing or unreadable: {exc}")
    try:
        f2 = r2_path.open("rb")
    except OSError as exc:
        f1.close()
        raise OpticalCheckError(f"R2 result file is missing or unreadable: {exc}")

    with f1, f2:
        pair_no = 0
        while True:
            raw1 = _next_raw_record(f1, "R1")
            raw2 = _next_raw_record(f2, "R2")
            if raw1 is None and raw2 is None:
                break
            if raw1 is None or raw2 is None:
                raise OpticalCheckError(
                    f"R1 and R2 contain different numbers of records "
                    f"(mismatch at pair {pair_no + 1})"
                )

            # The published files must already be aligned pair by pair.
            rid1 = normalize_read_id(raw1[0], "R1")
            rid2 = normalize_read_id(raw2[0], "R2")
            if rid1 != rid2:
                raise OpticalCheckError(
                    f"R1/R2 pair mismatch at pair {pair_no + 1}: "
                    f"'{rid1}' != '{rid2}'"
                )
            if rid1 in seen_ids:
                raise OpticalCheckError(f"duplicate read id '{rid1}' in results")
            seen_ids.add(rid1)

            pair_no += 1
            if pair_no > config.OPTICAL_MAX_PAIRS:
                raise OpticalCheckError(
                    f"result contains more than {config.OPTICAL_MAX_PAIRS} read pairs"
                )

            rec1 = _record_from_raw(raw1, "R1")
            rec2 = _record_from_raw(raw2, "R2")
            # Both ends describe the same physical cluster.  The textual ids
            # already matched, so this is a defensive consistency check.
            if (
                rec1.instrument != rec2.instrument
                or rec1.run != rec2.run
                or rec1.flowcell != rec2.flowcell
                or rec1.lane != rec2.lane
                or rec1.tile != rec2.tile
                or rec1.x != rec2.x
                or rec1.y != rec2.y
            ):
                raise OpticalCheckError(
                    f"R1/R2 coordinate mismatch at pair {pair_no}: '{rid1}'"
                )

            key = (
                rec1.instrument,
                rec1.run,
                rec1.flowcell,
                rec1.lane,
                rec1.tile,
                hashlib.blake2b(rec1.sequence, digest_size=16).digest(),
                hashlib.blake2b(rec2.sequence, digest_size=16).digest(),
            )
            group_id = groups.get(key)
            if group_id is None:
                group_id = len(group_members)
                groups[key] = group_id
                group_members.append([])
            group_members[group_id].append(len(pairs))
            pairs.append(
                OpticalPair(
                    read_id=rid1,
                    x=rec1.x,
                    y=rec1.y,
                    quality_total=rec1.quality_total + rec2.quality_total,
                )
            )

    if not pairs:
        raise OpticalCheckError("result files contain no read pairs")
    return pairs, group_members


def _distance_offsets(threshold: int) -> list[tuple[int, int]]:
    """Half-plane integer offsets (dx > 0, or dx == 0 and dy > 0) within range.

    A link between p and q is found when the half-plane endpoint is visited
    and q = p + offset is looked up, so the reverse orientation is covered by
    the other endpoint of the same edge.
    """
    limit_sq = threshold * threshold
    offsets: list[tuple[int, int]] = []
    for dx in range(0, threshold + 1):
        for dy in range(-threshold, threshold + 1):
            if (dx > 0 or dy > 0) and dx * dx + dy * dy <= limit_sq:
                offsets.append((dx, dy))
    return offsets


class _UnionFind:
    __slots__ = ("parent", "size")

    def __init__(self, n: int) -> None:
        self.parent = list(range(n))
        self.size = [1] * n

    def find(self, value: int) -> int:
        parent = self.parent
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(self, a: int, b: int) -> None:
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return
        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]


def cluster_pairs(
    pairs: Sequence[OpticalPair],
    group_members: Sequence[Sequence[int]],
    threshold: int,
) -> tuple[int, list[dict]]:
    """Return (duplicate_count, duplicate entries in R1 order)."""
    n = len(pairs)
    uf = _UnionFind(n)
    offsets = _distance_offsets(threshold)

    for members in group_members:
        if len(members) < 2:
            continue
        # Every pair at exactly the same pixel coordinate is at distance 0.
        by_coord: dict[tuple[int, int], list[int]] = {}
        for idx in members:
            pair = pairs[idx]
            by_coord.setdefault((pair.x, pair.y), []).append(idx)
        for idxs in by_coord.values():
            base = idxs[0]
            for other in idxs[1:]:
                uf.union(base, other)
        # Then link distinct coordinates within the pixel-distance threshold.
        # The threshold is at most 10, so enumerating integer offsets is both
        # exact and cheap (at most ~157 lookups per distinct coordinate).
        for (x, y), idxs in by_coord.items():
            base = idxs[0]
            for dx, dy in offsets:
                neighbour = by_coord.get((x + dx, y + dy))
                if neighbour is not None:
                    uf.union(base, neighbour[0])

    # Pick each component's winner: highest total quality, ties keep the pair
    # with the lower R1 ordinal (strict > leaves the earlier pair in place).
    winners: dict[int, int] = {}
    for idx in range(n):
        root = uf.find(idx)
        current = winners.get(root)
        if current is None or pairs[idx].quality_total > pairs[current].quality_total:
            winners[root] = idx

    duplicates: list[dict] = []
    for idx in range(n):
        keep_idx = winners[uf.find(idx)]
        if keep_idx != idx:
            duplicates.append(
                {
                    "read_id": pairs[idx].read_id,
                    "keep_read_id": pairs[keep_idx].read_id,
                }
            )
    return len(duplicates), duplicates


def analyze(
    r1_path: Path, r2_path: Path, threshold: int
) -> dict:
    """Run the full check and return the report body (without task metadata)."""
    if not isinstance(threshold, int) or not (
        config.OPTICAL_MIN_DISTANCE <= threshold <= config.OPTICAL_MAX_DISTANCE
    ):
        raise OpticalCheckError(
            f"pixel distance threshold must be between "
            f"{config.OPTICAL_MIN_DISTANCE} and {config.OPTICAL_MAX_DISTANCE}"
        )
    pairs, group_members = _load_pairs(r1_path, r2_path)
    duplicate_count, duplicates = cluster_pairs(pairs, group_members, threshold)
    return {
        "total_pairs": len(pairs),
        "duplicate_pairs": duplicate_count,
        "duplicates": duplicates,
    }
