"""Self-tests for the optical-duplicate check API.

    .venv/bin/python -I tests/test_optical.py

Every test run gets an isolated FASTQ_DATA_DIR under tests/output.  Reports
are produced through the real HTTP surface (httpx ASGI transport), and the
result files are hashed before and after to prove they are never modified.
"""

from __future__ import annotations

import asyncio
import hashlib
import io
import json
import os
import shutil
import sys
import tempfile
import time
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import httpx  # noqa: E402

from app import config, db, worker  # noqa: E402
from app.main import app  # noqa: E402

RESULTS: list[tuple[str, bool, str]] = []


def record(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" - {detail}" if detail else ""))


def fq(records: list[tuple[str, bytes, bytes]], suffix: str = "") -> io.BytesIO:
    out = bytearray()
    for rid, seq, qual in records:
        out += b"@" + rid.encode() + suffix.encode() + b"\n"
        out += seq + b"\n+\n" + qual + b"\n"
    return io.BytesIO(bytes(out))


def rid(inst: str, run: str, fc: str, lane: int, tile: int, x: int, y: int) -> str:
    return f"{inst}:{run}:{fc}:{lane}:{tile}:{x}:{y}"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


async def wait_terminal(client: httpx.AsyncClient, task_id: str,
                        timeout: float = 60.0) -> dict:
    deadline = time.time() + timeout
    while time.time() < deadline:
        body = (await client.get(f"/tasks/{task_id}")).json()
        if body["status"] in ("succeeded", "failed", "cancelled"):
            return body
        await asyncio.sleep(0.05)
    raise AssertionError(f"task {task_id} did not finish in {timeout}s")


def publish_via_db(task_id: str, r1: bytes, r2: bytes) -> None:
    """Create a task that looks exactly like a finished re-pairing task."""
    base = config.TASKS_DIR / task_id / "output"
    base.mkdir(parents=True)
    (base / "R1.fastq").write_bytes(r1)
    (base / "R2.fastq").write_bytes(r2)
    conn = db.connect()
    try:
        db.create_task(conn, task_id, "r1.fastq", "r2.fastq", len(r1), len(r2))
        db.update_task(
            conn,
            task_id,
            status="succeeded",
            stage=config.STAGE_DONE,
            progress=1.0,
            message="paired",
            result_r1="output/R1.fastq",
            result_r2="output/R2.fastq",
        )
    finally:
        conn.close()


def make_task(r1_rows, r2_rows, task_id: str | None = None) -> str:
    task_id = task_id or uuid.uuid4().hex
    r1 = fq(r1_rows, "/1").getvalue()
    r2 = fq(r2_rows, "/2").getvalue()
    publish_via_db(task_id, r1, r2)
    return task_id


async def check(client: httpx.AsyncClient, task_id: str, threshold: int):
    return await client.post(
        "/optical-checks",
        json={"task_id": task_id, "distance_threshold": threshold},
    )


async def main_async() -> int:
    tmp_root = Path(tempfile.mkdtemp(prefix="optical-test-",
                                     dir=str(ROOT / "tests" / "output")))
    data_dir = tmp_root / "data"
    os.environ["FASTQ_DATA_DIR"] = str(data_dir)
    config.DATA_DIR = data_dir
    config.DB_PATH = data_dir / "tasks.db"
    config.TASKS_DIR = data_dir / "tasks"
    config.PROGRESS_UPDATE_EVERY_BYTES = 64 * 1024

    db._initialized = False
    db.init_db()
    worker.recover_on_startup()

    transport = httpx.ASGITransport(app=app)
    async with httpx.AsyncClient(transport=transport, base_url="http://test",
                                 timeout=120) as client:
        # ---------------------------------------------------------------- 1
        # Clustering: distance 1, distance exactly at the threshold, and a
        # pair one pixel beyond it.  Coordinates are Euclidean, not Chebyshev.
        seq1, seq2 = b"A" * 10, b"T" * 10
        inst = ("INST", "RUN", "FC")
        q_hi, q_mid, q_lo = b"J" * 10, b"I" * 10, b"@" * 10  # Q41/40/31
        rows = [
            (rid(*inst, 1, 100, 10, 10), seq1, q_lo),
            (rid(*inst, 1, 100, 10, 11), seq1, q_hi),   # distance 1
            (rid(*inst, 1, 100, 13, 14), seq1, q_mid),  # distance 5 from (10,10)
            (rid(*inst, 1, 100, 20, 20), seq1, q_mid),  # isolated
        ]
        rows2 = [(r[0], seq2, r[2]) for r in rows]
        tid = make_task(rows, rows2)
        p1 = config.TASKS_DIR / tid / "output" / "R1.fastq"
        p2 = config.TASKS_DIR / tid / "output" / "R2.fastq"
        h1, h2 = digest(p1), digest(p2)

        r = await check(client, tid, 4)
        body = r.json()
        record("check returns 201", r.status_code == 201, str(r.status_code))
        record("threshold 4 does not link distance 5",
               body["duplicate_pairs"] == 1, str(body["duplicate_pairs"]))
        r = await check(client, tid, 5)
        body = r.json()
        record("threshold 5 links Euclidean distance 5",
               body["duplicate_pairs"] == 2, str(body["duplicate_pairs"]))
        keep = rid(*inst, 1, 100, 10, 11)
        record("highest quality in the cluster is kept",
               [d["keep_read_id"] for d in body["duplicates"]]
               == [keep, keep], json.dumps(body["duplicates"]))
        record("duplicates are reported in R1 order",
               [d["read_id"] for d in body["duplicates"]]
               == [rows[0][0], rows[2][0]], json.dumps(body["duplicates"]))
        record("total_pairs is the number of read pairs",
               body["total_pairs"] == 4, str(body["total_pairs"]))
        record("result files are untouched", digest(p1) == h1 and digest(p2) == h2)
        report_id = body["id"]
        record("check is persisted",
               (await client.get(f"/optical-checks/{report_id}")).json()
               == body)

        # Same coordinate, same ids, but different sequences: never linked.
        rows_b = [(rid(*inst, 2, 7, 50, 50), b"AAAA", b"J" * 4),
                  (rid(*inst, 2, 7, 50, 51), b"CCCC", b"@" * 4)]
        tid_b = make_task(
            rows_b, [(r[0], b"T" * 4, r[2]) for r in rows_b]
        )
        body = (await check(client, tid_b, 10)).json()
        record("different R1 sequences are separate clusters",
               body["duplicate_pairs"] == 0, json.dumps(body["duplicates"]))

        # Same id fields but different flowcell/run/tile are not comparable.
        rows_c = [(rid("INST", "RUN", "FC", 1, 100, 5, 5), seq1, q_hi),
                  (rid("INST", "RUN", "FC2", 1, 100, 5, 6), seq1, q_lo),
                  (rid("INST", "RUN", "FC", 2, 100, 5, 6), seq1, q_lo),
                  (rid("INST", "RUN", "FC", 1, 101, 5, 6), seq1, q_lo),
                  (rid("OTHER", "RUN", "FC", 1, 100, 5, 6), seq1, q_lo)]
        tid_c = make_task(rows_c, [(r[0], seq2, r[2]) for r in rows_c])
        body = (await check(client, tid_c, 10)).json()
        record("different flowcell/lane/tile/instrument are separate clusters",
               body["duplicate_pairs"] == 0, json.dumps(body["duplicates"]))
        rows_d = [(rid("INST", "RUN2", "FC", 1, 100, 5, 5), seq1, q_hi),
                  (rid("INST", "RUN3", "FC", 1, 100, 5, 6), seq1, q_lo)]
        tid_d = make_task(rows_d, [(r[0], seq2, r[2]) for r in rows_d])
        record("different run ids are separate clusters",
               (await check(client, tid_d, 10)).json()["duplicate_pairs"] == 0)

        # ---------------------------------------------------------------- 2
        # Connected components and tie-breaking by R1 order.  Chain A(0,0) -
        # B(4,0) - C(8,0) is one cluster at threshold 4, although A and C are
        # 8 apart.  All quality totals are equal, so the winner is the
        # earliest pair in R1.
        chain = [(rid(*inst, 9, 60, 0, 0), seq1, q_hi),
                 (rid(*inst, 9, 60, 4, 0), seq1, q_hi),
                 (rid(*inst, 9, 60, 8, 0), seq1, q_hi)]
        tid_e = make_task(chain, [(r[0], seq2, r[2]) for r in chain])
        body = (await check(client, tid_e, 4)).json()
        record("clusters are connected components",
               body["duplicate_pairs"] == 2, json.dumps(body["duplicates"]))
        record("tie breaks to the earliest pair in R1",
               all(d["keep_read_id"] == chain[0][0] for d in body["duplicates"]),
               json.dumps(body["duplicates"]))

        # Combined quality of both ends decides the winner.
        mixed = [(rid(*inst, 5, 300, 0, 0), b"AAAA", b"I" * 4),   # 4*40
                 (rid(*inst, 5, 300, 1, 0), b"AAAA", b"!" * 4)]   # 0
        mixed2 = [(mixed[0][0], b"TTTT", b"!" * 4),               # 0 -> total 160
                  (mixed[1][0], b"TTTT", b"I" * 4)]               # 160 -> total 160
        tid_f = make_task(mixed, mixed2)
        body = (await check(client, tid_f, 2)).json()
        record("equal combined totals keep the earlier pair",
               body["duplicates"] == [
                   {"read_id": mixed[1][0], "keep_read_id": mixed[0][0]}],
               json.dumps(body["duplicates"]))

        # ---------------------------------------------------------------- 3
        # Whole-run failures: no partial report, a failed record is kept, and
        # the result files stay byte-identical.
        def raw_task(r1: bytes, r2: bytes) -> str:
            task_id = uuid.uuid4().hex
            publish_via_db(task_id, r1, r2)
            return task_id

        good1 = fq([(rid(*inst, 7, 1, 1, 1), b"AAAA", b"~~~~")], "/1").getvalue()
        good2 = fq([(rid(*inst, 7, 1, 1, 1), b"TTTT", b"~~~~")], "/2").getvalue()

        broken = {
            "quality byte 32": (good1.replace(b"~~~~", b"  ~ "),
                                good2, 32),
            "quality byte 127": (good1.replace(b"~~~~", b"\x7f\x7f\x7f\x7f"),
                                 good2, None),
            "invalid coordinate": (
                fq([(rid(*inst, 7, 1, 1, "x"), b"AAAA", b"~~~~")], "/1").getvalue(),
                good2, None),
            "missing id fields": (
                b"@INST:RUN:FC:1:1\nAAAA\n+\n~~~~\n", good2, None),
            "corrupt plus line": (
                good1.replace(b"\n+\n", b"\n-\n"), good2, None),
            "truncated record": (good1[:-2], good2, None),
            "missing quality line": (b"@INST:RUN:FC:1:1:1:1\nAAAA\n+\n", good2, None),
            "record count mismatch": (
                good1 + fq([(rid(*inst, 8, 1, 1, 1), b"AAAA", b"~~~~")], "/1").getvalue(),
                good2, None),
            "R1/R2 id mismatch": (
                good1,
                fq([(rid(*inst, 7, 1, 1, 2), b"TTTT", b"~~~~")], "/2").getvalue(),
                None),
        }
        for name, (r1_bytes, r2_bytes, _) in broken.items():
            sub = raw_task(r1_bytes, r2_bytes)
            before = digest(config.TASKS_DIR / sub / "output" / "R1.fastq")
            resp = await check(client, sub, 2)
            detail = resp.json().get("detail")
            record(f"{name} fails the whole check",
                   resp.status_code == 422, str(resp.status_code))
            record(f"{name} returns no partial list",
                   "duplicates" not in resp.json(), str(detail)[:80])
            after = digest(config.TASKS_DIR / sub / "output" / "R1.fastq")
            record(f"{name} leaves the result untouched", before == after)

        # failed checks are persisted with their error
        conn = db.connect()
        try:
            row = conn.execute(
                "SELECT state, error FROM optical_checks WHERE task_id = ? "
                "ORDER BY created_at DESC LIMIT 1", (sub,)
            ).fetchone()
        finally:
            conn.close()
        record("failed check is persisted with an error",
               row["state"] == "failed" and bool(row["error"]),
               f"{row['state']}: {row['error']}")

        # over-long reads are rejected rather than silently truncated
        long_seq = b"A" * (config.OPTICAL_MAX_READ_LENGTH + 1)
        same_id = rid(*inst, 1, 9, 1, 1)
        tid_g = raw_task(
            fq([(same_id, long_seq, b"~" * len(long_seq))], "/1").getvalue(),
            fq([(same_id, b"T" * 4, b"~~~~")], "/2").getvalue(),
        )
        resp = await check(client, tid_g, 2)
        record("reads longer than the limit are rejected",
               resp.status_code == 422 and "limit" in json.dumps(resp.json()),
               str(resp.status_code))

        # more pairs than the documented maximum
        old_max = config.OPTICAL_MAX_PAIRS
        config.OPTICAL_MAX_PAIRS = 2
        try:
            over = [(rid(*inst, 1, 44, 1, i), b"AAAA", b"~~~~") for i in range(3)]
            tid_h = make_task(over, [(r[0], b"TTTT", r[2]) for r in over])
            resp = await check(client, tid_h, 2)
            record("more than the maximum pair count is rejected",
                   resp.status_code == 422 and "more than" in json.dumps(resp.json()),
                   str(resp.status_code))
        finally:
            config.OPTICAL_MAX_PAIRS = old_max

        # ---------------------------------------------------------------- 4
        # Request validation and task-state guards.
        tid_ok = make_task(rows, rows2)
        for bad in (0, 11, -1):
            resp = await check(client, tid_ok, bad)
            record(f"threshold {bad} rejected", resp.status_code == 422,
                   str(resp.status_code))
        resp = await client.post("/optical-checks", json={"task_id": tid_ok})
        record("missing threshold rejected", resp.status_code == 422,
               str(resp.status_code))
        resp = await client.post("/optical-checks", json={"task_id": 123,
                                                          "distance_threshold": 2})
        record("non-string task id rejected", resp.status_code == 422,
               str(resp.status_code))
        resp = await check(client, uuid.uuid4().hex, 2)
        record("unknown task -> 404", resp.status_code == 404, str(resp.status_code))

        bad_task = uuid.uuid4().hex
        conn = db.connect()
        try:
            db.create_task(conn, bad_task, "a", "b", 1, 1)
            db.update_task(conn, bad_task, status="failed")
        finally:
            conn.close()
        resp = await check(client, bad_task, 2)
        record("failed task -> 409", resp.status_code == 409, str(resp.status_code))
        resp = await check(client, bad_task, 2)
        record("rejected task produces no report",
               "duplicates" not in resp.json())

        # A succeeded task whose result file vanished must fail, not回 empty.
        missing = make_task(rows, rows2)
        (config.TASKS_DIR / missing / "output" / "R2.fastq").unlink()
        resp = await check(client, missing, 2)
        record("missing result file -> 409", resp.status_code == 409,
               str(resp.status_code))

        # Concurrent checks are refused: two worst-case runs could exceed the
        # 256 MiB cap. Simulate a run in progress by holding the slot.
        import app.main as main_module

        other = make_task(rows, rows2)
        acquired = main_module._optical_slot.acquire(blocking=False)
        try:
            resp = await check(client, other, 2)
            record("a second concurrent check -> 429", resp.status_code == 429,
                   str(resp.status_code))
            record("busy check produces no partial list",
                   "duplicates" not in resp.json())
        finally:
            if acquired:
                main_module._optical_slot.release()
        resp = await check(client, other, 2)
        record("check works again after the slot is free",
               resp.status_code == 201, str(resp.status_code))

        resp = await client.get(f"/optical-checks/{uuid.uuid4().hex}")
        record("unknown check -> 404", resp.status_code == 404,
               str(resp.status_code))

        # A check left 'running' by a crash is marked failed on restart.
        conn = db.connect()
        try:
            db.create_optical_check(conn, "stuck-check", tid, 2)
        finally:
            conn.close()
        conn = db.connect()
        try:
            db.abort_stale_optical_checks(conn)
            row = db.get_optical_check(conn, "stuck-check")
        finally:
            conn.close()
        record("stale running check is failed on restart",
               row["state"] == "failed" and "restart" in row["error"],
               f"{row['state']}: {row['error']}")

        # ---------------------------------------------------------------- 5
        # End-to-end: upload shuffled FASTQ, re-pair, then check the result.
        ids = [rid(*inst, 3, 77, 1, 1), rid(*inst, 3, 77, 1, 2),
               rid(*inst, 3, 77, 1, 3)]
        up_r1 = fq([(i, seq1, q_lo) for i in ids], "/1")
        up_r2 = fq([(i, b"T" * 10, b"@" * 10) for i in reversed(ids)], "/2")
        resp = await client.post(
            "/tasks",
            files={"r1": ("R1.fastq", up_r1, "application/octet-stream"),
                   "r2": ("R2.fastq", up_r2, "application/octet-stream")},
        )
        task = await wait_terminal(client, resp.json()["id"])
        record("upload -> re-pair still succeeds", task["status"] == "succeeded",
               task.get("message", ""))
        body = (await check(client, task["id"], 2)).json()
        record("optical check on a real re-paired task",
               body["state"] == "ready" and body["total_pairs"] == 3
               and body["duplicate_pairs"] == 2
               and [d["keep_read_id"] for d in body["duplicates"]]
               == [ids[0], ids[0]],
               json.dumps({k: body[k] for k in ("state", "total_pairs",
                                               "duplicate_pairs",
                                               "duplicates")}))

        # The OpenAPI document exposes the new endpoints for /docs.
        openapi = (await client.get("/openapi.json")).json()
        record("openapi lists POST /optical-checks",
               "/optical-checks" in openapi["paths"])
        record("openapi lists GET /optical-checks/{check_id}",
               "/optical-checks/{check_id}" in openapi["paths"])
        record("docs page is served", (await client.get("/docs")).status_code == 200)

    shutil.rmtree(tmp_root, ignore_errors=True)
    failed = [name for name, ok, _ in RESULTS if not ok]
    print("\n%d/%d checks passed" % (len(RESULTS) - len(failed), len(RESULTS)))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(asyncio.run(main_async()))
