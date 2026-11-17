"""FastAPI application: bounded uploads, status, cancellation, downloads."""

from __future__ import annotations

import json
import sqlite3
import tempfile
import threading
import uuid
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, File, HTTPException, UploadFile, status
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel, Field

from . import config, db, optical, worker
from .worker import UploadTooLarge, start_worker, task_dir


@asynccontextmanager
async def lifespan(_: FastAPI):
    db.init_db()
    worker.recover_on_startup()
    conn = db.connect()
    try:
        db.abort_stale_optical_checks(conn)
    finally:
        conn.close()
    yield


app = FastAPI(
    title="Paired FASTQ Service",
    version="1.1.0",
    description=(
        "Upload shuffled paired-end FASTQ files (R1/R2, uncompressed, "
        "max 256 MiB each). The service re-pairs them by read id, delivers "
        "both mates in R1 order and preserves the original bytes. "
        "It also exposes an optical-duplicate check over the results of a "
        "successful task."
    ),
    lifespan=lifespan,
)


class TaskOut(BaseModel):
    id: str
    status: str
    stage: str
    progress: float
    message: str
    r1_name: str
    r2_name: str
    r1_bytes: int
    r2_bytes: int
    created_at: str
    updated_at: str
    download_r1: str | None = None
    download_r2: str | None = None


def _serialize(row: sqlite3.Row) -> TaskOut:
    download_r1 = download_r2 = None
    if row["status"] == "succeeded" and row["result_r1"]:
        download_r1 = f"/tasks/{row['id']}/download/r1"
        download_r2 = f"/tasks/{row['id']}/download/r2"
    return TaskOut(
        id=row["id"],
        status=row["status"],
        stage=row["stage"],
        progress=round(row["progress"], 4),
        message=row["message"],
        r1_name=row["r1_name"],
        r2_name=row["r2_name"],
        r1_bytes=row["r1_bytes"],
        r2_bytes=row["r2_bytes"],
        created_at=row["created_at"],
        updated_at=row["updated_at"],
        download_r1=download_r1,
        download_r2=download_r2,
    )


async def _spool_bounded(upload: UploadFile, side: str, dest: Path) -> int:
    """Stream an upload to ``dest`` without holding it in memory.

    SpooledTemporaryFile rolls over to a real file after
    ``UPLOAD_SPOOL_THRESHOLD`` and the copy proceeds in fixed chunks; the
    per-file 256 MiB limit is enforced while streaming.
    """
    spool = tempfile.SpooledTemporaryFile(max_size=config.UPLOAD_SPOOL_THRESHOLD)
    total = 0
    try:
        while True:
            chunk = await upload.read(1024 * 1024)
            if not chunk:
                break
            total += len(chunk)
            if total > config.MAX_FILE_BYTES:
                raise UploadTooLarge(side, config.MAX_FILE_BYTES)
            spool.write(chunk)
        spool.rollover()
        spool.seek(0)
        with dest.open("wb") as out:
            while True:
                chunk = spool.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
    finally:
        spool.close()
    return total


@app.post("/tasks", response_model=TaskOut, status_code=status.HTTP_202_ACCEPTED)
async def create_task(
    r1: UploadFile = File(..., alias="r1", description="R1 FASTQ, uncompressed, <= 256 MiB"),
    r2: UploadFile = File(..., alias="r2", description="R2 FASTQ, uncompressed, <= 256 MiB"),
) -> TaskOut:
    task_id = uuid.uuid4().hex
    base = task_dir(task_id)
    base.mkdir(parents=True, exist_ok=False)
    r1_dest = base / "input_R1.fastq"
    r2_dest = base / "input_R2.fastq"
    try:
        size1 = await _spool_bounded(r1, "R1", r1_dest)
        size2 = await _spool_bounded(r2, "R2", r2_dest)
    except UploadTooLarge as exc:
        worker.cleanup_task_dir(task_id)
        raise HTTPException(
            status_code=getattr(
                status, "HTTP_413_CONTENT_TOO_LARGE",
                status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
            ),
            detail=str(exc),
        )
    except Exception:
        worker.cleanup_task_dir(task_id)
        raise
    finally:
        await r1.close()
        await r2.close()

    conn = db.connect()
    try:
        db.create_task(
            conn,
            task_id,
            r1_name=Path(r1.filename or "r1.fastq").name,
            r2_name=Path(r2.filename or "r2.fastq").name,
            r1_bytes=size1,
            r2_bytes=size2,
        )
        row = db.get_task(conn, task_id)
    finally:
        conn.close()

    worker.register_cancel_event(task_id)
    start_worker(task_id)
    return _serialize(row)


@app.get("/tasks/{task_id}", response_model=TaskOut)
def get_task(task_id: str) -> TaskOut:
    conn = db.connect()
    try:
        row = db.get_task(conn, task_id)
        if row is None:
            raise HTTPException(status_code=404, detail="task not found")
        return _serialize(row)
    finally:
        conn.close()


@app.post("/tasks/{task_id}/cancel", response_model=TaskOut)
def cancel_task(task_id: str) -> TaskOut:
    conn = db.connect()
    try:
        previous = db.try_set_cancelling(conn, task_id)
        if previous is None:
            raise HTTPException(status_code=404, detail="task not found")
        row = db.get_task(conn, task_id)
        if row["status"] in ("succeeded", "failed", "cancelled"):
            # The task reached a terminal state first; tell the caller plainly
            # that it cannot be cancelled.
            raise HTTPException(
                status_code=409,
                detail=f"task is already {row['status']} and cannot be cancelled",
            )
        # Flip the cooperative flag and return without waiting for the worker.
        worker.request_cancel(task_id)
        return _serialize(row)
    finally:
        conn.close()


def _result_path(task_id: str, side: int) -> Path:
    conn = db.connect()
    try:
        row = db.get_task(conn, task_id)
        if row is None:
            raise HTTPException(status_code=404, detail="task not found")
        if row["status"] != "succeeded":
            raise HTTPException(
                status_code=409,
                detail=f"task is {row['status']}; results are not published",
            )
        rel = row["result_r1"] if side == 1 else row["result_r2"]
        if not rel:
            raise HTTPException(status_code=404, detail="result file missing")
    finally:
        conn.close()

    base = task_dir(task_id).resolve()
    path = (base / rel).resolve()
    # Guard against path traversal even though the rel value is server-generated.
    if base not in path.parents or not path.is_file():
        raise HTTPException(status_code=404, detail="result file missing")
    return path


@app.get("/tasks/{task_id}/download/r1")
def download_r1(task_id: str) -> FileResponse:
    path = _result_path(task_id, 1)
    return FileResponse(
        path, media_type="text/plain", filename="R1.fastq",
        headers={"Content-Length": str(path.stat().st_size)},
    )


@app.get("/tasks/{task_id}/download/r2")
def download_r2(task_id: str) -> FileResponse:
    path = _result_path(task_id, 2)
    return FileResponse(
        path, media_type="text/plain", filename="R2.fastq",
        headers={"Content-Length": str(path.stat().st_size)},
    )


@app.get("/health")
def health() -> JSONResponse:
    return JSONResponse({"status": "ok"})


# ---------------------------------------------------------------------------
# Optical-duplicate check
# ---------------------------------------------------------------------------


class OpticalCheckRequest(BaseModel):
    task_id: str = Field(
        ...,
        description="Id of a successful re-pairing task.",
        examples=["0f3a...c9"],
    )
    distance_threshold: int = Field(
        ...,
        ge=config.OPTICAL_MIN_DISTANCE,
        le=config.OPTICAL_MAX_DISTANCE,
        description="Maximum Euclidean pixel distance (1..10) for a link.",
        examples=[2],
    )


class OpticalDuplicate(BaseModel):
    read_id: str = Field(..., description="Read id to remove.")
    keep_read_id: str = Field(..., description="Read id whose pair is kept.")


class OpticalCheckOut(BaseModel):
    id: str
    task_id: str
    state: str = Field(..., description="ready or failed.")
    distance_threshold: int
    total_pairs: int = Field(..., description="Read pairs examined.")
    duplicate_pairs: int = Field(..., description="Pairs flagged for removal.")
    duplicates: list[OpticalDuplicate] = Field(
        default_factory=list,
        description="Flagged pairs in R1 order, each with its retained pair.",
    )
    error: str = ""
    created_at: str
    updated_at: str


def _serialize_optical(row: sqlite3.Row) -> OpticalCheckOut:
    duplicates: list[OpticalDuplicate] = []
    if row["state"] == "ready" and row["report"]:
        body = json.loads(row["report"])
        duplicates = [OpticalDuplicate(**item) for item in body]
    return OpticalCheckOut(
        id=row["id"],
        task_id=row["task_id"],
        state=row["state"],
        distance_threshold=row["threshold"],
        total_pairs=row["total_pairs"],
        duplicate_pairs=row["duplicate_pairs"],
        duplicates=duplicates,
        error=row["error"],
        created_at=row["created_at"],
        updated_at=row["updated_at"],
    )


def _result_file_or_409(task_id: str, side: int) -> Path:
    """Locate a published result file, refusing anything not fully succeeded.

    This deliberately re-checks the database *and* the filesystem so a task
    whose results were removed or truncated externally fails the whole check
    instead of yielding a partial list.
    """
    conn = db.connect()
    try:
        row = db.get_task(conn, task_id)
        if row is None:
            raise HTTPException(status_code=404, detail="task not found")
        if row["status"] != "succeeded":
            raise HTTPException(
                status_code=409,
                detail=(
                    f"task is {row['status']}; only a successful re-pairing "
                    "task can be checked"
                ),
            )
        rel = row["result_r1"] if side == 1 else row["result_r2"]
        if not rel:
            raise HTTPException(
                status_code=409, detail="result files are missing for this task"
            )
    finally:
        conn.close()
    base = task_dir(task_id).resolve()
    path = (base / rel).resolve()
    if base not in path.parents or not path.is_file():
        raise HTTPException(
            status_code=409, detail=f"result R{side} file is missing"
        )
    return path


_optical_slot = threading.Lock()
"""One check at a time keeps the process inside its 256 MiB cap."""


def _run_optical_check(
    check_id: str, r1_path: Path, r2_path: Path, threshold: int
) -> sqlite3.Row:
    """Analyse, persist the outcome and return the stored report row."""
    try:
        result = optical.analyze(r1_path, r2_path, threshold)
    except (optical.OpticalCheckError, OSError) as exc:
        conn = db.connect()
        try:
            db.finish_optical_check(
                conn, check_id, state="failed", error=str(exc)
            )
        finally:
            conn.close()
        raise HTTPException(
            status_code=422,
            detail={
                "message": f"optical check failed as a whole: {exc}",
                "check_id": check_id,
            },
        ) from exc

    conn = db.connect()
    try:
        db.finish_optical_check(
            conn,
            check_id,
            state="ready",
            total_pairs=result["total_pairs"],
            duplicate_pairs=result["duplicate_pairs"],
            report=json.dumps(result["duplicates"], separators=(",", ":")),
        )
        row = db.get_optical_check(conn, check_id)
    finally:
        conn.close()
    return row


@app.post(
    "/optical-checks",
    response_model=OpticalCheckOut,
    status_code=status.HTTP_201_CREATED,
    summary="Detect optical duplicates in a successful task's results",
)
def create_optical_check(body: OpticalCheckRequest) -> OpticalCheckOut:
    """Run the whole check synchronously and return the complete report.

    Any missing file, corrupt record, invalid coordinate or quality byte
    outside ASCII 33..126 fails the entire check; no partial list is returned
    and the result files are never modified.  Only one check runs at a time so
    the process stays inside its 256 MiB cap.
    """
    r1_path = _result_file_or_409(body.task_id, 1)
    r2_path = _result_file_or_409(body.task_id, 2)

    check_id = uuid.uuid4().hex
    conn = db.connect()
    try:
        db.create_optical_check(conn, check_id, body.task_id, body.distance_threshold)
    finally:
        conn.close()

    # Only one check runs at a time.  A 50,000-pair worst case plus the API
    # process baseline stays under 256 MiB; two concurrent checks would not.
    if not _optical_slot.acquire(blocking=False):
        conn = db.connect()
        try:
            db.finish_optical_check(
                conn, check_id, state="failed",
                error="another optical check is already running",
            )
        finally:
            conn.close()
        raise HTTPException(
            status_code=429,
            detail={
                "message": "another optical check is already running; retry later",
                "check_id": check_id,
            },
        )
    try:
        row = _run_optical_check(
            check_id, r1_path, r2_path, body.distance_threshold
        )
    finally:
        _optical_slot.release()
    return _serialize_optical(row)


@app.get(
    "/optical-checks/{check_id}",
    response_model=OpticalCheckOut,
    summary="Fetch a previously produced optical-duplicate report",
)
def get_optical_check(check_id: str) -> OpticalCheckOut:
    conn = db.connect()
    try:
        row = db.get_optical_check(conn, check_id)
        if row is None:
            raise HTTPException(status_code=404, detail="optical check not found")
        return _serialize_optical(row)
    finally:
        conn.close()
