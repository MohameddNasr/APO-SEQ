"""Pipeline state tracking for APO-SEQ resumability."""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class StepStatus:
    """State for one pipeline step."""

    name: str
    status: str
    started_at: str | None = None
    completed_at: str | None = None
    error: str | None = None


def state_path(output_directory: Path, filename: str = "pipeline_state.json") -> Path:
    return output_directory / filename


def load_state(path: Path) -> dict[str, Any]:
    """Load existing state or create an empty state document."""

    if not path.exists():
        return {"steps": {}, "created_at": utc_now(), "updated_at": utc_now()}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def save_state(state: dict[str, Any], path: Path) -> None:
    """Write state to disk."""

    state["updated_at"] = utc_now()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(state, handle, indent=2, sort_keys=True)
        handle.write("\n")


def mark_started(state: dict[str, Any], path: Path, step: str) -> None:
    state.setdefault("steps", {})[step] = {
        "status": "running",
        "started_at": utc_now(),
        "completed_at": None,
        "error": None,
    }
    save_state(state, path)


def mark_completed(state: dict[str, Any], path: Path, step: str) -> None:
    existing = state.setdefault("steps", {}).get(step, {})
    existing["status"] = "completed"
    existing["completed_at"] = utc_now()
    existing["error"] = None
    state["steps"][step] = existing
    save_state(state, path)


def mark_skipped(state: dict[str, Any], path: Path, step: str) -> None:
    state.setdefault("steps", {})[step] = {
        "status": "skipped",
        "started_at": None,
        "completed_at": utc_now(),
        "error": None,
    }
    save_state(state, path)


def mark_failed(state: dict[str, Any], path: Path, step: str, error: Exception) -> None:
    existing = state.setdefault("steps", {}).get(step, {})
    existing["status"] = "failed"
    existing["completed_at"] = utc_now()
    existing["error"] = str(error)
    state["steps"][step] = existing
    save_state(state, path)


def is_completed(state: dict[str, Any], step: str) -> bool:
    return state.get("steps", {}).get(step, {}).get("status") == "completed"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()
