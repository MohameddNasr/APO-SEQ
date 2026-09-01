"""Reusable command planning and execution helpers."""

from __future__ import annotations

import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path

from aposeq.exceptions import ApoSeqError


@dataclass(frozen=True)
class PlannedCommand:
    """One executable command with optional redirected stdout."""

    step: str
    argv: tuple[str, ...]
    stdout_path: Path | None = None


def format_command(command: PlannedCommand) -> str:
    """Format a command for logs and manifests."""

    command_text = " ".join(shlex.quote(part) for part in command.argv)
    if command.stdout_path is not None:
        command_text = f"{command_text} > {shlex.quote(str(command.stdout_path))}"
    return command_text


def execute_command(command: PlannedCommand, context: str) -> None:
    """Execute one planned command."""

    try:
        if command.stdout_path is None:
            subprocess.run(command.argv, check=True)
        else:
            command.stdout_path.parent.mkdir(parents=True, exist_ok=True)
            with command.stdout_path.open("w", encoding="utf-8") as stdout:
                subprocess.run(command.argv, check=True, stdout=stdout)
    except FileNotFoundError as exc:
        raise ApoSeqError(f"Required executable was not found: {command.argv[0]}") from exc
    except subprocess.CalledProcessError as exc:
        raise ApoSeqError(f"Command failed for {context}: {format_command(command)}") from exc
