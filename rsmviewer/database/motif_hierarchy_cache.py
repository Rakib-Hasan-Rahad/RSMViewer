"""Text-keyed hierarchy and saved-query cache for RSMViewer."""

from __future__ import annotations

import json
import sqlite3
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

_MAX_AGE_SECONDS = 30 * 24 * 60 * 60
_DEFAULT_JACCARD_THRESHOLD = 0.60


def residue_key(residues: List[Tuple[str, int]]) -> str:
    """Return a stable order-independent residue key."""
    return ";".join(f"{chain}:{number}" for chain, number in sorted(set(residues)))


def _parse_residue_key_set(value: str):
    result = set()
    for token in value.split(";"):
        chain, _, number = token.partition(":")
        try:
            result.add((chain, int(number)))
        except ValueError:
            continue
    return result


def _jaccard(first, second) -> float:
    if not first or not second:
        return 0.0
    union = first | second
    return len(first & second) / len(union) if union else 0.0


class MotifHierarchyCache:
    """Per-structure residue hierarchy keyed only by canonical source names."""

    def __init__(self, db_path: str):
        self.db_path = db_path
        Path(db_path).parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(db_path)
        self._init_schema()
        self.purge_expired()

    def _init_schema(self) -> None:
        existing = {
            row[1]
            for row in self._conn.execute("PRAGMA table_info(motif_hierarchy)").fetchall()
        }
        if existing and ("source_id" in existing or "source_key" not in existing):
            self._migrate_legacy_hierarchy()
        self._conn.execute(
            """
            CREATE TABLE IF NOT EXISTS motif_hierarchy (
                pdb_id TEXT NOT NULL,
                residue_key TEXT NOT NULL,
                custom_id TEXT NOT NULL,
                source_key TEXT NOT NULL,
                source_label TEXT NOT NULL,
                motif_label TEXT NOT NULL,
                hierarchy_level INTEGER NOT NULL DEFAULT 1,
                rank INTEGER NOT NULL,
                cached_at REAL NOT NULL,
                PRIMARY KEY (pdb_id, residue_key, source_key, hierarchy_level)
            )
            """
        )
        self._conn.execute(
            "CREATE TABLE IF NOT EXISTS motif_id_sequence (pdb_id TEXT PRIMARY KEY, next_seq INTEGER NOT NULL)"
        )
        self._conn.execute(
            """
            CREATE TABLE IF NOT EXISTS query_alias (
                alias TEXT PRIMARY KEY,
                pdb_id TEXT NOT NULL,
                motif_filter TEXT,
                db_expression TEXT,
                residue_keys TEXT NOT NULL,
                labels_json TEXT,
                source_command TEXT,
                created_at REAL NOT NULL
            )
            """
        )
        self._conn.commit()

    def _migrate_legacy_hierarchy(self) -> None:
        """Rebuild old integer-keyed cache tables using canonical source keys."""
        self._conn.execute("ALTER TABLE motif_hierarchy RENAME TO motif_hierarchy_legacy")
        self._conn.execute(
            """
            CREATE TABLE motif_hierarchy (
                pdb_id TEXT NOT NULL,
                residue_key TEXT NOT NULL,
                custom_id TEXT NOT NULL,
                source_key TEXT NOT NULL,
                source_label TEXT NOT NULL,
                motif_label TEXT NOT NULL,
                hierarchy_level INTEGER NOT NULL DEFAULT 1,
                rank INTEGER NOT NULL,
                cached_at REAL NOT NULL,
                PRIMARY KEY (pdb_id, residue_key, source_key, hierarchy_level)
            )
            """
        )
        columns = {row[1] for row in self._conn.execute("PRAGMA table_info(motif_hierarchy_legacy)").fetchall()}
        if "source_key" in columns:
            source_expr = "source_key"
        else:
            source_expr = "CASE source_id WHEN 3 THEN 'RNA3DMotifAtlas' WHEN 4 THEN 'Rfam' WHEN 5 THEN 'FR3D' WHEN 7 THEN 'RNAMotifScanX' ELSE 'legacy:' || source_id END"
        self._conn.execute(
            f"""
            INSERT OR IGNORE INTO motif_hierarchy
                (pdb_id, residue_key, custom_id, source_key, source_label, motif_label, hierarchy_level, rank, cached_at)
            SELECT pdb_id, residue_key, custom_id, {source_expr}, source_label, motif_label,
                   hierarchy_level, rank, cached_at
            FROM motif_hierarchy_legacy
            """
        )
        self._conn.execute("DROP TABLE motif_hierarchy_legacy")
        self._conn.commit()

    def close(self) -> None:
        """Close the cache connection when the owner is finished."""
        if getattr(self, "_conn", None) is not None:
            self._conn.close()
            self._conn = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def purge_expired(self, max_age_seconds: int = _MAX_AGE_SECONDS) -> int:
        cursor = self._conn.execute(
            "DELETE FROM motif_hierarchy WHERE cached_at < ?",
            (time.time() - max_age_seconds,),
        )
        self._conn.commit()
        return cursor.rowcount

    def reset_pdb(self, pdb_id: str) -> None:
        self._conn.execute("DELETE FROM motif_hierarchy WHERE pdb_id = ?", (pdb_id,))
        self._conn.execute("DELETE FROM motif_id_sequence WHERE pdb_id = ?", (pdb_id,))
        self._conn.commit()

    def clear_all_hierarchy_data(self) -> None:
        self._conn.execute("DELETE FROM motif_hierarchy")
        self._conn.execute("DELETE FROM motif_id_sequence")
        self._conn.commit()

    def _next_custom_id(self, pdb_id: str) -> str:
        row = self._conn.execute(
            "SELECT next_seq FROM motif_id_sequence WHERE pdb_id = ?", (pdb_id,)
        ).fetchone()
        sequence = row[0] if row else 1
        self._conn.execute(
            "INSERT INTO motif_id_sequence VALUES (?, ?) ON CONFLICT(pdb_id) DO UPDATE SET next_seq = ?",
            (pdb_id, sequence + 1, sequence + 1),
        )
        return f"{pdb_id}{sequence:04d}"

    def get_or_assign_custom_id(
        self, pdb_id: str, r_key: str, jaccard_threshold: float = _DEFAULT_JACCARD_THRESHOLD
    ) -> str:
        row = self._conn.execute(
            "SELECT custom_id FROM motif_hierarchy WHERE pdb_id = ? AND residue_key = ? LIMIT 1",
            (pdb_id, r_key),
        ).fetchone()
        if row:
            return row[0]
        incoming = _parse_residue_key_set(r_key)
        for existing_key, existing_id in self._conn.execute(
            "SELECT DISTINCT residue_key, custom_id FROM motif_hierarchy WHERE pdb_id = ?",
            (pdb_id,),
        ):
            existing = _parse_residue_key_set(existing_key)
            if existing and (
                incoming.issubset(existing)
                or existing.issubset(incoming)
                or _jaccard(incoming, existing) >= jaccard_threshold
            ):
                return existing_id
        return self._next_custom_id(pdb_id)

    def record(
        self,
        pdb_id: str,
        r_key: str,
        source_key: str,
        source_label: str,
        motif_label: str,
        rank: int,
        hierarchy_level: int = 1,
        jaccard_threshold: float = _DEFAULT_JACCARD_THRESHOLD,
    ) -> str:
        source_key = str(source_key).strip()
        if not source_key:
            raise ValueError("source_key cannot be empty")
        custom_id = self.get_or_assign_custom_id(pdb_id, r_key, jaccard_threshold)
        self._conn.execute(
            """
            INSERT INTO motif_hierarchy
                (pdb_id, residue_key, custom_id, source_key, source_label, motif_label, hierarchy_level, rank, cached_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(pdb_id, residue_key, source_key, hierarchy_level) DO UPDATE SET
                custom_id=excluded.custom_id, source_label=excluded.source_label,
                motif_label=excluded.motif_label, rank=excluded.rank, cached_at=excluded.cached_at
            """,
            (pdb_id, r_key, custom_id, source_key, source_label, motif_label, hierarchy_level, rank, time.time()),
        )
        self._conn.commit()
        return custom_id

    def get_hierarchy_for_pdb(self, pdb_id: str) -> Dict[str, List[Dict]]:
        rows = self._conn.execute(
            "SELECT residue_key, custom_id, source_key, source_label, motif_label, hierarchy_level, rank FROM motif_hierarchy WHERE pdb_id = ? ORDER BY residue_key, source_key, hierarchy_level",
            (pdb_id,),
        ).fetchall()
        result: Dict[str, List[Dict]] = {}
        for r_key, custom_id, source_key, source_label, motif_label, level, rank in rows:
            result.setdefault(r_key, []).append({
                "custom_id": custom_id,
                "source_key": source_key,
                "source_label": source_label,
                "motif_label": motif_label,
                "hierarchy_level": level,
                "rank": rank,
            })
        return result

    def alias_exists(self, alias: str) -> bool:
        return self._conn.execute("SELECT 1 FROM query_alias WHERE alias = ?", (alias,)).fetchone() is not None

    def create_alias(self, alias: str, pdb_id: str, residue_keys: List[str], motif_filter: str = "", db_expression: str = "", labels_snapshot: Optional[list] = None, source_command: str = "") -> None:
        self._conn.execute(
            "INSERT INTO query_alias VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            (alias, pdb_id, motif_filter, db_expression, json.dumps(residue_keys), json.dumps(labels_snapshot or []), source_command, time.time()),
        )
        self._conn.commit()

    def get_alias(self, alias: str) -> Optional[Dict]:
        row = self._conn.execute(
            "SELECT alias, pdb_id, motif_filter, db_expression, residue_keys, labels_json, source_command, created_at FROM query_alias WHERE alias = ?",
            (alias,),
        ).fetchone()
        if not row:
            return None
        name, pdb_id, motif_filter, expression, keys, labels, command, created = row
        return {"alias": name, "pdb_id": pdb_id, "motif_filter": motif_filter, "db_expression": expression, "residue_keys": json.loads(keys), "labels_snapshot": json.loads(labels) if labels else [], "source_command": command, "created_at": created}

    def list_aliases(self) -> List[Dict]:
        rows = self._conn.execute("SELECT alias, pdb_id, motif_filter, db_expression, source_command, created_at FROM query_alias ORDER BY created_at DESC").fetchall()
        return [{"alias": a, "pdb_id": p, "motif_filter": m, "db_expression": e, "source_command": c, "created_at": t} for a, p, m, e, c, t in rows]

    def delete_alias(self, alias: str) -> bool:
        cursor = self._conn.execute("DELETE FROM query_alias WHERE alias = ?", (alias,))
        self._conn.commit()
        return cursor.rowcount > 0

    def reset_all_aliases(self) -> None:
        self._conn.execute("DELETE FROM query_alias")
        self._conn.commit()


_instance: Optional[MotifHierarchyCache] = None


def get_hierarchy_cache(db_path: Optional[str] = None) -> MotifHierarchyCache:
    global _instance
    if _instance is None or getattr(_instance, "_conn", None) is None:
        db_path = db_path or str(Path(__file__).parent / "motif_hierarchy_cache.sqlite3")
        _instance = MotifHierarchyCache(db_path)
    return _instance


def close_hierarchy_cache() -> None:
    """Close and drop the singleton so the next call reopens a fresh connection."""
    global _instance
    if _instance is not None:
        try:
            _instance.close()
        except Exception:
            pass
        _instance = None
