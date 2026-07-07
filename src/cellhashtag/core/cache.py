"""SQLite-backed evidence cache for LATS search results."""

import sqlite3
import json
import time
import os
from pathlib import Path
from typing import Optional


class EvidenceCache:
    CACHE_DIR = Path.home() / ".cellhashtag" / "cache"
    DB_PATH = CACHE_DIR / "evidence.db"

    def __init__(self):
        self._init_db()

    def _init_db(self):
        self.CACHE_DIR.mkdir(parents=True, exist_ok=True)
        with sqlite3.connect(self.DB_PATH) as conn:
            conn.execute("""
                CREATE TABLE IF NOT EXISTS evidence (
                    key TEXT PRIMARY KEY,
                    value TEXT NOT NULL,
                    created_at REAL NOT NULL,
                    ttl REAL NOT NULL
                )
            """)
            conn.commit()

    def get(self, key: str) -> Optional[dict]:
        with sqlite3.connect(self.DB_PATH) as conn:
            row = conn.execute(
                "SELECT value, created_at, ttl FROM evidence WHERE key = ?", (key,)
            ).fetchone()
            if row is None:
                return None
            value, created, ttl = row
            if time.time() - created > ttl:
                conn.execute("DELETE FROM evidence WHERE key = ?", (key,))
                conn.commit()
                return None
            return json.loads(value)

    def put(self, key: str, data: dict, ttl: float = 3600.0):
        with sqlite3.connect(self.DB_PATH) as conn:
            conn.execute(
                "INSERT OR REPLACE INTO evidence (key, value, created_at, ttl) VALUES (?, ?, ?, ?)",
                (key, json.dumps(data), time.time(), ttl),
            )
            conn.commit()

    def clear_expired(self):
        now = time.time()
        with sqlite3.connect(self.DB_PATH) as conn:
            conn.execute("DELETE FROM evidence WHERE created_at + ttl < ?", (now,))
            conn.commit()
