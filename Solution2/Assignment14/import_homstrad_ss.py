#!/usr/bin/env python3

import os
import sys
import argparse
from pathlib import Path
from typing import List, Dict, Optional, Tuple

# Pure-Python MySQL client. Install via: pip3 install --user pymysql
import pymysql

SS_HEADER = "secondary structure and phi angle"


def parse_tem_secondary_structure(tem_path: Path) -> List[str]:
    ss_list: List[str] = []

    with tem_path.open("r", encoding="utf-8", errors="replace") as f:
        lines = [line.rstrip("\n") for line in f]

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # Look for ">P1;..."
        if line.startswith(">P1;"):
            # Next line may be "sequence" OR "secondary structure and phi angle" etc.
            if i + 1 < len(lines) and lines[i + 1].strip() == SS_HEADER:
                # Collect subsequent lines until we hit a line containing '*'
                j = i + 2
                parts: List[str] = []
                while j < len(lines):
                    s = lines[j].strip()
                    if not s:
                        j += 1
                        continue
                    if "*" in s:
                        # include everything before '*'
                        before_star = s.split("*", 1)[0]
                        parts.append(before_star)
                        break
                    parts.append(s)
                    j += 1

                aligned_ss = "".join(parts).strip()
                if aligned_ss:
                    ss_list.append(aligned_ss)

                i = j + 1
                continue

        i += 1

    return ss_list


def db_connect(host: str, port: int, user: str, password: str, dbname: str):
    return pymysql.connect(
        host=host,
        port=port,
        user=user,
        password=password,
        database=dbname,
        charset="utf8mb4",
        autocommit=False,
        cursorclass=pymysql.cursors.DictCursor,
    )


def ensure_secondary_structure_pk(conn):
    with conn.cursor() as cur:
        cur.execute("SHOW INDEX FROM secondary_structure;")
        idx = cur.fetchall()

    # Check if there's already a unique index covering (alignment_id, sequence_id)
    # We accept PK or UNIQUE.
    unique_keys: Dict[str, List[str]] = {}
    for row in idx:
        if row.get("Non_unique") == 0:
            key_name = row.get("Key_name")
            col_name = row.get("Column_name")
            unique_keys.setdefault(key_name, []).append(col_name)

    for cols in unique_keys.values():
        if cols[:2] == ["alignment_id", "sequence_id"]:
            return  # good enough

    try:
        with conn.cursor() as cur:
            cur.execute(
                "ALTER TABLE secondary_structure "
                "ADD PRIMARY KEY (alignment_id, sequence_id);"
            )
        conn.commit()
        print("[OK] Added PRIMARY KEY (alignment_id, sequence_id) to secondary_structure.")
    except Exception as e:
        conn.rollback()
        print(
            "[WARN] Could not add PRIMARY KEY (alignment_id, sequence_id). "
            "UPSERT may fail if duplicates exist.\n"
            f"       Error: {e}"
        )


def get_alignment_id(conn, alignment_name: str) -> Optional[int]:
    with conn.cursor() as cur:
        cur.execute("SELECT alignment_id FROM alignment WHERE name=%s", (alignment_name,))
        row = cur.fetchone()
    return None if row is None else int(row["alignment_id"])


def get_alignment_sequence_ids(conn, alignment_id: int) -> List[int]:
    # Use row_id ordering (as you have row_id in alignment_row)
    with conn.cursor() as cur:
        cur.execute(
            "SELECT sequence_id FROM alignment_row WHERE alignment_id=%s ORDER BY row_id",
            (alignment_id,),
        )
        rows = cur.fetchall()
    return [int(r["sequence_id"]) for r in rows]


def upsert_secondary_structure(conn, alignment_id: int, sequence_id: int, aligned_ss: str):
    sql = (
        "INSERT INTO secondary_structure (alignment_id, sequence_id, aligned_ss) "
        "VALUES (%s, %s, %s) "
        "ON DUPLICATE KEY UPDATE aligned_ss=VALUES(aligned_ss)"
    )
    with conn.cursor() as cur:
        cur.execute(sql, (alignment_id, sequence_id, aligned_ss))


def process_family_folder(conn, folder: Path, dry_run: bool = False) -> Tuple[int, int]:
    folder_name = folder.name
    tem_path = folder / f"{folder_name}.tem"
    if not tem_path.exists():
        return (0, 0)

    alignment_id = get_alignment_id(conn, folder_name)
    if alignment_id is None:
        print(f"[WARN] No alignment found with name='{folder_name}' (skipping)")
        return (0, 1)

    ss_list = parse_tem_secondary_structure(tem_path)
    if not ss_list:
        print(f"[WARN] No SS blocks found in {tem_path} (skipping)")
        return (0, 1)

    seq_ids = get_alignment_sequence_ids(conn, alignment_id)
    if not seq_ids:
        print(f"[WARN] No alignment_row entries for alignment_id={alignment_id} (skipping)")
        return (0, 1)

    if len(ss_list) != len(seq_ids):
        print(
            f"[WARN] Count mismatch in '{folder_name}': "
            f"{len(ss_list)} SS blocks in .tem but {len(seq_ids)} alignment_row entries. "
            "Will import up to min(counts)."
        )

    n = min(len(ss_list), len(seq_ids))
    changed = 0

    if dry_run:
        print(f"[DRY] Would import {n} secondary_structure rows for '{folder_name}' (alignment_id={alignment_id})")
        return (0, 0)

    for i in range(n):
        upsert_secondary_structure(conn, alignment_id, seq_ids[i], ss_list[i])
        changed += 1

    return (changed, 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="Path to HOMSTRAD_DB/2019_Feb folder")
    ap.add_argument("--host", default=os.environ.get("DB_HOST", "127.0.0.1"))
    ap.add_argument("--port", type=int, default=int(os.environ.get("DB_PORT", "3306")))
    ap.add_argument("--user", default=os.environ.get("DB_USER", "bioprakt17"))
    ap.add_argument("--password", default=os.environ.get("DB_PASS", ""))
    ap.add_argument("--db", dest="dbname", default=os.environ.get("DB_NAME", "bioprakt17"))
    ap.add_argument("--dry-run", action="store_true", help="Parse and match but do not write to DB")
    args = ap.parse_args()

    root = Path(args.root)
    if not root.exists() or not root.is_dir():
        print(f"ERROR: root '{root}' does not exist or is not a directory", file=sys.stderr)
        sys.exit(1)

    if not args.password:
        print(
            "ERROR: No DB password provided. Set DB_PASS env var or pass --password.",
            file=sys.stderr,
        )
        sys.exit(1)

    conn = db_connect(args.host, args.port, args.user, args.password, args.dbname)
    try:
        ensure_secondary_structure_pk(conn)

        total_changed = 0
        total_warn = 0
        folders = [p for p in root.iterdir() if p.is_dir()]

        print(f"[OK] Found {len(folders)} family folders under {root}")

        for folder in sorted(folders):
            changed, warn = process_family_folder(conn, folder, dry_run=args.dry_run)
            total_changed += changed
            total_warn += warn

            if not args.dry_run:
                conn.commit()

        print(f"\n[DONE] Upserted {total_changed} secondary_structure rows. Warnings: {total_warn}")

    except KeyboardInterrupt:
        print("\n[ABORT] KeyboardInterrupt - rolling back current transaction.")
        conn.rollback()
        sys.exit(2)
    except Exception as e:
        print(f"\n[ERROR] {e}")
        conn.rollback()
        sys.exit(3)
    finally:
        conn.close()


if __name__ == "__main__":
    main()
