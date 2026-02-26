#!/usr/bin/env python3
"""
Import HOMSTRAD .ali files into MySQL (MariaDB).

- No biopython/pandas.
- Parses each family folder: <family>/<family>.ali
- Stores:
  * family (name)
  * alignment (name=family, type='msa')
  * sequences (ungapped)
  * source_entry mapping: source='pdb', source_key='<pdb>_<chain>'
  * family_member mapping: (family_id, source, source_key)
  * alignment_row aligned_seq (with gaps)

Usage example:
  python3 import_homstrad.py \
    --host mysql2-ext.bio.ifi.lmu.de \
    --db bioprakt17 --user bioprakt17 \
    --homstrad_dir /mnt/extstud/praktikum/bioprakt/Data/HOMSTRAD \
    --family 5_3_exonuclease

For all families:
  python3 import_homstrad.py ... --homstrad_dir .../HOMSTRAD
"""

import argparse
import os
from dataclasses import dataclass
from getpass import getpass
from typing import Iterator, Optional, Tuple, List

import mysql.connector


@dataclass
class AliEntry:
    pdb: str
    chain: str
    aligned_seq: str  # includes '-' gaps, no '*'


def parse_structurex_line(line: str) -> Tuple[str, str]:
    """
    Example:
      structureX:1bgx:   1 :T: 285 :T:TAQ DNA polymerase:Thermus aquaticus: 2.30: 18.9

    Split by ':' gives:
      [ 'structureX', '1bgx', '   1 ', 'T', ' 285 ', 'T', 'TAQ DNA polymerase', ... ]

    We take pdb = parts[1], chain = parts[3]
    """
    parts = line.strip().split(":")
    # Expect at least: structureX, pdb, start, chain, ...
    pdb = parts[1].strip().lower() if len(parts) > 1 else ""
    chain = parts[3].strip() if len(parts) > 3 else ""

    # Fallback: sometimes chain field is missing/empty
    if not chain:
        chain = "-"   # placeholder for "unknown chain"

    if not pdb:
        raise ValueError(f"Cannot parse PDB from structureX line: {line!r}")

    return pdb, chain

def iter_ali_entries(ali_path: str) -> Iterator[AliEntry]:
    """
    Iterates over entries in a HOMSTRAD .ali file.
    Rules:
      - ignore comment lines starting with 'C;'
      - entry starts at line '>P1;...'
      - next non-empty line should be 'structureX:...'
      - sequence lines until a line containing '*' (end marker)
    """
    with open(ali_path, "r", encoding="utf-8", errors="replace") as f:
        state = "WAIT_HEADER"
        pdb = chain = ""
        seq_chunks: List[str] = []

        for raw in f:
            line = raw.rstrip("\n")

            if not line:
                continue
            if line.startswith("C;"):
                continue

            if state == "WAIT_HEADER":
                if line.startswith(">P1;"):
                    state = "READ_META"
                    pdb = chain = ""
                    seq_chunks = []
                # else ignore
                continue

            if state == "READ_META":
                # Usually structureX line right after header
                if line.lower().startswith("structurex:"):
                    pdb, chain = parse_structurex_line(line)
                    state = "READ_SEQ"
                else:
                    # Sometimes there might be other metadata lines; keep waiting for structureX
                    continue
                continue

            if state == "READ_SEQ":
                # alignment lines may contain spaces; remove spaces
                clean = line.replace(" ", "")
                if "*" in clean:
                    before_star = clean.split("*", 1)[0]
                    if before_star:
                        seq_chunks.append(before_star)
                    aligned_seq = "".join(seq_chunks)
                    if not pdb or not aligned_seq:
                        raise ValueError(f"Incomplete entry in {ali_path}: pdb={pdb}, chain={chain}")
                    yield AliEntry(pdb=pdb, chain=chain, aligned_seq=aligned_seq)
                    state = "WAIT_HEADER"
                else:
                    seq_chunks.append(clean)
                continue


# --------------------------- DB helpers ---------------------------


def get_or_create_family_id(cur, family_name: str) -> int:
    cur.execute("INSERT IGNORE INTO family (name) VALUES (%s)", (family_name,))
    cur.execute("SELECT family_id FROM family WHERE name=%s", (family_name,))
    row = cur.fetchone()
    if not row:
        raise RuntimeError(f"Could not fetch family_id for {family_name}")
    return int(row[0])


def get_or_create_alignment_id(cur, family_name: str) -> int:
    cur.execute(
        "INSERT IGNORE INTO alignment (name, type) VALUES (%s, %s)",
        (family_name, "msa"),
    )
    cur.execute(
        "SELECT alignment_id FROM alignment WHERE name=%s AND type=%s",
        (family_name, "msa"),
    )
    row = cur.fetchone()
    if not row:
        raise RuntimeError(f"Could not fetch alignment_id for {family_name}")
    return int(row[0])


def get_sequence_id_if_exists(cur, source: str, source_key: str) -> Optional[int]:
    cur.execute(
        "SELECT sequence_id FROM source_entry WHERE source=%s AND source_key=%s",
        (source, source_key),
    )
    row = cur.fetchone()
    return int(row[0]) if row else None


def create_sequence_and_source_entry(cur, source: str, source_key: str, ungapped_seq: str) -> int:
    # Insert sequence
    cur.execute(
        "INSERT INTO sequence (seq, length, organism_id) VALUES (%s, %s, NULL)",
        (ungapped_seq, len(ungapped_seq)),
    )
    seq_id = cur.lastrowid
    if not seq_id:
        # Fallback: fetch last id
        cur.execute("SELECT LAST_INSERT_ID()")
        seq_id = int(cur.fetchone()[0])

    # Link source entry
    cur.execute(
        "INSERT IGNORE INTO source_entry (source, source_key, sequence_id) VALUES (%s, %s, %s)",
        (source, source_key, seq_id),
    )
    return int(seq_id)


def ensure_family_member(cur, family_id: int, source: str, source_key: str) -> None:
    cur.execute(
        "INSERT IGNORE INTO family_member (family_id, source, source_key) VALUES (%s, %s, %s)",
        (family_id, source, source_key),
    )


def insert_alignment_row(cur, alignment_id: int, sequence_id: int, aligned_seq: str) -> None:
    # Requires a UNIQUE key on (alignment_id, sequence_id) to avoid duplicates.
    cur.execute(
        "INSERT IGNORE INTO alignment_row (alignment_id, sequence_id, aligned_seq) VALUES (%s, %s, %s)",
        (alignment_id, sequence_id, aligned_seq),
    )


# --------------------------- main import ---------------------------

def import_family(conn, homstrad_dir: str, family_name: str) -> None:
    ali_path = os.path.join(homstrad_dir, family_name, f"{family_name}.ali")
    if not os.path.isfile(ali_path):
        print(f"[SKIP] No .ali file: {ali_path}")
        return

    with conn.cursor(buffered=True) as cur:
        family_id = get_or_create_family_id(cur, family_name)
        alignment_id = get_or_create_alignment_id(cur, family_name)

        n_entries = 0
        for entry in iter_ali_entries(ali_path):
            source = "pdb"
            source_key = f"{entry.pdb}_{entry.chain}"

            # sequence (ungapped)
            ungapped = entry.aligned_seq.replace("-", "")
            seq_id = get_sequence_id_if_exists(cur, source, source_key)
            if seq_id is None:
                seq_id = create_sequence_and_source_entry(cur, source, source_key, ungapped)

            # membership + aligned row
            ensure_family_member(cur, family_id, source, source_key)
            insert_alignment_row(cur, alignment_id, seq_id, entry.aligned_seq)

            n_entries += 1

        conn.commit()
        print(f"[OK] Imported family '{family_name}' with {n_entries} sequences.")


def iter_family_names(homstrad_dir: str) -> List[str]:
    fams = []
    for name in os.listdir(homstrad_dir):
        p = os.path.join(homstrad_dir, name)
        if os.path.isdir(p):
            fams.append(name)
    fams.sort()
    return fams


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", required=True)
    ap.add_argument("--db", required=True)
    ap.add_argument("--user", required=True)
    ap.add_argument("--password", default=None, help="If omitted, read from env MYSQL_PW or prompt.")
    ap.add_argument("--homstrad_dir", required=True)
    ap.add_argument("--family", default=None, help="Import only one family folder (e.g. 5_3_exonuclease).")
    args = ap.parse_args()

    pw = args.password or os.environ.get("MYSQL_PW") or getpass("MySQL password: ")

    conn = mysql.connector.connect(
        host=args.host,
        user=args.user,
        password=pw,
        database=args.db,
    )

    try:
        if args.family:
            import_family(conn, args.homstrad_dir, args.family)
        else:
            for fam in iter_family_names(args.homstrad_dir):
                import_family(conn, args.homstrad_dir, fam)
    finally:
        conn.close()


if __name__ == "__main__":
    main()

