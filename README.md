# Protein Secondary Structure Prediction & Sequence Alignment

Group project from the Bioinformatics programming lab (Programmierpraktikum) at LMU Munich / TUM, winter semester 2025/26.

The main project (`alignment-and-gor/`) implements two classic bioinformatics methods from scratch in Java and evaluates them on standard benchmarks:

1. **Pairwise sequence alignment** with affine gap costs (Gotoh), benchmarked on HOMSTRAD structural alignments
2. **Protein secondary structure prediction** with the GOR family of methods (GOR I, III, IV and the MSA-based GOR V), evaluated on CB513 and a CASP test set

## Results: secondary structure prediction

Mean Q3 accuracy (share of residues correctly predicted as helix, strand or coil), 5-fold cross-validation on CB513:

| Method | Q3 (%) |
|---|---|
| **GOR V (multiple sequence alignments) + GOR I** | **67.1** |
| GOR V (multiple sequence alignments) + GOR III | 66.2 |
| GOR I | 61.4 |
| GOR III | 58.9 |

Using multiple sequence alignments (GOR V) improves Q3 by about 6 percentage points over single-sequence GOR I. For reference, PSIPRED reaches 79.9 % on the CASP test set, where our GOR I reaches 60.1 %.

## Project structure

```
alignment-and-gor/
├── Alignment/            # Gotoh, global, local and free-shift alignment (Java)
│   └── Silja/            # HOMSTRAD parser, benchmark and gap-parameter tuning (Python)
└── GOR/                  # GOR I/III/IV/V training, prediction and evaluation (Java)
tools/                    # cross-validation, hyperparameter tuning, CASP evaluation (Python)
sequence-analysis/        # sequence tools: ORF finder, translation, PROSITE scan, CRISPR search, k-mers, web CGI
databases-and-structures/ # SQL databases, PDB structures, contact maps, HOMSTRAD import, secondary structure
```

## How to run

GOR prediction (from `alignment-and-gor/GOR`):

```bash
javac -d out $(find src -name "*.java")
java -cp out de.lmu.ifi.bio.gor.cli.Train   # see the class for arguments
java -cp out de.lmu.ifi.bio.gor.cli.Predict
```

Alignment benchmark:

```bash
python alignment-and-gor/Alignment/Silja/run_benchmark.py --quick
```

Scripts that use a database read their credentials from environment variables (`DB_HOST`, `DB_USER`, `DB_PASS`).

## Team

Silja Ohlsen and Zeno Fiocco (sequence alignment) · Lara Groß and Anna Niederlechner (GOR)
Supervision: Prof. Dr. Volker Heun, Felix Offensperger, Armin Hadziahmetovic

The web interface for this project lives in [protein-alignment-server](https://github.com/SiljaOh/protein-alignment-server).
