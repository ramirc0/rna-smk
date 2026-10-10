# Config, sample sheet, and shared helpers.

from pathlib import Path

import polars as pl
from snakemake.exceptions import WorkflowError


configfile: "config/config.yaml"


RESULTS = config["outdir"]
RUN_ID = config.get("run_id", "default")
OUTDIR = f"{RESULTS}/{RUN_ID}"          # per-run
REFDIR = f"{RESULTS}/refs"              # run-independent
LOGDIR = f"logs/{RUN_ID}"
BENCHDIR = f"benchmarks/{RUN_ID}"
REFERENCE_FASTA = config["references"]["fasta"]
REFERENCE_GTF = config["references"]["gtf"]

STAR_INDEX = f"{REFDIR}/index"
RSEM_INDEX = f"{REFDIR}/rsem_index/rsem"
RSEM_INDEX_FILES = multiext(
    RSEM_INDEX,
    ".chrlist",
    ".grp",
    ".idx.fa",
    ".n2g.idx.fa",
    ".seq",
    ".ti",
    ".transcripts.fa",
)

# Plain patterns, not multiext(): expand() rejects multiext-flagged paths.
SIGNAL_TRACKS = expand(
    f"{OUTDIR}/{{sample}}/signals/{{sample}}_genome_{{suffix}}.bw",
    suffix=(
        ["minusAll", "minusUniq", "plusAll", "plusUniq"]
        if config["sequencing"]["is_stranded"]
        else ["all", "uniq"]
    ),
    allow_missing=True,
)

# The workspace path is relative to the rule file. Every rule file sits in rules/.
SOFTWARE_ENV = pixi(workspace="../envs", env=config["pixi_env"], locked=True)


_SHEET = Path(config["samples"])
if not _SHEET.exists():
    raise WorkflowError(
        f"Sample sheet not found: {_SHEET}\n"
        "Columns: sample_id, r1, r2; see config/samples.example.tsv."
    )

# All-string read so empty cells are "" (absent), not a type error.
_manifest = pl.read_csv(_SHEET, separator="\t", infer_schema_length=0, comment_prefix="#")


def _column(name):
    """Map each sample to its value in sheet column `name`.

    Parameters
    ----------
    name : str
        Sample sheet column.

    Returns
    -------
    dict of str to str or None
        Value per `sample_id`, with empty cells as None. Empty when the sheet
        has no such column.
    """
    if name not in _manifest.columns:
        return {}
    return {
        row["sample_id"]: (row[name] or None)
        for row in _manifest.iter_rows(named=True)
    }


SAMPLES = _manifest.get_column("sample_id").to_list()
R1_OF = _column("r1")
R2_OF = _column("r2")

if not SAMPLES:
    raise WorkflowError(f"Sample sheet {_SHEET} has no rows.")


wildcard_constraints:
    sample=r"[A-Za-z0-9][A-Za-z0-9_.-]*",
