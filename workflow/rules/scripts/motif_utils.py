"""Shared, deterministic helpers for site-level GRIP-seq motif analysis."""
import gzip
import hashlib
import heapq
import math
import re
from pathlib import Path

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT",
    "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
IUPAC_BY_BASES = {frozenset(v): k for k, v in IUPAC.items() if k != "U"}


def revcomp(sequence):
    return sequence.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def text_open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def stable_hash(value, seed=17):
    return int.from_bytes(hashlib.blake2b(f"{seed}|{value}".encode(), digest_size=8).digest(), "big")


def locus_key(chrom, position, strand):
    return f"{chrom}:{position}:{strand}"


def locus_partition(locus, fraction=0.2, seed=17, block_size=10000, guard=100):
    """Genomic blocks shared by both strands, with unused boundary guards.

    Entire discovery and foreground/background validation windows remain within
    a block, so neighboring overlapping site sequences cannot leak across splits.
    """
    if block_size <= 2 * guard or guard < 0 or not 0 <= fraction < 1:
        raise ValueError("holdout block_size > 2*guard >= 0 and fraction in [0,1) required")
    chrom, position, strand = locus.rsplit(":", 2)
    position = int(position)
    within_block = position % block_size
    if within_block < guard or within_block >= block_size - guard:
        return "boundary_guard"
    key = f"{chrom}:block:{position // block_size}"
    return "holdout" if stable_hash(key, seed) / 2**64 < fraction else "discovery"


def is_holdout(locus, fraction=0.2, seed=17, block_size=10000, guard=100):
    return locus_partition(locus, fraction, seed, block_size, guard) == "holdout"


def fasta_records(path):
    name, seq = None, []
    with text_open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper()
                name, seq = line[1:], []
            elif line:
                if name is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                seq.append(line)
    if name is not None:
        yield name, "".join(seq).upper()


def record_locus(name):
    match = re.search(r"(?:^|\|)locus=([^|\s]+)", name)
    if match is None:
        raise ValueError("FASTA headers must include |locus=chrom:crosslink_pos0:RNAstrand")
    return match.group(1)


def select_discovery_sequences(paths, max_sequences=2000, holdout_fraction=0.2, seed=17, holdout_block_size=10000, holdout_guard=100):
    """Deduplicate loci, reject conflicting sequence centers, and hash-sample.

    Discovery/holdout is assigned to genomic blocks shared across RNA strands.
    A first streaming pass detects loci with incompatible sequences or centers;
    a second pass selects a bounded set independently of pool file order.
    Equal sequences at different loci remain distinct genomic observations.
    """
    if max_sequences < 1 or not 0 <= holdout_fraction < 1:
        raise ValueError("max_sequences >= 1 and 0 <= holdout_fraction < 1 required")
    seen, ambiguous = {}, set()
    stats = {"input_records": 0, "duplicate_loci": 0, "holdout_loci": 0,
             "ambiguous_sequences": 0, "unique_loci": 0, "holdout_boundary_guard_loci": 0}
    for path in paths:
        for name, seq in fasta_records(path):
            stats["input_records"] += 1
            key = record_locus(name)
            center = re.search(r"(?:^|\|)center=([^|\s]+)", name)
            fingerprint = stable_hash(seq + "|" + (center.group(1) if center else ""), 0)
            if key in seen:
                stats["duplicate_loci"] += 1
                if seen[key] != fingerprint:
                    ambiguous.add(key)
            else:
                seen[key] = fingerprint
    stats["unique_loci"] = len(seen)
    stats["ambiguous_duplicate_loci"] = len(ambiguous)
    processed, heap = set(), []
    for path in paths:
        for name, seq in fasta_records(path):
            key = record_locus(name)
            if key in processed or key in ambiguous:
                continue
            processed.add(key)
            if not seq or set(seq) - set("ACGT"):
                stats["ambiguous_sequences"] += 1
                continue
            partition = locus_partition(key, holdout_fraction, seed, holdout_block_size, holdout_guard)
            if partition == "boundary_guard":
                stats["holdout_boundary_guard_loci"] += 1
                continue
            if partition == "holdout":
                stats["holdout_loci"] += 1
                continue
            center = re.search(r"(?:^|\|)center=([^|\s]+)", name)
            chrom, position_text, strand = key.rsplit(":", 2)
            position = int(position_text)
            sequence_center = int(center.group(1).rsplit(":", 1)[1]) if center else position
            block_start = (position // holdout_block_size) * holdout_block_size
            if sequence_center - len(seq) // 2 < block_start or sequence_center + (len(seq) + 1) // 2 > block_start + holdout_block_size:
                stats["sequence_center_outside_locus_holdout_block"] = stats.get("sequence_center_outside_locus_holdout_block", 0) + 1
                continue
            rank = stable_hash("discovery|" + key, seed)
            canonical_name = f"locus={key}" + (f"|center={center.group(1)}" if center else "")
            item = (-rank, key, canonical_name, seq)
            if len(heap) < max_sequences:
                heapq.heappush(heap, item)
            elif item > heap[0]:
                heapq.heapreplace(heap, item)
    records = [(item[2], item[3]) for item in sorted(heap, reverse=True)]
    stats.update({"discovery_sequences": len(records), "holdout_fraction": holdout_fraction,
                  "holdout_block_size": holdout_block_size, "holdout_guard": holdout_guard,
                  "seed": seed, "max_sequences": max_sequences})
    return records, stats


def motif_regex(motif):
    motif = motif.upper().replace("U", "T")
    if not motif or set(motif) - set(IUPAC):
        raise ValueError(f"Invalid IUPAC motif: {motif}")
    return re.compile("".join("[" + IUPAC[base] + "]" for base in motif))


def motif_matches(sequence, center_index, motif, anchor_pos, window):
    """Return all motif-anchor offsets, including overlapping matches.

    Positive offsets point to RNA 3'; anchor_pos is 1 based in the motif.
    center_index identifies the physically defined crosslink, never a motif snap.
    """
    if not 1 <= anchor_pos <= len(motif):
        raise ValueError("Motif anchor_pos must be within motif")
    pattern = motif_regex(motif)
    result = set()
    for offset in range(-window, window + 1):
        start = center_index + offset - anchor_pos + 1
        end = start + len(motif)
        if start >= 0 and end <= len(sequence) and pattern.fullmatch(sequence[start:end]):
            result.add(offset)
    return result


def bh_adjust(pvalues):
    indexed = sorted((p, i) for i, p in enumerate(pvalues) if p is not None and math.isfinite(p))
    result, running = [None] * len(pvalues), 1.0
    for rank0 in range(len(indexed) - 1, -1, -1):
        value, index = indexed[rank0]
        running = min(running, value * len(indexed) / (rank0 + 1))
        result[index] = running
    return result


def parse_meme_consensuses(path, probability_mass=0.8):
    """Read MEME DNA PWMs and return declared IUPAC consensus scan models.

    This is not a PWM/FIMO significance scan. At each position the smallest
    base set spanning probability_mass is used (ties included). Distributions
    are evaluated on genomic loci withheld from MEME discovery.
    """
    lines = Path(path).read_text().splitlines()
    motifs, motif_id = [], None
    for i, line in enumerate(lines):
        if line.startswith("MOTIF "):
            motif_id = line.split()[1]
        elif "letter-probability matrix:" in line:
            width = re.search(r"\bw\s*=\s*(\d+)", line)
            if width is None or motif_id is None:
                continue
            pwm = []
            for row in lines[i + 1:i + 1 + int(width.group(1))]:
                values = [float(value) for value in row.split()[:4]]
                if len(values) != 4 or sum(values) <= 0:
                    raise ValueError(f"Invalid MEME PWM in {path}")
                total = sum(values)
                pwm.append([value / total for value in values])
            if len(pwm) != int(width.group(1)):
                raise ValueError(f"Incomplete MEME PWM in {path}")
            consensus = ""
            for values in pwm:
                order = sorted(range(4), key=lambda j: (-values[j], j))
                chosen, cumulative = [], 0.0
                for j in order:
                    chosen.append(j)
                    cumulative += values[j]
                    if cumulative >= probability_mass:
                        threshold = values[j]
                        chosen.extend(k for k in order if k not in chosen and abs(values[k] - threshold) < 1e-9)
                        break
                consensus += IUPAC_BY_BASES[frozenset("ACGT"[j] for j in chosen)]
            if set(consensus) == {"N"}:
                continue
            motifs.append({"name": f"MEME_{motif_id}", "sequence": consensus,
                           "anchor_pos": (len(consensus) + 1) // 2,
                           "source": "de_novo_consensus_holdout", "probability_mass": probability_mass})
    return motifs
