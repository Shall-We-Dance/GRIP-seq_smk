"""Exact sparse interval arithmetic for compatible BigWigs (no genome-sized arrays)."""
import heapq
import math
from pathlib import Path
import pyBigWig


def interval_events(reader, chrom, index):
    length = reader.chroms(chrom)
    for left in range(0, length, 1000000):
        right = min(left + 1000000, length)
        for start, end, value in reader.intervals(chrom, left, right) or ():
            if not math.isfinite(value):
                continue
            yield max(start, left), index, value
            yield min(end, right), index, -value


def combine_tracks(paths, destination, operation="mean", pseudocount=1.0):
    if not paths:
        raise ValueError("At least one input BigWig is required")
    if operation not in {"mean", "log2", "ratio", "subtract", "add"}:
        raise ValueError(f"Unsupported exact operation: {operation}")
    if operation != "mean" and len(paths) != 2:
        raise ValueError("A comparison needs exactly two input tracks")
    if operation in {"log2", "ratio"} and pseudocount <= 0:
        raise ValueError("A positive pseudocount is required")
    readers = []
    out = None
    try:
        readers = [pyBigWig.open(str(path)) for path in paths]
        chroms = readers[0].chroms()
        if any(reader.chroms() != chroms for reader in readers[1:]):
            raise ValueError("BigWig chromosome names and lengths must match exactly")
        Path(destination).parent.mkdir(parents=True, exist_ok=True)
        out = pyBigWig.open(str(destination), "w")
        out.addHeader(list(chroms.items()))
        for chrom, length in chroms.items():
            values = [0.0] * len(readers)
            events = heapq.merge(*(interval_events(r, chrom, i) for i, r in enumerate(readers)))
            pending = None
            buffer = []

            def emit(start, end):
                nonlocal pending
                if end <= start:
                    return
                if operation == "mean":
                    value = math.fsum(values) / len(values)
                elif operation == "subtract":
                    value = values[0] - values[1]
                elif operation == "add":
                    value = values[0] + values[1]
                else:
                    a, b = values[0] + pseudocount, values[1] + pseudocount
                    if a <= 0 or b <= 0:
                        raise ValueError("Log/ratio inputs plus pseudocount must be positive")
                    value = a / b
                    if operation == "log2":
                        value = math.log2(value)
                if abs(value) < 1e-12:
                    value = 0.0
                if pending and pending[1] == start and pending[2] == value:
                    pending = (pending[0], end, value)
                else:
                    if pending:
                        buffer.append(pending)
                    pending = (start, end, value)
                if len(buffer) >= 10000:
                    flush()

            def flush():
                if buffer:
                    out.addEntries([chrom] * len(buffer), [x[0] for x in buffer],
                                   ends=[x[1] for x in buffer], values=[float(x[2]) for x in buffer])
                    buffer.clear()

            previous = 0
            for position, index, delta in events:
                if position != previous:
                    emit(previous, position)
                    previous = position
                values[index] += delta
                if abs(values[index]) < 1e-12:
                    values[index] = 0.0
            emit(previous, length)
            if pending:
                buffer.append(pending)
            flush()
    finally:
        if out is not None:
            out.close()
        for reader in readers:
            reader.close()
