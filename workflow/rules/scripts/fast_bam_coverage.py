"""deepTools bamCoverage with vectorized run-length output for 1-bp scaling.

Counting, read filters, CPM denominator, parallel chunking, blacklist handling,
and BigWig serialization remain deepTools operations. The scalar output loop and redundant ndarray/list round-trip are replaced;
the read-counting method itself is unchanged. Other modes fall back.
"""
import numpy as np
from deeptools import writeBedGraph as wg
from deeptools import utilities
ORIGINAL = wg.WriteBedGraph.writeBedGraph_worker

def write_worker(self, chrom, start, end, func_to_call, func_args, bed_regions_list=None):
    if (self.binLength != 1 or self.smoothLength or self.skipZeroOverZero
            or func_to_call is not wg.scaleCoverage):
        return ORIGINAL(self, chrom, start, end, func_to_call, func_args, bed_regions_list)
    if start > end:
        raise ValueError('start exceeds end')
    # In this exact single-BAM, contiguous 1-bp mode, deepTools calls this
    # same counting routine, then needlessly converts every base to a Python
    # list and back to an ndarray. Retain its ndarray without that round-trip.
    if (len(getattr(self, 'bamFilesList', [])) == 1 and self.stepSize == 1
            and not self.save_data and bed_regions_list is None):
        from deeptools import bamHandler
        with bamHandler.openBam(self.bamFilesList[0]) as bam:
            counts = self.get_coverage_of_region(bam, chrom, [(start, end, 1)])
        coverage = np.asarray(counts, dtype='float64').reshape(-1, 1)
    else:
        coverage, _ = self.count_reads_in_region(chrom, start, end)
    values = func_args['scaleFactor'] * coverage[:, 0]
    filename = utilities.getTempFileName(suffix='.bg')
    with open(filename, 'w') as stream:
        if len(values):
            cuts = np.r_[0, np.flatnonzero(values[1:] != values[:-1]) + 1, len(values)]
            for a, b in zip(cuts[:-1], cuts[1:]):
                value = values[a]
                if not np.isnan(value) and start + int(a) != end:
                    stream.write(f'{chrom}\t{start + int(a)}\t{min(start + int(b), end)}\t{value:g}\n')
    return chrom, start, end, filename

def main():
    wg.WriteBedGraph.writeBedGraph_worker = write_worker
    from deeptools.bamCoverage import main as bamcoverage_main
    bamcoverage_main()

if __name__ == '__main__':
    main()
