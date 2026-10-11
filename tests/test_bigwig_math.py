"""Synthetic signal arithmetic and reference-compatibility regression checks."""
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import pyBigWig
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/rules/scripts'))
from bigwig_math import combine_tracks


class BigWigMathTests(unittest.TestCase):
    def test_exact_mean_ratio_log_and_subtraction(self):
        with tempfile.TemporaryDirectory() as d:
            paths=[]
            arrays=[]
            for i, intervals in enumerate(([(1,4,2.),(6,8,4.)],[(2,5,6.),(8,10,2.)])):
                p=Path(d)/f'{i}.bw'; paths.append(p)
                a=np.zeros(12)
                with pyBigWig.open(str(p),'w') as b:
                    b.addHeader([('chrT',12)])
                    for start,end,value in intervals:
                        b.addEntries(['chrT'],[start],ends=[end],values=[value]); a[start:end]=value
                arrays.append(a)
            a,b=arrays
            for op,expected in [('mean',(a+b)/2),('log2',np.log2((a+1)/(b+1))),('ratio',(a+1)/(b+1)),('subtract',a-b)]:
                with self.subTest(operation=op):
                    out=Path(d)/f'{op}.bw'; combine_tracks(paths,out,op)
                    with pyBigWig.open(str(out)) as bw:
                        actual=np.nan_to_num(bw.values('chrT',0,12,numpy=True))
                        np.testing.assert_allclose(actual,expected,rtol=1e-6,atol=1e-6)
                        self.assertAlmostEqual(bw.stats('chrT',0,12,exact=True)[0],float(expected.mean()),places=6)

    def test_empty_inputs_and_reference_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError): combine_tracks([],Path(d)/'out.bw')
            paths=[]
            for i,size in enumerate((10,11)):
                p=Path(d)/f'{i}.bw'; paths.append(p)
                with pyBigWig.open(str(p),'w') as b: b.addHeader([('chrT',size)])
            with self.assertRaisesRegex(ValueError,'chromosome'):
                combine_tracks(paths,Path(d)/'out.bw')

    def test_contiguous_intervals_and_zero_gap(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.bw'
            with pyBigWig.open(str(p),'w') as b:
                b.addHeader([('chrT',10)])
                b.addEntries(['chrT']*3,[0,2,6],ends=[2,4,10],values=[1.,2.,3.])
            out=Path(d)/'out.bw'; combine_tracks([p],out)
            with pyBigWig.open(str(out)) as b:
                np.testing.assert_array_equal(np.nan_to_num(b.values('chrT',0,10,numpy=True)),[1,1,2,2,0,0,3,3,3,3])


if __name__=='__main__': unittest.main()
