import importlib.util,sys,tempfile,unittest,gzip
from pathlib import Path
R=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('existing_tests',R/'tests/test_r2_signal.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
class LengthFilterTests(unittest.TestCase):
 def test_only_reference_aligned_query_bases_count(self):
  read=module.read
  values=[read(10,'19M'),read(10,'10M2I9M'),read(10,'10M100N10M'),read(10,'20M',True),read(10,'20M10S'),read(10,'2S20M')]
  with tempfile.TemporaryDirectory() as tmp:
   stats,args=module.ExtractionTests().run_signal(tmp,values,min_aligned_bases=20)
   self.assertEqual(stats['accepted_read2'],3)
   self.assertEqual(stats['excluded_R2_aligned_M_EQ_X_bases_below_minimum'],2)
   self.assertEqual(stats['excluded_five_prime_unaligned'],1)
   with gzip.open(args['out_crosslink_bed'],'rt') as f:rows=[r.split() for r in f]
   self.assertEqual([(int(r[1]),r[5],int(r[4])) for r in rows],[(9,'+',2),(30,'-',1)])
if __name__=='__main__':unittest.main()
