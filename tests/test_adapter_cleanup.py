import json,subprocess,tempfile,unittest,shutil
from pathlib import Path
R=Path(__file__).resolve().parents[1]
class AdapterCleanupTests(unittest.TestCase):
 def test_residual_adapters_trim_without_moving_five_prime_endpoint(self):
  self.assertIsNotNone(shutil.which('cutadapt'), 'Install cutadapt to run adapter tests')
  a='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA';b='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
  rna='GCGGTTCAACTGATCCAGTTCGATCGTTAGC'
  pairs=[(rna+a,rna+b),(rna+a,'AACCGGTT'+b),(rna,rna),(rna+a[:8],rna+b[:8]),(rna+a[:15]+'A'+a[16:],rna+b[:15]+'A'+b[16:]),(rna+a[:13]+'TTTTTT',rna+b[:13]+'TTTTTT')]
  with tempfile.TemporaryDirectory() as tmp:
   tmp=Path(tmp)
   for j in range(2):
    (tmp/f'in{j}.fq').write_text(''.join(f'@read{i}\n{pair[j]}\n+\n'+('I'*len(pair[j]))+'\n' for i,pair in enumerate(pairs)))
   cmd=[shutil.which('cutadapt'),'-a',a,'-A',b,'-a','AGATCGGAAGAGC;min_overlap=13;max_error_rate=0','-A','AGATCGGAAGAGC;min_overlap=13;max_error_rate=0','--times','2','-e','0.1','-O','8','-m','20:20','--pair-filter=any','--json',str(tmp/'qc.json'),'-o',str(tmp/'out0.fq'),'-p',str(tmp/'out1.fq'),str(tmp/'in0.fq'),str(tmp/'in1.fq')]
   subprocess.run(cmd,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,check=True)
   for j in range(2):
    lines=(tmp/f'out{j}.fq').read_text().splitlines()
    self.assertEqual(lines[::4],['@read0','@read2','@read3','@read4','@read5'])
    self.assertEqual(lines[1::4],[rna]*5)
if __name__=='__main__':unittest.main()
