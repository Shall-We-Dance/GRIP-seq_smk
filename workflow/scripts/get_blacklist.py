# workflow/scripts/get_blacklist.py
import os, gzip, shutil, urllib.request
from snakemake import snakemake  # just to satisfy some linters (not required)

genome = snakemake.config["genome"]
cache_dir = snakemake.config["blacklist"]["cache_dir"]
os.makedirs(cache_dir, exist_ok=True)

url_dict = {
    "mm10": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz",
    "mm39": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm39-blacklist.v2.bed.gz",
    "hg19": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz",
    "hg38": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz",
    "ce11": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce11-blacklist.bed.gz",
    "ce10": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/ce10-blacklist.bed.gz",
    "dm3": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm3-blacklist.bed.gz",
    "dm6": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/dm6-blacklist.bed.gz",
}
if genome not in url_dict:
    raise ValueError(f"Unsupported genome: {genome}. Available: {', '.join(url_dict.keys())}")

url = url_dict[genome]
gz_path = os.path.join(cache_dir, os.path.basename(url))
bed_path = snakemake.output.bed

if not os.path.exists(gz_path):
    print(f"Downloading blacklist for {genome}: {url}")
    urllib.request.urlretrieve(url, gz_path)

print(f"Unzipping blacklist to: {bed_path}")
with gzip.open(gz_path, "rb") as f_in, open(bed_path, "wb") as f_out:
    shutil.copyfileobj(f_in, f_out)
