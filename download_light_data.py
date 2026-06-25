from pathlib import Path
import requests
import urllib3


ROOT = Path(__file__).resolve().parents[1]
BASE = "https://genedev.bionet.nsc.ru/ftp/by_User/DashaPanchenko/OMICS_course_spring_2026"


def download(url, out):
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and out.stat().st_size > 0:
        print("skip", out)
        return
    tmp = out.with_suffix(out.suffix + ".part")
    print("download", out)
    with requests.get(url, stream=True, verify=False, timeout=60) as r:
        r.raise_for_status()
        with tmp.open("wb") as f:
            for chunk in r.iter_content(1024 * 1024):
                if chunk:
                    f.write(chunk)
    tmp.replace(out)


def main():
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    for sample in ["MoPh7", "MoPh11", "MoPh14", "MoPh15"]:
        download(
            f"{BASE}/day2/tracks/rnaseq/{sample}.rnaseq.STAR.bw",
            ROOT / f"course/day2_RNA_practice/results/tracks/rnaseq/{sample}.rnaseq.STAR.bw",
        )
        download(
            f"{BASE}/day2/tracks/cage/{sample}.cage.STAR.bw",
            ROOT / f"course/day2_RNA_practice/results/tracks/cage/{sample}.cage.STAR.bw",
        )
    download(
        f"{BASE}/day2/tracks/bwa/MoPh7.rnaseq.BWA.q30.bw",
        ROOT / "course/day2_RNA_practice/results/tracks/bwa/MoPh7.rnaseq.BWA.q30.bw",
    )
    download(
        f"{BASE}/genome/chm13v2.0_main_protein_coding_chrNames.gtf.gz",
        ROOT / "course/day5_omics_practice/data/annotation/chm13v2.0_main_protein_coding_chrNames.gtf.gz",
    )


if __name__ == "__main__":
    main()

