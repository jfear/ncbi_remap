import os
from pathlib import Path
import shutil

CLEAN_UP = os.environ.get("CLEAN_UP", False)

def main():
    for srx_pth in Path("../output/rnaseq-wf/samples").iterdir():
        if not srx_pth.is_dir():
            continue

        srx = srx_pth.name

        # UCSC
        first = srx_pth / f"{srx}.first.bw"
        second = srx_pth / f"{srx}.second.bw"
        move(first, second, "ucsc_bigwigs")

        # FlyBase
        first = srx_pth / f"{srx}.flybase.first.bw"
        second = srx_pth / f"{srx}.flybase.second.bw"
        move(first, second, "flybase_bigwigs")


def move(first: Path, second: Path, dir_name: str):
    output_path = Path("../output/rnaseq-wf") / dir_name
    output_path.mkdir(exist_ok=True)
    first_out = output_path / first.name
    second_out = output_path / second.name

    if first.exists() and second.exists() and CLEAN_UP:
        shutil.move(first, first_out)
        shutil.move(second, second_out)
    elif first.exists() and second.exists():
        print(first, "-->", first_out)
        print(second, "-->", second_out)
    else:
        print("Missing:", first, second, sep="\t")


if __name__ == "__main__":
    main()
