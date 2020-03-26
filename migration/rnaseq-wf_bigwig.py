from pathlib import Path
import shutil

CLEAN_UP = False

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
    if CLEAN_UP:
        shutil.move(first, output_path / first.name)
        shutil.move(second, output_path / second.name)
    else:
        print(first, "-->", output_path / first.name)
        print(second, "-->", output_path / second.name)


if __name__ == "__main__":
    main()
