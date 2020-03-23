from pathlib import Path


def main():
    for pth in Path("../output/prealn-wf/samples").iterdir():
        if not pth.is_dir():
            continue

        remove_parquet(pth)

        if is_empty(pth):
            print(f"Removing:\t{pth}")
            pth.rmdir()


def remove_parquet(pth):
    for file_pth in pth.iterdir():
        if file_pth.is_dir():
            continue

        if file_pth.suffix == ".parquet":
            file_pth.unlink()

        print(file_pth)


def is_empty(pth):
    return len(list(pth.iterdir())) == 0


if __name__ == "__main__":
    main()
