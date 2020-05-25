import shutil
from pathlib import Path


def main():
    if "prealn-wf" in snakemake.output[0]:
        prealn_cleanup()
    elif "rnaseq-wf" in snakemake.output[0]:
        rnaseq_cleanup()


def prealn_cleanup():
    sample_id = snakemake.params.sample_id  # SRR#####
    sample_folder = snakemake.params.sample_folder # ../output/prealn-wf/samples/SRR#####
    done_queue = snakemake.params.done_queue  # ../output/prealn-wf/done.txt

    remove_folder(sample_folder)
    append_sample_to_done_queue(sample_id, done_queue)
    touch_snakemake_done_file()


def rnaseq_cleanup():
    sample_id = snakemake.params.sample_id  # SRX######
    sample_folder = snakemake.params.sample_folder # ../output/rnaseq-wf/samples/SRX######
    done_queue = snakemake.params.done_queue  # ../output/rnaseq-wf/done.txt
    srrs = snakemake.params.srrs # [SRR######, SRR######,]

    remove_fastqs(srrs)
    keep_bams_in_rnaseq_folder(sample_folder)
    append_sample_to_done_queue(sample_id, done_queue)
    touch_snakemake_done_file()


def remove_file(file_name: str):
    if Path(file_name).exists():
        Path(file_name).unlink()


def remove_folder(folder_name: str):
    if Path(folder_name).exists():
        shutil.rmtree(folder_name)


def remove_fastqs(srrs: list):
    for srr in srrs:
        remove_file(f"../output/fastq-wf/fastqs/{srr}_1.fastq")
        remove_file(f"../output/fastq-wf/fastqs/{srr}_2.fastq")
        remove_file(f"../output/fastq-wf/sra_cache/{srr}.sra")


def keep_bams_in_rnaseq_folder(folder_name: str):
    for pth in Path(folder_name).iterdir():
        if pth.is_dir():
            remove_folder(pth.as_posix())
        elif pth.suffix in [".bam", ".bai"]:
            continue
        else:
            remove_file(pth.as_posix())


def append_sample_to_done_queue(sample_id: str, done_queue: str):
    with open(done_queue, "a") as fh:
        fh.write(sample_id + "\n")


def touch_snakemake_done_file():
    Path(snakemake.output[0]).touch()


if __name__ == "__main__":
    if "snakemake" not in locals() or not hasattr(snakemake, "scriptdir"):
        from ncbi_remap.mock import MockSnake

        snakemake = MockSnake(
            output="../output/test/prelaln-wf/done/SRX123456",
            params=dict(
                sample_id="SRX123456",
                sample_folder="../output/test/prelaln-wf/samples/SRX123456",
                done_queue="../output/test/prelaln-wf/done.txt",
            ),
        )

        Path(snakemake.output[0]).parent.mkdir(exist_ok=True, parents=True)
        Path(snakemake.params.sample_folder).mkdir(exist_ok=True, parents=True)

    main()
