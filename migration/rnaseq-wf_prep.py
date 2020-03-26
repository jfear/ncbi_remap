"""Prep for RNA-Seq migration.

* Removes samples that are not in rnaseq.pkl

* Clears folder contents
    * ../output/rnaseq-wf/alignment_bad 
    * ../output/rnaseq-wf/atropos_bad
    * ../output/rnaseq-wf/done

* Creates directories
    * ../output/rnaseq-wf/atropos
    * ../output/rnaseq-wf/hisat2
    * ../output/rnaseq-wf/aln_stats
    * ../output/rnaseq-wf/{gene,junction,intergenic,segment,fusion}_counts
    * ../output/rnaseq-wf/{flybase,ucsc}_bigwigs
"""
import pickle
from pathlib import Path
import shutil


def main():
    remove_non_rnaseq()

    clear_folder("../output/rnaseq-wf/alignment_bad")
    clear_folder("../output/rnaseq-wf/atropos_bad")
    clear_folder("../output/rnaseq-wf/done")

    create_folder("../output/rnaseq-wf/atropos")
    create_folder("../output/rnaseq-wf/hisat2")
    create_folder("../output/rnaseq-wf/aln_stats")
    create_folder("../output/rnaseq-wf/gene_counts")
    create_folder("../output/rnaseq-wf/junction_counts")
    create_folder("../output/rnaseq-wf/intergenic_counts")
    create_folder("../output/rnaseq-wf/segment_counts")
    create_folder("../output/rnaseq-wf/fusion_counts")
    create_folder("../output/rnaseq-wf/flybase_bigwigs")
    create_folder("../output/rnaseq-wf/ucsc_bigwigs")

    cleanup_special_cases()

def remove_path(file_name):
    pth = Path(file_name)
    if pth.exists() and pth.is_file():
        pth.unlink()

def remove_folder(folder_name):
    pth = Path(folder_name)
    if pth.exists():
        shutil.rmtree(pth)

def remove_non_rnaseq():
    rnaseq = pickle.load(open("../output/library_strategy-wf/rnaseq.pkl", "rb"))
    workflow = {x.name for x in Path("../output/rnaseq-wf/samples").iterdir()}
    not_rnaseq = workflow - rnaseq
    for srx in not_rnaseq:
        remove_folder(f"../output/rnaseq-wf/samples/{srx}")


def clear_folder(folder_name):
    pth = Path(folder_name)
    for item in pth.iterdir():
        if item.is_file():
            item.unlink()
        else:
            print("This is not a file", item, sep="\t")


def create_folder(folder_name):
    pth = Path(folder_name)
    pth.mkdir(exist_ok=True)


def cleanup_special_cases():
    # Empty junction counts
    remove_folder("../output/rnaseq-wf/samples/SRX019645")



if __name__ == "__main__":
    main()
