from pymongo import MongoClient
import sqlite3
import pandas as pd
from pprint import pprint


def main():

    df = pd.merge(get_sra(), get_annot(), on="BioSample")

    apply_func(df, "BioProject", format_bioproject_accession)
    apply_func(df, "BioSample", format_biosample_accession)
    apply_func(df, "SRA Accession", format_sra_accession)
    apply_func(df, "Sample Description", format_sample_desc)
    apply_func(df, "Contact", lambda x: x[0] if len(x) > 0 else "")
    apply_func(df, "Papers", format_papers)

    (
        df.set_index(
            ["BioProject", "Project Title", "Project Description", "Submission Date"]
        ).to_excel("../output/notebook/2020-05-20_example_annot.xlsx")
    )


def get_sra():
    client = MongoClient()
    db = client["sramongo"]["ncbi"]

    df = (
        pd.DataFrame(
            db.aggregate(
                [
                    {"$match": {"library_strategy": "RNA-Seq"}},
                    {
                        "$group": {
                            "_id": "$BioSample.accn",
                            "BioProject": {"$first": "$BioProject.accn"},
                            "Project Title": {"$first": "$BioProject.title"},
                            "Project Description": {"$first": "$BioProject.description"},
                            "Submission Date": {"$first": "$sra_create_date"},
                            "SRA Accession": {"$addToSet": "$srx"},
                            "Sample Title": {"$first": "$BioSample.title"},
                            "Sample Description": {"$first": "$BioSample.attributes"},
                            "Papers": {"$first": "$papers"},
                            "Contact": {"$first": "$BioSample.contacts.email"},
                        }
                    },
                ],
                allowDiskUse=True
            )
        )
        .set_index("_id")
        .rename_axis("BioSample")
        .reset_index()
    )

    return df


def apply_func(df, col, func):
    df[col] = df[col].map(func)


def format_bioproject_accession(x):
    return f'=HYPERLINK("https://www.ncbi.nlm.nih.gov/bioproject/?term={x.upper()}", "{x.upper()}")'


def format_biosample_accession(x):
    return f'=HYPERLINK("https://www.ncbi.nlm.nih.gov/biosample/?term={x.upper()}", "{x.upper()}")'


def format_sra_accession(x):
    return ", ".join(
        [
            f'=HYPERLINK("https://www.ncbi.nlm.nih.gov/sra/?term={y.upper()}", "{y.upper()}")'
            for y in x
        ]
    )


def format_sample_desc(x):
    return "\n".join([f"{attr['name']}: {attr['value']}" for attr in x])


def format_papers(papers):
    if len(papers) == 0:
        return ""
    return ", ".join(
        set(
            [
                f'=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/{paper["accn"]}", "{paper["accn"]}")'
                for paper in papers
            ]
        )
    )


def get_annot():
    con = sqlite3.connect("../data/biometa.db")
    cur = con.cursor()
    cur.execute("SELECT * FROM biometa")
    df = pd.DataFrame(
        cur.fetchall(),
        columns=["BioSample", "sex", "devel_stage", "tissue", "cell_line", "perturbed", "complete"],
    ).drop("complete", axis=1)
    return df[df.tissue.str.contains("neuron")]


if __name__ == "__main__":
    main()
