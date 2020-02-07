"""Get things started."""
rule targets:
    input:
        "output/db_download.date",
        "output/prealn_queue.date",
        "output/srx2srr.csv"


rule sra2mongo:
    output: "output/db_download.date"
    shell: """
    sra2mongo --api $ENTREZ_API_KEY --query '"Drosophila melanogaster"[orgn]' && \
    date > {output[0]}
    """

rule srx2srr:
    output: "output/srx2srr.csv"
    script: "scripts/srx2srr.py"


rule initialize_prealn_queue:
    input: rules.sra2mongo.output[0]
    output: "output/prealn_queue.date"
    threads: 2
    shell: """
    ./prealn-wf/prealn-store.py init --append -j {threads} && \
    date > {output[0]}
    """


rule initialize_store:
    threads: 2
    shell: "./prealn-wf/prealn-store.py init -j {threads}"
