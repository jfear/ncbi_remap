"""Get things started."""

rule sra2mongo:
    output: 'output/db_download.date'
    shell: """
    sra2mongo --api $ENTREZ_API_KEY --query "Drosophila melanogaster"[orgn] && \
    date > {output[0]}
    """