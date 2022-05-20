rule alignment:
    input:
        r1="results/UMIex_tech_seq_removal/{sample}_1.fq.gz",
        r2="results/UMIex_tech_seq_removal/{sample}_2.fq.gz",
        fa=config["FASTA"]
    output:
        "results/alignment/{sample}.bam"
#    params:
    conda:
        "../envs/alignment.yml"
    threads:
        8
    log: "logs/alignment/{sample}.log"
    shell:
        "ngm "
        "-r {input.fa} "
        "-p -1 {input.r1} -2 {input.r2} "
        "--bam -o {output} "
        "-t {threads} "
        "--end-to-end "
        "--pair-score-cutoff 0.5 "
        "--sensitivity 0.3 "
        "--kmer 13 "
        "--kmer-skip 0 " 
        "--skip-save"
