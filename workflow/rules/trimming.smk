#trimmomatic 1 - read through removal
rule read_though_removal:
    input:
        r1="rawreads/{sample}_1.fq.gz",
        r2="rawreads/{sample}_2.fq.gz"
    output:
        r1="results/read_through_removal/{sample}_1.fq.gz",
        r2="results/read_through_removal/{sample}_2.fq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="results/read_through_removal/{sample}_1.unpaired.fq.gz",
        r2_unpaired="results/read_through_removal/{sample}_2.unpaired.fq.gz"
    log:
        "logs/read_through_removal/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=config["RTA_trimmer"],
        # optional parameters
        extra="",
        compression_level="-9"
    threads:
        4
    resources:
        mem_mb=1024
    wrapper:
        "v1.5.0/bio/trimmomatic/pe"

#UMI extraction & technical sequence removal
rule UMIex_tech_seq_removal:
    input:
        r1="results/read_through_removal/{sample}_1.fq.gz",
        r2="results/read_through_removal/{sample}_2.fq.gz"
    output:
        r1="results/UMIex_tech_seq_removal/{sample}_1.fq.gz",
        r2="results/UMIex_tech_seq_removal/{sample}_2.fq.gz"
#    container:
#        "docker://tgstoecker/py2_ips_tir_tag"
    conda:
        "../envs/trim_tag_python2_environment.yml"
    threads: 4
    shell:
        "python workflow/scripts/trim.tag.py {threads} {input.r1} {input.r2} {output.r1} {output.r2}" 
