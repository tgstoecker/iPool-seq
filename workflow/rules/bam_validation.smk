rule bam_validation:
    input:
        "results/alignment/{sample}.bam"
    output:
        "results/bam_validation/{sample}.invalid.txt"
#    params:
    conda:
        "../envs/picard.yml"
#in python "\" escapes so we add another "\" to escape that
# "|| true" at the end since output file is more than likely empty; this way force output status 0
    shell:
        "picard ValidateSamFile "
        "-INPUT {input} "
        "-MAX_OUTPUT 1000000000 "
        "-IGNORE RECORD_MISSING_READ_GROUP "
        "sed -n 's/^.*Read name \\([^,]*\\),.*$/\\1/p' | sort | uniq > {output} || true"
