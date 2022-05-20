rule invalid_reads_removal:
    input:
        bam="results/alignment/{sample}.bam",
        inv="results/bam_validation/{sample}.invalid.txt",
    output:
#    params:
    conda:
        "../envs/picard.yml"
    #if invalid reads exist remove them and sort; else just sort
    shell:
        """
        if [ $$(wc -l < "$${NGM_INVAL}") -gt 0 ]; then \
          picard FilterSamReads \
          -INPUT {input.bam} \
          -OUTPUT ngm.filtered.bam \
          -READ_LIST_FILE ngm.invalid.txt \
          -FILTER excludeReadList \
          -VALIDATION_STRINGENCY SILENT \
          -WRITE_READS_FILES false \
          -SORT_ORDER queryname
        else \
          picard SortSam \
          -INPUT {input}.bam \
          -OUTPUT ngm.sorted.bam \
          -SORT_ORDER queryname
