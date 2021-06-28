rule expression_analysis:
    input:
        expand("results/kallisto/{sample}", sample=samples['sample']),
        tx2gene="resources/tx2gene.tsv"
    output:
        "results/diffexpr/tabela.DE.kallisto.txt"
    conda:
        "../envs/diffexpr.yaml"
    params:
        samples_file=samples_file
    script:
        "../scripts/kallisto_deseq2.R"