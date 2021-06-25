rule expression_analysis:
    input:
        expand("results/kallisto/{sample}", sample=samples['sample']) 
    output:
        "results/diffexpr/tabela.DE.kallisto.txt"
    conda:
        "../envs/diffexpr.yaml"
    params:
        samples_file=samples_file,
        transcripts_to_genes=config['ref']['tx2gene']
    script:
        "../scripts/kallisto_deseq2.R"