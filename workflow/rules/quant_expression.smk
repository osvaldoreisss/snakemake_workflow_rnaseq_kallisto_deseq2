rule pseudomap_kallisto:
    input:
        fq1="results/concatenated/{sample}.R1.fq.gz",
        fq2="results/concatenated/{sample}.R2.fq.gz"
    output:
        directory("results/kallisto/{sample}")
    conda:
        "../envs/quantification.yaml"
    threads: threads
    params:    
        index=config["ref"]["index"],
        gtf=config['ref']['gtf'],
        extra=config['kallisto']['quant']['extra']
    log:
        "results/logs/kallisto/{sample}.log"
    shell:
        "kallisto quant -t {threads} {params.extra} " 
        "-g {params.gtf} "
        "-i {params.index} -o {output} {input.fq1} {input.fq2} > {log} 2>&1"