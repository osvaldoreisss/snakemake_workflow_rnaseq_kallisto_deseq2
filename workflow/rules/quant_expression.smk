rule create_kallisto_index:
    output:
        index="resources/kallisto_index",
        gtf="resources/transcripts.gtf.gz",
        cdna="resources/transcripts.fa.gz",
        tx2gene="resources/tx2gene.tsv"
    conda:
        "../envs/quantification.yaml"
    params:    
        gtf=config['ref']['gtf'],
        cdna=config['ref']['cdna']
    shell:
        """
        wget {params.cdna} -O {output.cdna}
        wget {params.gtf} -O {output.gtf}
        kallisto index -i {output.index} {output.cdna}
        zcat {output.gtf} | \
        awk '{{if($3=="transcript"){{print $10"\\t"$12}}}}' | \
        sed 's/;//g' | sed 's/\\"//g' > {output.tx2gene}
        """

rule pseudomap_kallisto:
    input:
        fq1="results/concatenated/{sample}.R1.fq.gz",
        fq2="results/concatenated/{sample}.R2.fq.gz",
        index="resources/kallisto_index",
        gtf="resources/transcripts.gtf.gz",
    output:
        directory("results/kallisto/{sample}")
    conda:
        "../envs/quantification.yaml"
    threads: threads
    params:    
        extra=config['kallisto']['quant']['extra']
    log:
        "results/logs/kallisto/{sample}.log"
    shell:
        "kallisto quant -t {threads} {params.extra} " 
        "-g {input.gtf} "
        "-i {input.index} -o {output} {input.fq1} {input.fq2} > {log} 2>&1"
