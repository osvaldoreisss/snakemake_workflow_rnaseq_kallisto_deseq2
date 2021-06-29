rule get_seq:
    input: 
        unpack(get_fastq)
    output: 
        fq1="results/fastq/{sample}.{run}.R1.fq.gz",
        fq2="results/fastq/{sample}.{run}.R2.fq.gz"
    log: 
        "results/logs/get_seq/{sample}.{run}.log"
    shell:
        """
        gunzip < {input.fq1} > {output.fq1}
        gunzip < {input.fq2} > {output.fq2}
        """


rule trim_galore:
    input: 
        "results/fastq/{sample}.{run}.R1.fq.gz",
        "results/fastq/{sample}.{run}.R2.fq.gz"
    output: 
        directory("results/quality_analysis/{sample}.{run}")
    conda:
        "../envs/trim_galore.yaml"
    threads: threads
    params:
        **config["params"]
    log: 
        "results/logs/trim_galore/{sample}-{run}.log"
    shell:
        """
        trim_galore {params.trim} --cores {threads} --output_dir {output} {input} 2> {log}
        """

rule concatenate:
    input: 
        lambda wildcards: \
            [f"results/quality_analysis/{wildcards.sample}.{run}" \
                for run in samples.loc[(wildcards.sample), ["run"]]['run'].dropna()
            ]
    output: 
        fq1="results/concatenated/{sample}.R1.fq.gz",
        fq2="results/concatenated/{sample}.R2.fq.gz"
    threads: 1
    log:
        "results/logs/concatenate/{sample}.log"
    shell:
        """
        set +e
        FQ1=`echo {input} | awk '{{for(X=1;X<=NF;X++){{OUT=OUT $X"/*val_1.fq.gz "}}}}END{{print OUT}}'`
        FQ2=`echo {input} | awk '{{for(X=1;X<=NF;X++){{OUT=OUT $X"/*val_2.fq.gz "}}}}END{{print OUT}}'`
        echo $FQ1
        echo $FQ2
        cat $FQ1 > {output.fq1} 2> {log}
        cat $FQ2 > {output.fq2} 2>> {log}
        """