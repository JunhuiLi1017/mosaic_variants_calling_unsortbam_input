rule fastp:
    input:
        unpack(get_fastq)
    output:
        temp(["{outpath}/{sample}/01_multiqc/fastp/{sample}.R1.fastq.gz", "{outpath}/{sample}/01_multiqc/fastp/{sample}.R2.fastq.gz"]) if paired_end else temp("{outpath}/{sample}/01_multiqc/fastp/{sample}.R1.fastq.gz"),
        j="{outpath}/{sample}/01_multiqc/fastp/{sample}.json", 
        h="{outpath}/{sample}/01_multiqc/fastp/{sample}.html"
    log:
        "{outpath}/{sample}/01_multiqc/logs/{sample}.fastp.log"
    params:
        trim_expr= "" if config['trim'] == False else f"-f {trim_value} -F {trim_value}"
    run:
        if paired_end:
            shell("fastp {params.trim_expr} -i {input.r1} -I {input.r2} -o {output[0]} -O {output[1]} -j {output.j} -h {output.h} > {log} 2>&1")
        else:
            shell("fastp {params.trim_expr} -i {input.r1} -o {output[0]} -j {output.j} -h {output.h} > {log} 2>&1")

rule fastqc:
    input:
        get_clean_fastq
    output:
        ["{outpath}/{sample}/01_multiqc/fastqc/{sample}.R1_fastqc.html", "{outpath}/{sample}/01_multiqc/fastqc/{sample}.R2_fastqc.html"] if paired_end else "{outpath}/{sample}/01_multiqc/fastqc/{sample}.R1_fastqc.html",
        ["{outpath}/{sample}/01_multiqc/fastqc/{sample}.R1_fastqc.zip", "{outpath}/{sample}/01_multiqc/fastqc/{sample}.R2_fastqc.zip"] if paired_end else "{outpath}/{sample}/01_multiqc/fastqc/{sample}.R1_fastqc.zip"
    conda:
        "../envs/multiqc_env.yaml"
    log:
        "{outpath}/{sample}/01_multiqc/logs/{sample}.fastqc.log"
    params:
        out_fastqc="{outpath}/{sample}/01_multiqc/fastqc"
    shell:
        "fastqc -o {params.out_fastqc} {input} > {log} 2>&1"

rule multiqc:
    input:
        get_multiqc_input
    output:
        "{outpath}/{sample}/01_multiqc/multiqc_report.html"
    conda:
        "../envs/multiqc_env.yaml"
    log:
        "{outpath}/{sample}/01_multiqc/logs/multiqc.log"
    shell:
        """
        multiqc -o {outpath}/{sample}/01_multiqc {outpath}/{sample}/01_multiqc/fastqc {outpath}/{sample}/02_map/bqsr_stat/ --force
        """