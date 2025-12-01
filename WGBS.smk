configfile : "config.yaml"

SAMPLES_NEST = config.get("samples", {})

SAMPLE2GROUP = {
    sample: group
    for group, sd in SAMPLES_NEST.items()
    for sample in sd.keys()
}

ALL_SAMPLES = SAMPLE2GROUP.keys()

ALL_GROUPS = set(SAMPLE2GROUP.values())

def fq_path(sample, read):
    grp = SAMPLE2GROUP.get(sample)
    return SAMPLES_NEST[grp][sample][read]

def grp_sample(group):
    SAMPLE=[]
    for sample in SAMPLE2GROUP:
        grp = SAMPLE2GROUP.get(sample)
        if grp==group:
            SAMPLE.append(sample)
    return SAMPLE

def make_output():
    output = []
    output += expand("data/{sample}/{sample}_R{Num}_val_{Num}_fastqc.html", sample=ALL_SAMPLES, Num=["1","2"])
    output.append("text/diff_simple_gr.txt")
    output += expand("data/multiqc_report/{sample}/multiqc_report.html", sample=ALL_SAMPLES)
    return output

rule all:
    input:
        make_output()

rule qc_before_filtered:
    output:
        html_R1=temp("qc/{sample}/{sample}_R1_fastqc.html"),
        html_R2=temp("qc/{sample}/{sample}_R2_fastqc.html"),
    params:
        FQ1 = lambda wc: fq_path(wc.sample, "R1"),
        FQ2 = lambda wc: fq_path(wc.sample, "R2"),
        threads=8,
    container:
        config["container"]
    shell:
        """
        fastqc -t {params.threads} {params.FQ1} {params.FQ2} -o qc/{wildcards.sample}
        """

rule trim_fastq:
    input:
        html_R1="qc/{sample}/{sample}_R1_fastqc.html",
        html_R2="qc/{sample}/{sample}_R2_fastqc.html",
    output:
        trimmed_FQ1=temp("qc/{sample}/{sample}_R1_val_1.fq.gz"),
        trimmed_FQ2=temp("qc/{sample}/{sample}_R2_val_2.fq.gz"),
    params:
        FQ1 = lambda wc: fq_path(wc.sample, "R1"),
        FQ2 = lambda wc: fq_path(wc.sample, "R2"),
    threads: 
        8
    container:
        config["container"]
    shell:
        """
        trim_galore --paired --fastqc --cores {threads} -o qc/{wildcards.sample} {params.FQ1} {params.FQ2}
        """

rule qc_after_filtered:
    input:
        trimmed_FQ1="qc/{sample}/{sample}_R1_val_1.fq.gz",
        trimmed_FQ2="qc/{sample}/{sample}_R2_val_2.fq.gz",
    output:
        html_R1="data/{sample}/{sample}_R1_val_1_fastqc.html",
        html_R2="data/{sample}/{sample}_R2_val_2_fastqc.html",
    threads: 
        8
    container:
        config["container"]
    shell:
        """
        fastqc -t {threads} qc/{wildcards.sample}/{wildcards.sample}_R1_val_1.fq.gz qc/{wildcards.sample}/{wildcards.sample}_R2_val_2.fq.gz -o data/{wildcards.sample}
        """

rule bitmapperBS_search:
    input:
        trimmed_FQ1="qc/{sample}/{sample}_R1_val_1.fq.gz",
        trimmed_FQ2="qc/{sample}/{sample}_R2_val_2.fq.gz",
    output:
        WGBS_bitmapperbs_bam=temp("mapped/{sample}.WGBS.bitmapperbs.bam"),
    params:
        BITMAPPERBSIDX=config["BITMAPPERBSIDX"],
    container:
        config["container"]
    shell:
        """
        bitmapperBS --search {params.BITMAPPERBSIDX} --sensitive -e 0.1 --seq1 qc/{wildcards.sample}/{wildcards.sample}_R1_val_1.fq.gz --seq2 qc/{wildcards.sample}/{wildcards.sample}_R2_val_2.fq.gz --bam -o {output.WGBS_bitmapperbs_bam}
        """

rule samtools_sort_bam_n:
    input:
        WGBS_bitmapperbs_bam="mapped/{sample}.WGBS.bitmapperbs.bam",
    output:
        WGBS_bitmapperbs_namesort_bam=temp("mapped/{sample}.WGBS.bitmapperbs.namesort.bam"),
    container:
        config["container"]
    shell:
        """
        samtools sort -n {input.WGBS_bitmapperbs_bam} > {output.WGBS_bitmapperbs_namesort_bam}
        """

rule samtools_fixmate_bam:
    input:
        WGBS_bitmapperbs_namesort_bam="mapped/{sample}.WGBS.bitmapperbs.namesort.bam",
    output:
        WGBS_bitmapperbs_fixmate_bam=temp("mapped/{sample}.WGBS.bitmapperbs.fixmate.bam"),
    container:
        config["container"]
    shell:
        """
        samtools fixmate -m {input.WGBS_bitmapperbs_namesort_bam} mapped/{wildcards.sample}.WGBS.bitmapperbs.fixmate.bam
        """

rule samtools_sort_bam:
    input:
        WGBS_bitmapperbs_fixmate_bam="mapped/{sample}.WGBS.bitmapperbs.fixmate.bam",
    output:
        WGBS_bitmapperbs_sort_bam="mapped/{sample}.WGBS.bitmapperbs.sort.bam",
    container:
        config["container"]
    shell:
        """
        samtools sort {input.WGBS_bitmapperbs_fixmate_bam} > {output.WGBS_bitmapperbs_sort_bam}
        """

rule samtools_index_bam:
    input:
        WGBS_bitmapperbs_sort_bam="mapped/{sample}.WGBS.bitmapperbs.sort.bam",
    output:
        WGBS_bitmapperbs_sort_bam_index=temp("mapped/{sample}.WGBS.bitmapperbs.sort.bam.bai"),
    container:|l
        config["container"]
    shell:
        """
        samtools index {input.WGBS_bitmapperbs_sort_bam}
        """

rule samtools_markdup_bam:
    input:
        WGBS_bitmapperbs_sort_bam_index="mapped/{sample}.WGBS.bitmapperbs.sort.bam.bai",
    output:
        WGBS_filter_bitmapperbs_bam="mapped/{sample}.WGBS_filter.bitmapperbs.bam",
    threads:
        8
    params:
        WGBS_bitmapperbs_sort_bam="mapped/{sample}.WGBS.bitmapperbs.sort.bam",
    container:
        config["container"]
    shell:
        """
        samtools markdup -r -@ {threads} {params.WGBS_bitmapperbs_sort_bam} mapped/{wildcards.sample}.WGBS_filter.bitmapperbs.bam
        """

rule samtools_index_bam2:
    input:
        WGBS_filter_bitmapperbs_bam="mapped/{sample}.WGBS_filter.bitmapperbs.bam",
    output:
        WGBS_filter_bitmapperbs_bam_index=temp("mapped/{sample}.WGBS_filter.bitmapperbs.bam.bai"),
    container:
        config["container"]
    shell:
        """
        samtools index {input.WGBS_filter_bitmapperbs_bam}
        """

rule MethylDackel_mbias:
    input:
        WGBS_filter_bitmapperbs_bam_index="mapped/{sample}.WGBS_filter.bitmapperbs.bam.bai",
    output:
        mbias=temp("text/mBias_table_MethylDackel_{sample}.txt"),
        OT_OB=temp("text/{sample}_options.log"),
    threads:
        4
    params:
        WGBS_filter_bitmapperbs_bam="mapped/{sample}.WGBS_filter.bitmapperbs.bam",
        BITMAPPERBSIDX=config["BITMAPPERBSIDX"],
    container:
        config["container"]
    shell:
        """
        MethylDackel mbias -@ {threads} --txt {params.BITMAPPERBSIDX} {params.WGBS_filter_bitmapperbs_bam} text/{wildcards.sample}.WGBS.bitmapperbs.mbias > text/mBias_table_MethylDackel_{wildcards.sample}.txt 2>text/{wildcards.sample}_options.log
        """

rule MethylDackel_extract:
    input:
        mbias="text/mBias_table_MethylDackel_{sample}.txt",
        OT_OB="text/{sample}_options.log",
    output:
        WGBS_filter_bitmapperbs_CpG_bedGraph=temp("mapped/{sample}.WGBS_filter.bitmapperbs_CpG.bedGraph"),
    params:
        WGBS_filter_bitmapperbs_bam="mapped/{sample}.WGBS_filter.bitmapperbs.bam",
        BITMAPPERBSIDX=config["BITMAPPERBSIDX"],
    container:
        config["container"]
    shell:
        """
        param=`cat text/{wildcards.sample}_options.log|sed \'s/Suggested inclusion options: //g\'|sed \'s/,0/,1\\/150/g\'|sed \'s/--OT 0,/--OT 1\\/150,/g\'|sed \'s/--OB 0,/--OB 1\\/150,/g\'`
        MethylDackel extract --mergeContext {params.BITMAPPERBSIDX} $param {params.WGBS_filter_bitmapperbs_bam}
        """

rule methylSig:
    input:
        expand(rules.MethylDackel_extract.output.WGBS_filter_bitmapperbs_CpG_bedGraph, sample=ALL_SAMPLES),
    output:
        diff_simple_gr="text/diff_simple_gr.txt",
    container:
        config["container"]
    shell:
        """
        Rscript scripts/methylSig_WGBS.r samplesheet.txt
        """

rule Final_MultiQC:
    input:
        mapped_bam_index="mapped/{sample}.WGBS_filter.bitmapperbs.bam.bai",
    output:
        multiqc_report_html="data/multiqc_report/{sample}/multiqc_report.html",
    container:
        config["container"] 
    shell:
        """
        multiqc data/ -o data/multiqc_report/{wildcards.sample} --force
        """