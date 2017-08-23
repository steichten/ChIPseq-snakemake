#hard-coded location of things - TO EDIT AS NEEDED
SAMPLES = glob_wildcards("raw/{S}_1.fastq.gz").S
REFERENCE = "/home/steve/genomes/TAIR10"

#software requirements
#subread
#bamtools
#fastqc
#trim_galore (cutadapt)


#bash safe mode
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

rule final:
    input:
        expand("cleaned/{sample}_1.fastq.gz", sample=SAMPLES),
        expand("cleaned/{sample}_2.fastq.gz", sample=SAMPLES),
        expand("qc/{sample}_1_fastqc.html", sample=SAMPLES),
        expand("qc/{sample}_2_fastqc.html", sample=SAMPLES),
        expand("out/{sample}.sorted.bam",sample=SAMPLES),
        expand("out/{sample}.sorted.bam.bai",sample=SAMPLES),

rule fastqc:
    input:
        r1="raw/{sample}_1.fastq.gz",
        r2="raw/{sample}_2.fastq.gz",
    output:
        q1="qc/{sample}_1_fastqc.html",
        q2="qc/{sample}_2_fastqc.html",
    shell:
        """
        fastqc {input.r1}
        fastqc {input.r2}
        cd raw
        rm {wildcards.sample}_1_fastqc.zip
        rm {wildcards.sample}_2_fastqc.zip
        mv {wildcards.sample}_1_fastqc.html ../{output.q1}
        mv {wildcards.sample}_2_fastqc.html ../{output.q2}
        """

rule trim:
    input:
        r1="raw/{sample}_1.fastq.gz",
        r2="raw/{sample}_2.fastq.gz",
    output:
        r1="cleaned/{sample}_1.fastq.gz",
        r2="cleaned/{sample}_2.fastq.gz",
        l1="logs/{sample}_1.fastq.gz_trimming_report.txt",
        l2="logs/{sample}_2.fastq.gz_trimming_report.txt",
    shell:
        """
        trim_galore --paired {input.r1} {input.r2}
        mv {wildcards.sample}_1_val_1.fq.gz {output.r1}
        mv {wildcards.sample}_2_val_2.fq.gz {output.r2}
        mv {wildcards.sample}_1.fastq.gz_trimming_report.txt {output.l1}
        mv {wildcards.sample}_2.fastq.gz_trimming_report.txt {output.l2}
        """

rule alignment:
    input:
        r1="cleaned/{sample}_1.fastq.gz",
        r2="cleaned/{sample}_2.fastq.gz",
    output:
        bam="out/{sample}.sorted.bam",
        bai="out/{sample}.sorted.bam.bai",
    params:
        ref=REFERENCE,
    shell:
        """
        subread-align -t 1 -T 8 -P 3 -d 50 -D 600 -i {params.ref} -r {input.r1} -R {input.r2} -o {output.bam}
        # -t 1 = DNA-seq data
        # -T 8 = 8 threads
        # -P 3 = Phread+33
        # d 50 = min frag length 50bp
        # -D 600 = max frag length 600bp
        samtools sort {output.bam}
        bamtools index -in {output.bam}
        """


#macs2 callpeak -t {input.chip} -c {input.in} -f BAM -g mm -n {sample} -B -q 0.01 --outdir out
#Rscript out/{sample}_model.r
#for FILE in *.bdg; do sort -k1,1 -k2,2n $FILE > ${FILE}.sort.bdg; done
#for FILE in $sort.bdg; do bedGraphToBigWig $FILE <chrom.sizes.file> {sample}.bigwig; done
