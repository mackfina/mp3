configfile: "./config/config_pipeline.json"
# edit this to change where megapipe gets the list of ids it processes
with open(config["logs_analysis"]+"missing_runs.txt", "r") as file:
    samples = file.read().splitlines()
    
#I read the file with the of the samples to analyze to get the list of samples
rule all:
    input:
        config["path_ref_genome"]+".bwt",
        expand("/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/bam/{sample}.duprem.bam.bai", sample=samples),
        expand("/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/fast-lineage-caller/{sample}.lineage", sample=samples),
        expand("/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/spades/scaffolds.fasta", sample=samples)

rule generate_idx_ref:
    input: config["path_ref_genome"]
    output: config["path_ref_genome"]+".bwt"
    log: config["logs_analysis"]+"indexed_ref.txt"
    shell: "bwa index {input} 2> {log}"

#make the trimmed output temp before running the real thing; temp(...) around the whole thing
rule trim_with_trimmomatic:
    input:
        r1 = config["fastq_dir"] + "{sample}/IlluminaWGS/FASTQs/{sample}_1.fastq.gz",
        r2 = config["fastq_dir"] + "{sample}/IlluminaWGS/FASTQs/{sample}_2.fastq.gz",
    output:
        r1 =temp(config["temp_dir"] + "{sample}/{sample}_1_trimmed.fastq.gz"),
        r2 =temp(config["temp_dir"] + "{sample}/{sample}_2_trimmed.fastq.gz"),
        # reads where trimming entirely removed the mate
        r1_unpaired = temp(config["temp_dir"] + "{sample}/{sample}_1_trimmed.unpaired.fastq.gz"),
        r2_unpaired = temp(config["temp_dir"] + "{sample}/{sample}_2_trimmed.unpaired.fastq.gz"),
    log:
        config["logs_analysis"]+"{sample}/trimmomatic.txt"
    params:
        # list of trimmers (see manual) http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
        trimmer=["ILLUMINACLIP:'data/adapters/adapters_combined_256.fasta':2:30:10:2:true SLIDINGWINDOW:4:20 MINLEN:75"],
    threads: 8
    wrapper:
        "0.38.0/bio/trimmomatic/pe"

rule spades_assemble:
    input:
        read1=config["temp_dir"] + "{sample}/{sample}_1_trimmed.fastq.gz",
        read2=config["temp_dir"] + "{sample}/{sample}_2_trimmed.fastq.gz"
    output:
        asm1="/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/spades/scaffolds.fasta",
        dir="/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/spades",
    log:
        config["logs_analysis"]+"{sample}/spades.txt"
    shell:
       "spades.py --isolate -1 {input.read1} -2 {input.read2} -o {output.dir}"        

rule align_to_ref:
    input:
        read1=config["temp_dir"] + "{sample}/{sample}_1_trimmed.fastq.gz",
        read2=config["temp_dir"] + "{sample}/{sample}_2_trimmed.fastq.gz"
    output:
        temp(config["temp_dir"]+ "{sample}/{sample}.sam")
    log:
        config["logs_analysis"]+"{sample}/align_to_ref.txt"
    shell:
        """
            minimap2 -ax sr {config[path_ref_genome]} {input.read1} {input.read2} > {output} 2> {log}
        """

rule sort_convert_tobam:
    input:
        config["temp_dir"]+ "{sample}/{sample}.sam"
    output:
        temp(config["temp_dir"]+ "{sample}/{sample}.bam")
    log:
        config["logs_analysis"]+"{sample}/sort_convert_tobam.txt"
    shell:
        """
        java -Xmx16G -jar {config[path_to_picard]} SortSam INPUT={input} OUTPUT={output} SORT_ORDER=coordinate > {log} 2>&1
        """

rule remove_duplicates:
    input:
        config["temp_dir"]+ "{sample}/{sample}.bam"
    output:
        outfile="/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/bam/{sample}.duprem.bam", metrics=config["temp_dir"]+"results/{sample}/bam/{sample}.metrics"
    log:
        config["logs_analysis"]+"{sample}/duprem.txt"
    shell:
        "java -Xmx32G -jar {config[path_to_picard]} MarkDuplicates I={input} O={output.outfile} REMOVE_DUPLICATES=true M={output.metrics} ASSUME_SORT_ORDER=coordinate > {log} 2>&1"

rule calculate_depth:
    input:
        "/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/bam/{sample}.duprem.bam"
    output:
        "/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/depth/{sample}.depth.gz", "/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/depth/{sample}_depth_OK"
    log:
        config["logs_analysis"]+"{sample}/calc_depth.txt"
    shell:
        """
        samtools depth -a {input} > /n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{wildcards.sample}/depth/{wildcards.sample}.depth 2> {log}
        ./bin/check_depth.py /n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{wildcards.sample}/depth/{wildcards.sample}.depth >> {log} 2>&1
        gzip /n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{wildcards.sample}/depth/{wildcards.sample}.depth >> {log} 2>&1
        """

rule indexing_bam:
    input:
        bam="/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/bam/{sample}.duprem.bam", depth_ok="/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/depth/{sample}_depth_OK"
    output:
        "/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/bam/{sample}.duprem.bam.bai"
    log:
        config["logs_analysis"]+"{sample}/indexing_bam.txt"
    shell:
        "samtools index {input.bam} > {log} 2>&1"

rule variant_calling:
    input:
        bam="/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/bam/{sample}.duprem.bam",depth_ok="/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/depth/{sample}_depth_OK",  bai="/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/bam/{sample}.duprem.bam.bai"
    output:
        "/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/pilon/{sample}.vcf", "/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/pilon/{sample}.fasta", "/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/pilon/{sample}_full.vcf.gz"
    log:
        config["logs_analysis"]+"{sample}/variant_calling.txt"
    shell:
        """
        pilon -Xmx32G --genome {config[path_ref_genome]} --bam {input.bam} --output /n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{wildcards.sample}/pilon/{wildcards.sample} --variant > {log} 2>&1
        mv /n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{wildcards.sample}/pilon/{wildcards.sample}.vcf /n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{wildcards.sample}/pilon/{wildcards.sample}_full.vcf >> {log} 2>&1
        ./bin/vcf_cutter.py /n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{wildcards.sample}/pilon/{wildcards.sample}_full.vcf /n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{wildcards.sample}/pilon/{wildcards.sample}.vcf >> {log} 2>&1
        gzip /n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{wildcards.sample}/pilon/{wildcards.sample}_full.vcf >> {log} 2>&1
        """

rule lineage_calling:
    input:
        "/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/pilon/{sample}.vcf"
    output:
        "/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/fast-lineage-caller/{sample}.lineage"
    log:
        config["logs_analysis"]+"{sample}/lineage_calling.txt"
        #Ideally could have everything we install into SW put into the megapipe repo and downloaded concurrently
    shell:
        "fast-lineage-caller-vcf {input} /home/maf7596/sw/fast-lineage-caller-vcf/snp_schemes/coll.tsv > {output} 2> {log}"

rule rem_temp_files:
    shell:
        """
        cd {config[temp_dir]}
        printf "current directory: ";pwd
        rm -rf -I *
        """

rule rem_results:
    shell:
        """
        cd results
        printf "current directory: ";pwd
        rm -rf -I *
        """
