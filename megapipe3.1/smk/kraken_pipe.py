configfile: "./config/config_pipeline.json"
# edit this to change where megapipe gets the list of ids it processes
with open("/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/logs/test_50_runs/kraken_rerun.txt", "r") as file:
    samples = file.read().splitlines()

#I read the file with the of the samples to analyze to get the list of samples
rule all:
    input:
        config["path_ref_genome"]+".bwt",
        expand("/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/kraken/{sample}.krkn", sample=samples),
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
rule kraken:
    input:
        r1 = config["temp_dir"] + "{sample}/{sample}_1_trimmed.fastq.gz",
        r2 = config["temp_dir"] + "{sample}/{sample}_2_trimmed.fastq.gz"
    output:
        r1 = "/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/{sample}/kraken/{sample}.krkn",
    shell:
       "kraken --paired --fastq-input --gzip-compressed --db '/n/data1/hms/dbmi/farhat/bin/kraken/tbdb' {input.r1} {input.r2} > {output.r1}"
