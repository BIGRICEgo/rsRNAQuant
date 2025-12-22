import os
import json

# Config File Path
CONFIGFILE = "config/config.yaml"
configfile: CONFIGFILE

# Global variables
RAW_DATA_DIR = config["raw_data_dir"]
RESULTS_DIR = config["results_dir"]
SCRIPTS_DIR = config["scripts_dir"]
ANNOTATIONS_DIR = config["annotations_dir"]
RRNA_LIST = ["5S", "5.8S", "12S", "16S", "18S", "28S", "45S"]
NCRNA_LIST = ["nuclear_tRNA", "mt_RNA", "yRNA", "miRNA", "snoRNA", "snRNA", "vaultRNA", "piRNA"]

# Sample related variables
GSE = config["project"]
SAMPLES = config["samples"]
GROUPS = config["groups"]

# check whether to run ncRNA mapping
NCRNA_FLAG = str(config.get("ncRNA", False)).strip().lower() in ("true", "1")
KEEP_FILES_FLAG = str(config.get("intermediate_file", False)).strip().lower() in ("true", "1")


rule all:
    input:
        f"{RESULTS_DIR}/{GSE}/stats/{GSE}_alignment_summary.csv",
        f"{RESULTS_DIR}/{GSE}/stats/{GSE}_rRNA_count_matrix.txt",
        f"{RESULTS_DIR}/{GSE}/stats/{GSE}_rRNA_expression_matrix.txt",
        expand(f"{RESULTS_DIR}/{GSE}/plots/{GSE}_{{rRNA}}_Coverage.png", rRNA=RRNA_LIST),
        f"{RESULTS_DIR}/{GSE}/plots/{GSE}_rRNA_Length_Distribution.png",
        f"{RESULTS_DIR}/{GSE}/plots/{GSE}_rRNA_subtype_HollowPieChart.png",
        f"{RESULTS_DIR}/{GSE}/plots/{GSE}_rRNA_End_StackBarPlot.png",
        f"{RESULTS_DIR}/{GSE}/{GSE}_report.html",
        f"{RESULTS_DIR}/{GSE}/logs/.all.done"


rule trim_reads:
    input:
        raw_data=f"{RAW_DATA_DIR}/{GSE}/{{sample}}.fastq.gz"
    output:
        trimmed=f"{RESULTS_DIR}/{GSE}/trimmed/{{sample}}_trimmed.fq.gz",
        report=f"{RESULTS_DIR}/{GSE}/trimmed/{{sample}}.fastq.gz_trimming_report.txt"
    log:
        f"{RESULTS_DIR}/{GSE}/logs/{{sample}}_trim.log"
    shell:
        """
        trim_galore --quality 30 --length 15 --max_length 50 -j 8 \
                    --output_dir {RESULTS_DIR}/{GSE}/trimmed/ \
                    --gzip {input.raw_data} > {log} 2>&1
        
        echo "Successfully completed trimming reads." >> {log}
        """

# Custom step: Filter reads against UniVec database
rule filter_reads:
    input:
        trimmed=f"{RESULTS_DIR}/{GSE}/trimmed/{{sample}}_trimmed.fq.gz"
    output:
        filtered=f"{RESULTS_DIR}/{GSE}/filtered/{{sample}}_univec_unmapped.fastq",
        stats=f"{RESULTS_DIR}/{GSE}/filtered/{{sample}}_match_univec.bowtie.stat"
    log:
        f"{RESULTS_DIR}/{GSE}/logs/{{sample}}_filter.log"
    shell:
        """
        bowtie -p 4 -x {ANNOTATIONS_DIR}/UniVec/UniVec \
               --un {output.filtered} -q {input.trimmed} \
               > /dev/null 2> {output.stats}
        
        echo "Successfully completed filtering reads." >> {log}
        """

rule count_unique_sequences:
    input:
        filtered=f"{RESULTS_DIR}/{GSE}/filtered/{{sample}}_univec_unmapped.fastq"
    output:
        unique_fa=f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_unique.fa",
        len_dist=f"{RESULTS_DIR}/{GSE}/filtered/{{sample}}_len_dist.txt"
    params:
        tmp_seq=f"{RESULTS_DIR}/{GSE}/filtered/{{sample}}_tmp.seq",
        tmp_count=f"{RESULTS_DIR}/{GSE}/filtered/{{sample}}_tmp.count"
    log:
        f"{RESULTS_DIR}/{GSE}/logs/{{sample}}_count_unique.log"
    shell:
        """
        bash {SCRIPTS_DIR}/count_unique_sequences.sh \
            {input.filtered} {output.unique_fa} {output.len_dist} \
            {params.tmp_seq} {params.tmp_count} {log}
        
        gzip {input.filtered}
        """

rule map_genome:
    input:
        unique_fa=f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_unique.fa"
    output:
        mapped=f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_genome_mapped.fa",
        unmapped=f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_genome_unmapped.fa",
        stats=f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_match_genome.bowtie.stat"
    log:
        f"{RESULTS_DIR}/{GSE}/logs/{{sample}}_map_genome.log"
    shell:
        """
        bowtie -v 1 -p 10 -x {ANNOTATIONS_DIR}/genome/hg38/genome -k 1 \
               --al {output.mapped} --un {output.unmapped} \
               -f {input.unique_fa} > /dev/null 2> {output.stats}
        
        echo "Successfully completed mapping to genome." >> {log}
        """

rule map_rRNA:
    input:
        genome_mapped = f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_genome_mapped.fa"
    output:
        match = f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_MG_match_{{rRNA}}_rRNA.fa",
        unmatch = f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_MG_unmatch_{{rRNA}}_rRNA.fa",
        stats = f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_MG_match_{{rRNA}}_rRNA.bowtie.stat",
        detail = f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_detail_match_{{rRNA}}_rRNA"    
    log:
        rRNA_log = f"{RESULTS_DIR}/{GSE}/logs/{{sample}}_map_rRNA_{{rRNA}}.log"
    threads: 8
    params:
        script_rRNA = f"{SCRIPTS_DIR}/map_rRNA.sh",
        output_dir = f"{RESULTS_DIR}/{GSE}"
    shell:
        r"""
            bash {params.script_rRNA} {wildcards.sample} {params.output_dir} {ANNOTATIONS_DIR} {threads} > {log.rRNA_log} 2>&1
            echo "Successfully completed mapping to rRNA." >> {log.rRNA_log}
        """

rule align_summary:
    input:
        rRNA_log = expand(f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_detail_match_{{rRNA}}_rRNA", sample=SAMPLES, rRNA=RRNA_LIST),
        ncRNA_log = expand(f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_detail_match_{{ncRNA}}_ncRNA", sample=SAMPLES, ncRNA=NCRNA_LIST if NCRNA_FLAG else [])
    output:
        summary = f"{RESULTS_DIR}/{GSE}/stats/{GSE}_alignment_summary.csv",
    log:
        f"{RESULTS_DIR}/{GSE}/logs/{GSE}_align_summary.log"
    params:
        script = f"{SCRIPTS_DIR}/align_summary.sh",
        output_dir = f"{RESULTS_DIR}/{GSE}",
        sample_list = " ".join(SAMPLES),
        ncrna_flag = "true" if NCRNA_FLAG else "false"
    shell:
        r"""
        bash {params.script} {GSE} {params.output_dir} "{params.sample_list}" {params.ncrna_flag} > {log} 2>&1
        """

rule process_rna_seq:
    input:
        unique_fa = expand(f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_unique.fa", sample=SAMPLES),
        rRNA_detail = expand(f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_detail_match_{{rRNA}}_rRNA", sample=SAMPLES, rRNA=RRNA_LIST),
        ncRNA_log = expand(f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_detail_match_{{ncRNA}}_ncRNA", sample=SAMPLES, ncRNA=NCRNA_LIST if NCRNA_FLAG else [])
    output:
        output_txt = f"{RESULTS_DIR}/{GSE}/stats/{{sample}}_output.txt",
        summary_txt = f"{RESULTS_DIR}/{GSE}/stats/{{sample}}_summary.txt", 
        length_txt = f"{RESULTS_DIR}/{GSE}/stats/{{sample}}_length_distribution.txt",
        rRNA_output_txt = f"{RESULTS_DIR}/{GSE}/stats/{{sample}}_rRNA_output.txt"
    log:
        f"{RESULTS_DIR}/{GSE}/logs/{{sample}}_process_rna_seq.log"
    params:
        output_dir = f"{RESULTS_DIR}/{GSE}",
        ncrna_flag = "true" if NCRNA_FLAG else "false"
    shell:
        """
        python {SCRIPTS_DIR}/process_rna_seq.py {wildcards.sample} {params.output_dir} {params.ncrna_flag} > {log} 2>&1
        """

rule calculate_rpm:
    input:
        rRNA_output_txt = expand(f"{RESULTS_DIR}/{GSE}/stats/{{sample}}_rRNA_output.txt", sample=SAMPLES),
        summary_txt = expand(f"{RESULTS_DIR}/{GSE}/stats/{{sample}}_summary.txt", sample=SAMPLES)
    output:
        rRNA_length_txt =f"{RESULTS_DIR}/{GSE}/stats/{{sample}}_rRNA_length_distribution.txt",
        rRNA_output_with_rpm = f"{RESULTS_DIR}/{GSE}/stats/{{sample}}_rRNA_output_with_rpm.txt"
    params:
        output_dir = f"{RESULTS_DIR}/{GSE}/stats"
    log:
        f"{RESULTS_DIR}/{GSE}/logs/{{sample}}_calculate_rpm.log"
    shell:
        """
        python {SCRIPTS_DIR}/RPM.py {GSE} {wildcards.sample} {params.output_dir} > {log} 2>&1
        """

rule count_matrix:
    input:
        rRNA_output_with_rpm = expand(f"{RESULTS_DIR}/{GSE}/stats/{{sample}}_rRNA_output_with_rpm.txt", sample=SAMPLES)
    output:
        rRNA_count_matrix = f"{RESULTS_DIR}/{GSE}/stats/{GSE}_rRNA_count_matrix.txt"
    params:
        output_dir = f"{RESULTS_DIR}/{GSE}/stats",
        samples = ",".join(SAMPLES),
        case_control = '"' + str(GROUPS).replace('"', '\\"') + '"'
    log:
        f"{RESULTS_DIR}/{GSE}/logs/{GSE}_calculate_CountMatrix.log"
    shell:
        """
        python {SCRIPTS_DIR}/CountMatrix.py {GSE} {params.output_dir} {params.samples} {params.case_control} > {log} 2>&1
        """

rule expr_matrix:
    input:
        rRNA_output_with_rpm = expand(f"{RESULTS_DIR}/{GSE}/stats/{{sample}}_rRNA_output_with_rpm.txt", sample=SAMPLES)
    output:
        rRNA_expression_matrix = f"{RESULTS_DIR}/{GSE}/stats/{GSE}_rRNA_expression_matrix.txt"
    params:
        output_dir = f"{RESULTS_DIR}/{GSE}/stats",
        samples = ",".join(SAMPLES),
        case_control = '"' + str(GROUPS).replace('"', '\\"') + '"' 
    log:
        f"{RESULTS_DIR}/{GSE}/logs/{GSE}_calculate_ExprMatrix.log"
    shell:
        """
        python {SCRIPTS_DIR}/ExprMatrix.py {GSE} {params.output_dir} {params.samples} {params.case_control} > {log} 2>&1
        """

rule plot_coverage:
    input:
        rRNA_expr_matrix = f"{RESULTS_DIR}/{GSE}/stats/{GSE}_rRNA_expression_matrix.txt"
    output:
        txts = expand(f"{RESULTS_DIR}/{GSE}/plots/{GSE}_{{rRNA}}_Coverage_Data.txt",rRNA=RRNA_LIST),
        plots = expand(f"{RESULTS_DIR}/{GSE}/plots/{GSE}_{{rRNA}}_Coverage.png",rRNA=RRNA_LIST),
        flag = temp(f"{RESULTS_DIR}/{GSE}/plots/{GSE}_coverage.finished")
    log:
        f"{RESULTS_DIR}/{GSE}/logs/{GSE}_plot_coverage.log"
    params:
        output_dir = f"{RESULTS_DIR}/{GSE}",
        samples = ','.join(SAMPLES),
        case_control = json.dumps(GROUPS).replace('"', '\\"')
    shell:
        """
        Rscript {SCRIPTS_DIR}/compare_Plot_rRNA_coverage.R -i {GSE} -o {params.output_dir} -s {params.samples} -g "{params.case_control}" > {log} 2>&1
        echo "Finished" > {output.flag}
        """

rule plot_length:
    input:
        rRNA_expr_matrix = f"{RESULTS_DIR}/{GSE}/stats/{GSE}_rRNA_expression_matrix.txt"
    output:
        length_txt = f"{RESULTS_DIR}/{GSE}/plots/{GSE}_rRNA_Length_distribution_RPM.txt",
        length_plot = f"{RESULTS_DIR}/{GSE}/plots/{GSE}_rRNA_Length_Distribution.png",
        subtype_txt = f"{RESULTS_DIR}/{GSE}/plots/{GSE}_rRNA_subtype_RPM.txt",
        subtype_plot = f"{RESULTS_DIR}/{GSE}/plots/{GSE}_rRNA_subtype_HollowPieChart.png",
        end_txt = f"{RESULTS_DIR}/{GSE}/plots/{GSE}_rRNA_End_RPM.txt",
        end_plot = f"{RESULTS_DIR}/{GSE}/plots/{GSE}_rRNA_End_StackBarPlot.png",
        flag = temp(f"{RESULTS_DIR}/{GSE}/plots/{GSE}_length.finished")
    log:
        f"{RESULTS_DIR}/{GSE}/logs/{GSE}_plot_length.log"
    params:
        output_dir = f"{RESULTS_DIR}/{GSE}",
        samples = ','.join(SAMPLES),
        case_control = json.dumps(GROUPS).replace('"', '\\"')
    shell:
        """
        Rscript {SCRIPTS_DIR}/rRNA_length_distribution.R \
            -i {GSE} \
            -o {params.output_dir} \
            -s {params.samples} \
            -g "{params.case_control}" >> {log} 2>&1
        echo "Finished" > {output.flag}
        """

rule plot_DA:
    input:
        matrix1 = f"{RESULTS_DIR}/{GSE}/stats/{GSE}_rRNA_count_matrix.txt",
        matrix2 = f"{RESULTS_DIR}/{GSE}/stats/{GSE}_rRNA_expression_matrix.txt",
        length_flag = f"{RESULTS_DIR}/{GSE}/plots/{GSE}_length.finished",
        coverage_flag = f"{RESULTS_DIR}/{GSE}/plots/{GSE}_coverage.finished"
    output:
        flag = temp(f"{RESULTS_DIR}/{GSE}/plots/{GSE}_DA.finished")
    log:
        f"{RESULTS_DIR}/{GSE}/logs/{GSE}_plot_DA.log"
    params:
        output_dir = f"{RESULTS_DIR}/{GSE}",
        samples = ','.join(SAMPLES),
        case_control = json.dumps(GROUPS).replace('"', '\\"'),
        n_case_samples = len(GROUPS["case"]),
        n_control_samples = len(GROUPS["control"])
    shell:
        r"""
        if [ {params.n_case_samples} -lt 3 ] && [ {params.n_control_samples} -lt 3 ]; then
            echo "[plot_DA] Not enough samples, skip DA analysis" >> {log}
            echo "SKIPPED due to insufficient samples" > {output.flag}
        else
            echo "[plot_DA] Running DA.R" >> {log}
            Rscript {SCRIPTS_DIR}/DA.R \
                -i {GSE} \
                -o {params.output_dir} \
                -s {params.samples} \
                -g "{params.case_control}" >> {log} 2>&1 || true

            echo "Finished" > {output.flag}
        fi
        """

rule export_result_report:
    input:
        alignment_summary=f"{RESULTS_DIR}/{GSE}/stats/{GSE}_alignment_summary.csv",
        matrix1 = f"{RESULTS_DIR}/{GSE}/stats/{GSE}_rRNA_count_matrix.txt",
        matrix2 = f"{RESULTS_DIR}/{GSE}/stats/{GSE}_rRNA_expression_matrix.txt",
        da_flag = f"{RESULTS_DIR}/{GSE}/plots/{GSE}_DA.finished"
    output:
        report = f"{RESULTS_DIR}/{GSE}/{GSE}_report.html"
    shell:
        """
        python {SCRIPTS_DIR}/generate_report.py {GSE} {CONFIGFILE} {RESULTS_DIR} {SCRIPTS_DIR} {output.report}
        """

rule merge_RNA_logs: 
    input: 
        report = f"{RESULTS_DIR}/{GSE}/{GSE}_report.html"
    output:
        flag = f"{RESULTS_DIR}/{GSE}/logs/.all.done"
    params:
        samples = ','.join(SAMPLES),
        ncrna_flag = "true" if NCRNA_FLAG else "false",
        files_flag = "keep" if KEEP_FILES_FLAG else "remove"
    shell:
        """
        bash {SCRIPTS_DIR}/re_logs.sh {GSE} {RESULTS_DIR} {params.samples} {params.ncrna_flag} {params.files_flag}
        """



if NCRNA_FLAG:
    # Custom step: Map reads to ncRNA
    # 如果ncrna_flag为true，则进行生成ncRNA mapping
    rule map_ncRNA:
        input:
            prev_unmapped = f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_MG_unmatch_45S_rRNA.fa"
        output: 
            stats = f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_MG_match_{{ncRNA}}_ncRNA.bowtie.stat",
            detail = f"{RESULTS_DIR}/{GSE}/mapped/{{sample}}_detail_match_{{ncRNA}}_ncRNA"
        threads: 8
        log:
            ncRNA_log = f"{RESULTS_DIR}/{GSE}/logs/{{sample}}_map_ncRNA_{{ncRNA}}.log"
        params:
            script_ncRNA = f"{SCRIPTS_DIR}/map_ncRNA.sh",
            output_dir = f"{RESULTS_DIR}/{GSE}",
            ncrna_flag = "true" if NCRNA_FLAG else "false"
        shell:
            """
                if [ "{params.ncrna_flag}" = "true" ]; then
                    bash {params.script_ncRNA} {wildcards.sample} {params.output_dir} {ANNOTATIONS_DIR} {threads} > {log.ncRNA_log} 2>&1
                    echo "Successfully completed mapping to ncRNA." >> {log.ncRNA_log}
                fi
            """