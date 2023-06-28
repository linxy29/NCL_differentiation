# PACKAGES
import snakemake, os
import numpy as np

# CONFIG
configfile:"config/config.json"
SOFTPOWERS = np.arange(start=int(config["MIN_SOFTPOWER"]), stop=int(config["MAX_SOFTPOWER"]), step=1, dtype=int).tolist()

# RULES
rule all:
    input:
        expand(os.path.join(config["OUTPUT_FOLDER"] + "_{softpower}", "heatmap.pdf"), softpower = SOFTPOWERS),
        os.path.join(config["WD"], "softpower.pdf")

rule heatmap:
    input:
        markers_path = os.path.join(config["WD"], "WGCNA/results/genesModulesCorrespondence_SF-{softpower}.tsv"),
        annotation_path = config["ANNOTATION"],
        expr_Mat_path = config["EXPR"],
        heatmap_config_path = os.path.join(os.getcwd(), "config/config_heatmap.json")
    output:
        file = os.path.join(config["OUTPUT_FOLDER"] + "_{softpower}", "heatmap.pdf")
    params:
        output_folder = config["OUTPUT_FOLDER"] + "_{softpower}"
    conda:
        "envs/heatmap.yaml"
    script:
        config["SCRIPT_HEATMAP"]

rule wgcna:
    input:
        WD = config["WD"],
        ANNOTATION = config["ANNOTATION"],
        EXPR = config["EXPR"]
    output:
        output_file = os.path.join(config["WD"], "WGCNA/results/genesModulesCorrespondence_SF-{softpower}.tsv"),
    params:
        MIN_SOFTPOWER = config["MIN_SOFTPOWER"],
        MAX_SOFTPOWER = config["MAX_SOFTPOWER"],
        MAX_BLOCKSIZE = config["MAX_BLOCKSIZE"],
        SOFTPOWER = "{softpower}"
    conda:
        "envs/wgcna.yaml"
    script:
        config["SCRIPT_WGCNA"]

rule softpower:
    input:
        WD = config["WD"],
        EXPR = config["EXPR"]
    output:
        os.path.join(config["WD"], "softpower.pdf")
    params:
        MIN_SOFTPOWER = config["MIN_SOFTPOWER"],
        MAX_SOFTPOWER = config["MAX_SOFTPOWER"],
        MAX_BLOCKSIZE = config["MAX_BLOCKSIZE"],
    conda:
        "envs/wgcna.yaml"
    script:
        config["SCRIPT_SOFTPOWER"]
