# PACKAGES
import snakemake, os
import numpy as np

# CONFIG
configfile:"config/config.json"
SOFTPOWERS = np.arange(start=int(config["MIN_SOFTPOWER"]), stop=int(config["MAX_SOFTPOWER"]), step=1, dtype=int).tolist()

# RULES
rule all:
    input:
        expand(os.path.join(config["OUTPUT_FOLDER"] + "_{softpower}", "_no_grey_heatmap.pdf"), softpower = SOFTPOWERS),

rule heatmap:
    input:
        markers_path = os.path.join(config["WD"], "WGCNA/results/genesModulesCorrespondence_SF-{softpower}_no_grey.tsv"),
        annotation_path = config["ANNOTATION"],
        expr_Mat_path = config["EXPR"],
        heatmap_config_path = os.path.join(os.getcwd(), "config/config_heatmap.json")
    output:
        file = os.path.join(config["OUTPUT_FOLDER"] + "_{softpower}", "_no_grey_heatmap.pdf")
    params:
        output_folder = config["OUTPUT_FOLDER"] + "_{softpower}"
    conda:
        "envs/heatmap.yaml"
    script:
        config["SCRIPT_HEATMAP"]

