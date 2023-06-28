#https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_SCENIC-protocol-CLI.ipynb

# PACKAGES import dependencies
import json
import scanpy
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE

# VARIABLES
print("test")
with open("/SCRATCH-BIRD/users/nantunes/embryo-scenic/loom_creation.py_settings.json", "r") as f:
    settings = json.load(f)

# set a working directory
os.chdir(settings["WD"])

adata = sc.read_text(settings["exprDat"], delimiter="\t",first_column_names=True).transpose()

row_attrs = {
    "Gene": np.array(adata.var.index) ,
}
col_attrs = {
    "CellID": np.array(adata.obs.index) ,
    "nGene": np.array(np.sum(adata.X.transpose()>0 , axis=0)).flatten()
}
print(row_attrs)
print(col_attrs)
lp.create(settings["loom"], adata.X.transpose(), row_attrs, col_attrs)

# https://github.com/aertslab/pySCENIC/issues/82
