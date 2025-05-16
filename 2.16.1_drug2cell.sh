import sys
import scanpy as sc
import drug2cell as d2c
import blitzgsea as blitz
import pandas as pd
import anndata as ad
import pooch
import numpy as np
import os
sc.settings.set_figure_params(dpi=800)

adata = sc.read('./files/PH_raw.h5ad')
adata.raw = sc.read('./files/PH_data.h5ad')
sc.pp.scale(adata)
d2c.score(adata, use_raw=True)
adata.uns['drug2cell']

##
d2c.score(adata, use_raw=True)
adata.uns['drug2cell']

sc.tl.rank_genes_groups(adata.uns['drug2cell'], method="wilcoxon", groupby="cell_type")
sc.pl.rank_genes_groups_dotplot(adata.uns['drug2cell'], swap_axes=True, dendrogram=False, n_genes=10, save=os.path.join(SAVE_PATH,'drug2cell_dotplat.pdf'))

plot_args = d2c.util.prepare_plot_args(adata.uns['drug2cell'], categories=["B01","B02","B03"])
sc.pl.dotplot(adata.uns['drug2cell'], groupby="cell_type", swap_axes=True, **plot_args, save=os.path.join(SAVE_PATH,'drug2cell_dotplat2.pdf'))

## save
adata.write('/home/chencx/figures/results/PH_drug.h5ad')
pd.DataFrame(adata.uns['drug2cell'].uns['rank_genes_groups']['names']).to_csv('./files/drug_names.csv', index=False)
pd.DataFrame(adata.uns['drug2cell'].uns['rank_genes_groups']['scores']).to_csv('./files/drug_scores.csv', index=False)
pd.DataFrame(adata.uns['drug2cell'].uns['rank_genes_groups']['pvals_adj']).to_csv('./files/drug_pvals_adj.csv', index=False)
pd.DataFrame(adata.uns['drug2cell'].uns['rank_genes_groups']['pvals']).to_csv('./files/drug_pvals.csv', index=False)
pd.DataFrame(adata.uns['drug2cell'].uns['rank_genes_groups']['logfoldchanges']).to_csv('./files/drug_logFC.csv', index=False)

#### target genes ----
import pickle
with open('./files/drug-target_dicts.pkl', 'rb') as file:
    drug_target_dicts = pickle.load(file)

key_count = len(drug_target_dicts)
print(f"There are {key_count} keys")
keys = list(drug_target_dicts.keys())
print(f"They are：{keys}")

## 
data = []
for drug_class, compounds in drug_target_dicts.items():
    for compound, targets in compounds.items():
        for gene in targets:
            data.append({'Gene': gene, 'Drug': compound, 'Category': drug_class})

df = pd.DataFrame(data)
df.to_csv('./files/drug-target_dicts.csv', index=False)

