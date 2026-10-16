

###################################################################################
# clustering schema - cluster sequences using seurat_dataset, then match back to wide versions of donor-id sequences

import pandas as pd
import numpy as np
import os
import networkx as nx
from importlib.resources import files
from tcrdistgpu.distance import TCRgpu

INPUT_DIR = #"" #edit as needed
OUTPUT_DIR = #"" #edit as needed

#donors 1-5 
df = pd.read_csv(os.path.join(INPUT_DIR, 'donors_1-5_wide_version_cmv_121825.csv'))
#or - donor 0 
df = pd.read_csv(os.path.join(INPUT_DIR, 'donor_0_wide_version_cmv_121825.csv'))


#check for maits & inkts and drop 
x = list()
maitv = df['v_a_gene'] == "TRAV1-2*01" #1124 
j_a_genes_to_check = ["TRAJ33*01", "TRAJ20*01", "TRAJ12*01"]
maitj = df['j_a_gene'].isin(j_a_genes_to_check) 

for i,j in zip(maitv, maitj):
    if i and j:
        x.append("True")
    else:
        x.append("False")
df['mait'] = x
df['mait'].value_counts() ##223 MAITs 

x = list()
inktv = df['v_a_gene'] == "TRAV10*01" 
inktvb = df['v_b_gene'] == "TRBV25-1*01" 
inktj = df['j_a_gene'] == "TRAJ18*01"

for i,j,k in zip(inktv, inktvb, inktj):
  if i and j and k:
    x.append("True")
  else:
    x.append("False")
df['inkt'] = x
df['inkt'].value_counts() ##3 iNKTs

print(pd.crosstab(df['mait'], df['inkt']))
df['mait'] = df['mait'].map({'False': False, 'True': True})
df['inkt'] = df['inkt'].map({'False': False, 'True': True})
df = df[(df['mait'] == False) & (df['inkt'] == False)]
df #9334, 19
df = df.reset_index(drop = True)

df_tcr = df[['cdr3_a_aa', 'cdr3_b_aa', 'v_a_gene', 'v_b_gene', 'j_a_gene',
       'j_b_gene', 'donor']]
df_tcr.drop_duplicates() #no change
df_tcr.index #fine 

chain = "ab"
max_dist = 64

tg = TCRgpu(mode = 'cpu',# If you only have a cpu set mode to "cpu"
  chunk_size = 100,  # Limit rows processed at once to prevent RAM overflow
  cdr3a_col = 'cdr3_a_aa',
  cdr3b_col = 'cdr3_b_aa',
  va_col = 'v_a_gene',
  vb_col = 'v_b_gene')
encoding = tg.encode_tcrs(df_tcr) # check v and j genes if this errors 
csr_dists_max_dist = tg.compute_csr_experimental(encoding, encoding, max_dist = max_dist)

G = nx.from_scipy_sparse_array(csr_dists_max_dist)

## add in cluster labels and save csv - 
components = list(nx.connected_components(G))
components_sorted = sorted(components, key=len, reverse=True)

# assign cluster IDs by size
cluster_labels = {}
for cid, comp in enumerate(components_sorted):
    for node in comp:
        cluster_labels[node] = cid

# write back into data; -1 if node not in any component (e.g. singletons you kept)
df["cluster"] = [cluster_labels.get(i, "") for i in range(len(df))]
print(f"Found {len(components)} clusters")

df.to_csv(os.path.join(OUTPUT_DIR, 'cmv_d1-5_dataset_paired_clusters_et64_042826.csv'), index = False)


## graph with cluster labels - 
import matplotlib.pyplot as plt
import networkx as nx
from tcrdist.html_colors import get_html_colors
from tcrdistgpu.draw import draw_categories

data = df.reset_index(drop=True)

# options to graph by different stim conditions 
data['count'] = df["IE1_Late"] + df['IE1_Early'] 
data['count'] = df["pp65_Late"] + df['pp65_Early']
data['count'] = df["Bulk_Late"] + df['Bulk_Early']

# mod for graphing 
chain = "ab"
max_dist = 64

tg = TCRgpu(mode = 'cpu',# If you only have a cpu set mode to "cpu"
  chunk_size = 100,  # Limit rows processed at once to prevent RAM overflow
  cdr3a_col = 'cdr3_a_aa',
  cdr3b_col = 'cdr3_b_aa',
  va_col = 'v_a_gene',
  vb_col = 'v_b_gene')
encoding = tg.encode_tcrs(data) # check v and j genes if this errors 
csr_dists_max_dist = tg.compute_csr_experimental(encoding, encoding, max_dist = max_dist)

G = nx.from_scipy_sparse_array(csr_dists_max_dist)

## add in cluster labels and save csv - 
components = list(nx.connected_components(G))
components_sorted = sorted(components, key=len, reverse=True)

# assign cluster IDs by size
cluster_labels = {}
for cid, comp in enumerate(components_sorted):
    for node in comp:
        cluster_labels[node] = cid

# write back into data; -1 if node not in any component (e.g. singletons you kept)
df["cluster"] = [cluster_labels.get(i, "") for i in range(len(df))]
print(f"Found {len(components)} clusters")

G.remove_edges_from(nx.selfloop_edges(G)) ### remove self-loop edges 

# drop single nodes (just dots, not barbells)
G_filtered = G.copy()
for component in list(nx.connected_components(G_filtered)):
  if len(component) < 2: # can remove singletons or small connected compoents
    G_filtered.remove_nodes_from(component)
print("Re-filtered graph:", G_filtered)

### use mod above for graphing to generate G_filtered (just for visualization)
pos = nx.spring_layout(G_filtered, seed=1)  

# Layout 
pos = nx.nx_agraph.graphviz_layout(G_filtered, prog="neato")

plt.figure(figsize=(8, 8))

#  choose what to color by
color_var = "donor"   # or "cluster" or "cell_type"

# build node property dicts directly from `data`
nodes = list(G_filtered.nodes())
node_cts = {n: max(0, int(data.loc[n, 'count'])) for n in nodes} ### min vs max if expanded cells
color_vals = pd.Series([data.loc[n, color_var] for n in nodes])

# map category -> color
colors = get_html_colors(color_vals.nunique())
var_to_color = {v: c for v, c in zip(color_vals.value_counts().index, colors)}
#var_to_color = {'IE1': '#1f77b4', "not": '#aec7e8', 'pp65': '#ff7f0e','cmv_bulk': '#ffbb78'}
node_to_color_var = {n: data.loc[n, color_var] for n in nodes}

# connected components of the *filtered* graph
# already called above 

for comp in components:
    g = G_filtered.subgraph(comp)
    comp_nodes = list(g.nodes())

    # pick node within this component
    first_node = max(comp_nodes, key=lambda n: g.degree[n])
    comp_cluster = int(data.loc[first_node, "cluster"])

    my_labels = {n: "" for n in comp_nodes}
    my_labels[first_node] = str(comp_cluster)

    nx.draw(g, pos,
        edge_color="grey",
        alpha=0.7,
        node_color=[var_to_color[node_to_color_var[n]] for n in g.nodes()],
        node_size=[node_cts[n] for n in g.nodes()],
        vmin=0.0, vmax=1.0,
        with_labels=True,
        labels=my_labels,
        font_color="black",
        font_size=10)

plt.axis("off")
plt.tight_layout()
plt.savefig("cmv_d1-5_paired_clusters_bulk_et64_042826.png", dpi=300)
plt.clf()
