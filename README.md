# cdr3art

![cluster_0](https://user-images.githubusercontent.com/46639063/231309188-225288d0-8afd-4156-b74c-888916124148.png)

**cdr3art - concise visualization of T cell receptor clusters**

This R script allows a quick way to produce publication ready vector graphics with tidyverse and ggplot extensions, 
where it is often informative to see the gene usage pattern alongside CDR3 junction.

Critically, it allows for generation of multiple clusters of T cell receptors per page. 

Input files are produced by the Python package tcrdist3, but could be user supplied

## R Dependencies 

* dplyr
* ggseqlogo
* ggsankey
* ggplot2
* gridExtra
* purrr
* readr

## Commandline arguments 

```{verbatim}
cdr3art - concise visualization of T cell receptor clusters

flags:
  -h, --help            show this help message and exit

optional arguments:
  -x, --opts            RDS file containing argument values
  -i, --input           path to input .csv [default:
                        data/motif_graphic_instructions.csv]
  -o, --outdir          output folder [default: ./output]
  -c, --combine         combine graphics into a single file [default:
                        TRUE]
  -k, --k_per_page      number of graphics per page [default: 10]
  --individual          make individual graphics one per cluster
                        [default: FALSE]
  -w, --width           graphic width (inches) [default: 4]
  --height              graphic height (inches) [default: 2]
  -f, --font_size       font size of gene names [default: 4]
  -a, --axis_font_size  font size of logo axis text) [default: 6]
```

## Examples


### how can I visualize single chain motifs?

Here we generate pages with 5 clusters per page, with width of 4 inches and height 2 inches per diagram.

```{verbatim}
mkdir beta_chain_example
Rscript cdr3art.R --input data/motif_graphic_instructions.csv \
  --outdir beta_chain_example \
  --individual FALSE \
  --combine TRUE \
  --k_per_page 5 \
  --width 4 \
  --height 2 
```

### How can I visualize paired chain motifs?

Use `--combine' set to TRUE, to generate pages with multiple diagrams per page. 

```{verbatim}
mkdir alpha_beta_chain_example
Rscript cdr3art.R --input data/motif_graphic_instructions_alpha_beta_example.csv \
  --outdir alpha_beta_chain_example \
  --individual FALSE \
  --combine TRUE \
  --k_per_page 8 \
  --width 8 \
  --height 2.0 \
  --font_size 1
```


### How can I make one motif per file?

Use `--individual' set to TRUE, to generate a graphic for each cluster_id in the input.

```{verbatim}
mkdir alpha_beta_chain_example
Rscript cdr3art.R --input data/motif_graphic_instructions_alpha_beta_example.csv \
  --outdir alpha_beta_chain_example \
  --individual TRUE \
  --combine FALSE \
  --k_per_page 8 \
  --width 8.5 \
  --height 1.5 \
  --font_size 1.25
```


```{verbatim}
mkdir just_one_cluster
Rscript cdr3art.R --input data/cluster_0.csv \
  --outdir just_one_cluster \
  --individual TRUE \
  --combine FALSE \
  --k_per_page 8 \
  --width 11 \
  --height 3 \
  --font_size 1.25
```

### TODO

* Add functionality to let the user specify the order of the clusters in the motif folio. 
* Allow user to pass custom cluster names with a 'cluster_name' column.


## How do I create the input files with tcrdist3?

Creating background subtracted motifs is a subject in it's own right, 
but here is some template code that can complete this task. See tcrdist3 doc for more details. 

Here is an example where each connected component within a TCRb sequence 
similarity network is assigned a cluster_id.


```bash
#!python -c "from tcrsampler.setup_db import install_all_next_gen; install_all_next_gen(dry_run = False)"
```

```python
import pandas as pd
from tcrdist.repertoire import TCRrep
from tcrdist.public import _neighbors_fixed_radius
import pandas as pd
import os
import networkx as nx
from tcrdist.html_colors import get_html_colors
from matplotlib import pyplot as plt
from IPython.core.display import display, HTML
from palmotif import compute_pal_motif, svg_logo
from tcrsampler.sampler import TCRsampler   
from tcrdist.regex import  _matrix_to_regex

df = pd.read_csv(YOUR_CORRECTLY_FORMATED_INPUT_DATA)

tr = TCRrep(cell_df = df 
            organism = 'human', 
            chains = ['beta'])
# Compute Network
# <tr.pw_alpha_beta> is paired chain TCRdist.
edge_threshold =12 # edit distance of 1 - 1.5, equivailent
# <network> initialize a list to populate with edges between TCRs.
network = list()
for i,n in enumerate(_neighbors_fixed_radius(tr.pw_beta, edge_threshold)):
    for j in n:
        if i != j:
          row = list()
          row.append(i)
          row.append(j)
          row.append(len(n)-1)
          row.append((tr.pw_beta )[i,j])
          for col in tr.clone_df.columns:
            row.append(tr.clone_df[col].iloc[i])
            row.append(tr.clone_df[col].iloc[j])
          network.append(row)                      # 'K_neighbors' - number of neighbors
import numpy as np
cols = ['node_1', 'node_2','K_neighbors','dist'] + np.concatenate([(f"{x}_1", f"{x}_2") for x in tr.clone_df.columns.to_list()]).tolist()
df_net = pd.DataFrame(network, columns = cols)
df_net['weight'] = (edge_threshold - df_net['dist'])/edge_threshold

G = nx.from_pandas_edgelist(pd.DataFrame({'source' : df_net['node_1'],
                                          'target' : df_net['node_2'],
                                          'weight' : df_net['weight'],
                                          'count1' : df_net['count_1'],
                                          'count2' : df_net['count_2']}))

cc_to_node = dict()
node_to_cc = dict()
for i,x in enumerate(sorted(nx.connected_components(G), key=len, reverse=True)):
  cc_to_node[i] = x
  for j in x:
    node_to_cc[j] = i
tr.clone_df['i'] = range(0, tr.clone_df.shape[0])
tr.clone_df['graph_cc'] = tr.clone_df['i'].apply(lambda i : node_to_cc.get(i))

tr.clone_df['cluster_id_cat'] =  pd.Series(tr.clone_df['graph_cc'].to_list(), dtype="category").\
  cat.reorder_categories(tr.clone_df['graph_cc'].value_counts().index)
  
  
  
# This step can take up to 1 minute, so be patient
ts_beta = TCRsampler(default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv')
ts_beta.build_background(max_rows = 50, stratify_by_subject = True)
ts_alpha = TCRsampler(default_background = 'olga_human_alpha_t.sampler.tsv')
ts_alpha.build_background(max_rows = 50, stratify_by_subject = True)

def custom_logo_routine(cluster_df, label,  alpha = False, beta = True,):
  if alpha:
    # Identify the gene usage pattern for the alpha and beta chains
    gene_usage_alpha = cluster_df.groupby(['v_a_gene','j_a_gene']).\
        size().\
        reset_index().\
        to_dict('split')['data']  
    # Sample a background for the alpha chain
    sampled_rep_alpha = ts_alpha.sample(gene_usage_alpha, flatten = True, depth = 100)
    # Remove any None values that could have been generated by unknown pairs
    sampled_rep_alpha =  [x for x in sampled_rep_alpha if x is not None]
    # Compute background subtracted positionally aligned motif logo
    motif, stat = compute_pal_motif(
        seqs = cluster_df['cdr3_a_aa'].to_list(),
        refs = sampled_rep_alpha, 
        centroid = cluster_df['cdr3_a_aa'].value_counts().index[0])
    # Compute raw positionally aligned motif logo
    background_subtracted_svg_alpha = svg_logo(motif, return_str= True).\
      replace('height="100%"', 'height="10%"').\
      replace('width="100%"', 'width="10%"')
    # Compute raw positionally aligned motif logo
    motif_raw, stat_raw = compute_pal_motif(
        seqs = cluster_df['cdr3_a_aa'].to_list(),
        refs = None,
        centroid = cluster_df['cdr3_a_aa'].value_counts().index[0])
      
    mr_a = motif_raw.copy()
    m_a = motif.copy()
    raw_svg_alpha = svg_logo(motif_raw, return_str= True).\
      replace('height="100%"', 'height="10%"').\
      replace('width="100%"', 'width="10%"')
    # Generate a regeular expression from the pal motif matrix
    regex_alpha =  _matrix_to_regex(motif_raw,max_ambiguity=5, ntrim=0, ctrim=0)
  else:
      mr_a = None
      m_a = None

  if beta:
    gene_usage_beta = cluster_df.groupby(['v_b_gene','j_b_gene']).\
      size().\
      reset_index().\
      to_dict('split')['data']

    # Reapeat for the the BETA chain 
    # Sample a background for the beta chain
    sampled_rep_beta = ts_beta.sample(gene_usage_beta, flatten = True, depth = 100)
    # remove any None values that could have been generated by unknown pairs
    sampled_rep_beta =  [x for x in sampled_rep_beta if x is not None]

    # Compute background subtracted positionally aligned motif logo
    motif, stat = compute_pal_motif(
        seqs = cluster_df['cdr3_b_aa'].to_list(),
        refs = sampled_rep_beta, 
        centroid = cluster_df['cdr3_b_aa'].value_counts().index[0])
    background_subtracted_svg_beta = svg_logo(motif, return_str= True).\
        replace('height="100%"', 'height="10%"').\
        replace('width="100%"', 'width="10%"')
    # Compute raw positionally aligned motif logo
    motif_raw, stat_raw = compute_pal_motif(
        seqs = cluster_df['cdr3_b_aa'].to_list(),
        refs = None,
        centroid = cluster_df['cdr3_b_aa'].value_counts().index[0])
    raw_svg_beta = svg_logo(motif_raw, return_str= True).\
        replace('height="100%"', 'height="10%"').\
        replace('width="100%"', 'width="10%"')
    regex_beta =  _matrix_to_regex(motif_raw, max_ambiguity=5, ntrim=0, ctrim=0)
    mr_b = motif_raw.copy()
    m_b = motif.copy()
  else:
    mr_b = None
    m_b = None
  return(mr_a,m_a, mr_b, m_b)
external_frames = list()
for i,cluster_df in tr.clone_df.groupby('cluster_id_cat'):
  if cluster_df.shape[0] > 2:
    
    cluster_id = cluster_df['cluster_id_cat'].iloc[0]
    mr_a ,m_a, mr_b, m_b= custom_logo_routine(cluster_df = cluster_df, label = f"Cluster {cluster_id}", alpha = False) 

    mr_b.columns = [f"x{i}" for i in range(len(mr_b.columns))]
    m_b.columns  = [f"x{i}" for i in range(len(m_b.columns))]
    m_b['type']  = 'subtracted'
    mr_b['type'] = 'raw'
    mr_b['chain'] = "b"
    m_b['chain']  = "b"
    #gu = cluster_df[['v_a_gene','j_a_gene','v_b_gene','j_b_gene']].copy()
    gu = cluster_df[['v_b_gene','j_b_gene']].copy()
    gu['type'] = 'gene_usage'
    gu['chain']  = "b"
    external = pd.concat([mr_b,m_b,gu]).reset_index()
    external['cluster_id'] = cluster_id
    #external.to_csv(f'transfers/external_test_{cluster_id}.csv')
    external_frames.append(external)

pd.concat(external_frames).to_csv("motif_graphic_instructions.csv", index = False)
````








