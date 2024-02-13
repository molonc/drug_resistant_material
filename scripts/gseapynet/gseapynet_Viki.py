## improved by Viki on 15 Dec 2020
import pandas
import logging
import gseapy as gp
import networkx as nx
import collections
import math
from itertools import combinations
import matplotlib.pyplot as plt
import matplotlib
import pygraphviz as pgv

import matplotlib.gridspec as gridspec
import numpy
import sys
from gseapy.plot import barplot
import argparse

### Change color map
cmap = matplotlib.cm.coolwarm



parser = argparse.ArgumentParser(description='Run ranked GSEA with network output.')
parser.add_argument('--csv', type=str, help='CSV with DE genes.')
parser.add_argument('--sig', type=float, help='Threshold for signficance in DE input.')
parser.add_argument('--sig-col', type=str, help='Name of column with signficance in input csv.')
parser.add_argument('--minfc', type=float, help='Threshold for fold change column in csv input.')
parser.add_argument('--minfc-col', type=str, help='Name of column with fold change in input csv.')
parser.add_argument('--pathway-fdr', type=float, help='FDR threshold for GSEA pathway result.')
parser.add_argument('--gene-col', type=str, help='Name of column for gene symbol or Id.')
parser.add_argument('--png', type=str, help='Name of png output file. (end in .svg for svg result).')
parser.add_argument('--gmt', type=str, help='Either a custom gmt file or the name of an enrichr library (https://amp.pharm.mssm.edu/Enrichr/#stats)')
# MA: 18 Dec 2020, adding a title
parser.add_argument('--title', type=str, help='Title for the plot')
parser.add_argument('--figw', type=str, help='Figure width in inches')
parser.add_argument('--figh', type=str, help='Figure height in inches')
args = parser.parse_args()


plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams["font.family"] = "Helvetica"
plt.rc('font', size=6)          # controls default text sizes
plt.rc('axes', titlesize=7)     # fontsize of the axes title
plt.rc('axes', labelsize=7)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=6)    # fontsize of the tick labels
plt.rc('ytick', labelsize=6)    # fontsize of the tick labels
plt.rc('legend', fontsize=7)    # legend fontsize
plt.rc('figure', titlesize=10)  # fontsize of the figure title

logger = logging.getLogger("Pathway Analysis")
logger.setLevel(logging.INFO)
logging.basicConfig(level=logging.INFO)

csv = args.csv
fdr = args.pathway_fdr
minfc = args.minfc
adjp = args.minfc_col
minedge = 1
png = args.png

logger.info("Reading CSV file")
deg = pandas.read_csv(csv)

deg = deg[abs(deg[args.sig_col]) < float(args.sig)]
#fig = plt.figure(figsize=(10,7))
print("Figure sizes")
print(args.figw)
print(args.figh)
print("fdr")
print(fdr)
fig = plt.figure(figsize=(float(args.figw),float(args.figh)))
gs = gridspec.GridSpec(nrows=2, ncols=2, width_ratios=[30,1], height_ratios=[30, 1])
axes = []
scales = []
axes.append(fig.add_subplot(gs[0, 0]))
axes.append(fig.add_subplot(gs[1, 0]))
axes.append(fig.add_subplot(gs[0, 1]))
margin=0.05
fig.subplots_adjust(margin, margin, 1.-margin, 1.-margin)

deg = deg[abs(deg[args.minfc_col]) > float(minfc)]
data = deg[[args.gene_col,args.minfc_col]]

data = data.sort_values(args.minfc_col)

pre_res = gp.prerank(rnk=data, gene_sets=args.gmt,
                    processes=4,
                    permutation_num=100,
                    outdir="prerank", format='png', seed=6,min_size=1,max_size=30000000000)

print("Printing head of pre_res")
print(pre_res.res2d.sort_index().head())
result = pre_res.res2d
result = result[result["fdr"]<float(args.pathway_fdr)] #### ADD ME BACK
#
result.to_csv(png.replace(".png",".csv").replace(".svg",".csv").replace(".pdf",".csv"),header=True)
pathways = result.index
nodes = collections.defaultdict(dict)
edges = dict()


minnes = min(result["nes"])
maxnes = max(result["nes"])

print("Pathways")
print(pathways)

for term in pathways:
    hallmark = term
    print("Genes")
    print(result.loc[hallmark]["genes"])
    term = term.replace("HALLMARK_","").replace("_"," ")
    nodes[term]["genes"] = result.loc[hallmark]["genes"].split(";")
    nodes[term]["nes"] = result.loc[hallmark]["nes"]
G = nx.petersen_graph()
node_color = []
node_order = []
node_size = []
for term in nodes:
    node_order.append(term)
    node_color.append(float(nodes[term]["nes"]))
    node_size.append(abs(float(nodes[term]["nes"] * 50.0)))
    G.add_node(term, nes=float(nodes[term]["nes"]))

if len(node_color) == 1:
    print("No Enriched Pathways! Try relaxing the parameters.")
    exit(0)

edge_order = []
edge_color = []
overlaps = []
for edge in combinations(node_order,2):
    overlap = set(nodes[edge[0]]["genes"]).intersection(nodes[edge[1]]["genes"])
    if len(overlap) > 0:
        overlaps.append(len(overlap))
        G.add_edge(edge[0], edge[1], weight=len(overlap))
        print(edge[0],edge[1],len(overlap), overlap)
        edge_color.append(len(overlap))
        edge_order.append(edge)

G.remove_nodes_from(list(nx.isolates(G)))
if len(overlaps) == 0:
    maxvedge = 1
else:
    maxvedge = max(overlaps)
_node_order = []
_node_size = []
_node_color = []
for nord, nsiz, ncol in zip(node_order, node_size, node_color):
    if nord in G.nodes():
        _node_order.append(nord)
        _node_size.append(nsiz)
        _node_color.append(ncol)

node_order = _node_order
node_size = _node_size
node_color = _node_color
print(len(node_order), len(node_size), len(G.nodes()))
#k=1/math.sqrt(len(G.nodes()))
pos = nx.nx_agraph.graphviz_layout(G, prog="neato")
# MA: 18 Dec 2020: using a title given as argument
#axes[0].set_title(args.png.replace(".png","").replace(".svg",""), fontsize=11)
axes[0].set_title(args.title, fontsize=10, fontfamily="Helvetica")
axes[0].axis('equal')



nx.draw(G,pos,ax=axes[0], nodelist=node_order, node_size=node_size, vmin=minnes, vmax=maxnes,edgelist=edge_order, node_color=node_color, edge_color=edge_color, cmap=cmap, edge_cmap=plt.cm.Greys, edge_vmin=0, edge_vmax=maxvedge, with_labels=True, font_size=6, font_family="Helvetica")

# MA 2 Oct 2020: color bar doesn't show the right colors
#norm = matplotlib.colors.Normalize(vmin=-8, vmax=18)
norm = matplotlib.colors.Normalize(vmin=minnes, vmax=maxnes)
cb1 = matplotlib.colorbar.ColorbarBase(axes[1], cmap=cmap,
                                        norm=norm,
                                        orientation='horizontal')
axes[1].set_title("Normalized Enrichment Score",fontsize=7, fontfamily="Helvetica")
if overlaps == []:
    overlaps.append(12)

cmap2 = matplotlib.cm.binary
norm2 = matplotlib.colors.Normalize(vmin=0, vmax=max(overlaps))
cb2 = matplotlib.colorbar.ColorbarBase(axes[2], cmap=cmap2,
                                        norm=norm2,
                                        orientation='vertical', drawedges=False)
cb2.set_ticks(list(range(0,maxvedge,1)))
axes[2].set_title("Common Genes",fontsize=7, fontfamily="Helvetica")
plt.tight_layout()
#plt.constrained_layout()
#plt.savefig(png,dpi=1200)
## savig it as pdf
plt.savefig(png)
