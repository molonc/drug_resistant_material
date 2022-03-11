# library(devtools)
# devtools::install(pkg = "/home/htran/storage/install_software/scshaper/")
# devtools::install(pkg = "/home/htran/storage/install_software/dyno/")
# Load scShaper to use the prosstt1 dataset, which comes with the package
# BiocManager::install("dynwrap")
# BiocManager::install("dyneval")
# BiocManager::install("dyntoy")
# BiocManager::install("processx")
library(processx)
suppressPackageStartupMessages(library(scShaper))
suppressPackageStartupMessages(library(dynwrap))
suppressPackageStartupMessages(library(dyneval))
suppressPackageStartupMessages(library(tidyverse))
library(dyndimred)

data("prosstt1")

dataset <- wrap_expression(
  counts = prosstt1$rna.raw,
  expression = prosstt1$rna.normalized
)

set.seed(123456)
model <- infer_trajectory(dataset, ti_scShaper(k.range = 2:100,
                                               num.pcs = 50,
                                               span = 0.1),verbose = TRUE)

suppressPackageStartupMessages(library(dyno))
# Generate MDS dimensionality reduction for visualization (MDS is default of dyno)
model <- model %>% add_dimred(dyndimred::dimred_mds, expression_source = dataset$expression)
# Visualize the trajectory
?plot_dimred
head(get_dimred(model), 5)
plot_dimred(model,color_cells ='milestone')
# The importance for the whole trajectory
plot_dendro(model)
plot_dendro(model, "pseudotime")
overall_feature_importance <- calculate_overall_feature_importance(model, expression_source=dataset$expression)
View(head(overall_feature_importance))
# Take top 100 genes
prosstt1_features <- overall_feature_importance %>% 
  top_n(100, importance) %>% 
  pull(feature_id)

plot_heatmap(
  model, 
  expression_source = dataset$expression, 
  features_oi = prosstt1_features
)
plot_heatmap(model, expression_source = dataset)


definition <- definition(
  method = def_method(
    id = "scShaper"
  ),
  parameters = def_parameters(
    dynparam::integer_parameter(
      id = "k.range",
      default = 2:100,
      distribution = dynparam::uniform_distribution(2, 1000),
      description = "The range of k values used in estimating discrete pseudotimes"
    ),
    dynparam::numeric_parameter(
      id = "span",
      default = 0.1,
      distribution = dynparam::uniform_distribution(0, 1),
      description = "The parameter that determines the degree of smoothing in LOESS"
    ),
    dynparam::integer_parameter(
      id = "num.pcs",
      default = 50,
      distribution = dynparam::uniform_distribution(4, 100),
      description = "Before t-SNE how many principal components to use."
    ),
    dynparam::numeric_parameter(
      id = "perplexity",
      default = 30,
      distribution = dynparam::uniform_distribution(2, 100),
      description = "Perplexity of t-SNE"
    )
  ),
  wrapper = def_wrapper(
    input_required = "expression",
    input_optional = "start_id"
  )
)


run_fun <- function(expression, priors, parameters, seed, verbose) {
  
  suppressMessages(sce <- SingleCellExperiment(assays = list(counts = t(expression),logcounts = t(expression))))
  sce <- RunscShaper(sce,k.range = parameters$k.range,
                     span=parameters$span,
                     num.pcs = parameters$num.pcs,
                     tsne.perplexity = parameters$perplexity)
  pseudotime <- sce@metadata$scshaper$continuous.pseudotime
  
  names(pseudotime) <- rownames(expression)
  
  # flip pseudotimes using start_id
  if (!is.null(priors$start_id)) {
    if(mean(pseudotime[start_id]) > 0.5) {
      pseudotime <- 1-pseudotime
    }
  }
  
  dynwrap::wrap_data(cell_ids = rownames(expression)) %>%
    dynwrap::add_linear_trajectory(pseudotime = pseudotime)
}

ti_scShaper <- create_ti_method_r(definition, run_fun, package_loaded = c("SingleCellExperiment","dplyr","scShaper"))



