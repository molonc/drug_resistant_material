## Sep 27 2020, maling feather files for SA609 new untreated line
# 4 Oct 2020, rerunning all the feather files to make sure I have the correct versions

## SA906b
for filename in ../../data-SA906b/dlp/csv/*.csv; do
  sample=$(basename "$filename" .csv)
  echo "Sample $sample"
  outfile="../../data-SA906b/dlp/$sample.feather"

  if [ -f "$outfile" ]; then
    echo "$outfile exists"
  else
    Rscript ../parse_cnv/R/format_segments.R \
      --segment_file $filename \
      --outfname $outfile \
      >& ../../data-SA906b/dlp/csv/$sample.log
  fi    
done   

## SA906a
for filename in ../../data-SA906a/dlp/csv_x10_x15/*.csv; do
  sample=$(basename "$filename" .csv)
  echo "Sample $sample"
  outfile="../../data-SA906a/dlp/$sample.feather"

  if [ -f "$outfile" ]; then
    echo "$outfile exists"
  else
    Rscript ../parse_cnv/R/format_segments.R \
      --segment_file $filename \
      --outfname $outfile \
      >& ../../data-SA906a/dlp/csv_x10_x15/$sample.log
  fi    
done   

# SA609
# Run format_segments to find gene by cell in the dlp data
for filename in ../../data-SA609/dlp/csv-v6/*.csv; do
  sample=$(basename "$filename" .csv)
  echo "Sample $sample"
  outfile="../../data-SA609/dlp/csv-v6/$sample.feather"

  if [ -f "$outfile" ]; then
    echo "$outfile exists"
  else
    Rscript ../parse_cnv/R/format_segments.R \
      --segment_file $filename \
      --outfname $outfile \
      >& ../../data-SA609/dlp/csv-v6/$sample.log
  fi    
done    

### SA535
for filename in ../../data-SA535/dlp-cisplatin_v7/csv/*.csv; do
  sample=$(basename "$filename" .csv)
  echo "Sample $sample"
  outfile="../../data-SA535/dlp-cisplatin_v7/$sample.feather"

  if [ -f "$outfile" ]; then
    echo "$outfile exists"
  else
    Rscript ../parse_cnv/R/format_segments.R \
      --segment_file $filename \
      --outfname $outfile \
      >& ../../data-SA535/dlp-cisplatin_v7/log/$sample.log
  fi    
done  

## SA1035
for filename in ../../data-SA1035/dlp/csv-v6/*.csv; do
  sample=$(basename "$filename" .csv)
  echo "Sample $sample"
  outfile="../../data-SA1035/dlp/$sample.feather"

  if [ -f "$outfile" ]; then
    echo "$outfile exists"
  else
    Rscript ../parse_cnv/R/format_segments.R \
      --segment_file $filename \
      --outfname $outfile \
      >& ../../data-SA1035/dlp/$sample.log
  fi    
done    



## MA 13 Apr 2020
# Also make the gene_clone_cn files with purity threshold 0 (don't throw our any genes)
# These will not be used by clonealign, but will be used in the differential expression analysis
for filename in ../../data/dlp/gene_cell_cn/*.feather; do
  sample=$(basename "$filename" .feather)
  echo "Sample $sample"
  outfig="../../data/dlp/gene_clone_cn_ignore_filters/$sample.png"
  outfile="../../data/dlp/gene_clone_cn_ignore_filters/$sample.tsv"

  if [ -f "$outfile" ]; then
    echo "$outfile exists"
  else
    Rscript ../parse_cnv/R/select_genes-same-timepoint.R \
      --gene_cn $filename \
      --pct_pure 0.0 \
      --ignore_filters Yes \
      --outfig $outfig \
      --outfname $outfile \
      >& ../../data/dlp/gene_clone_cn_ignore_filters/$sample.log
  fi    
done    

