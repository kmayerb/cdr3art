# April 11, 2023
# Author: kmayerb
# Commandline program for concise diagrams of T cell receptor junctions

library("argparser")
source('R/sankey_and_logo.R')
parser <- arg_parser(description='cdr3art - concise visualization of T cell receptor clusters')
parser <- add_argument(parser, arg="--input", type="character", help = "path to input .csv", default = 'data/motif_graphic_instructions.csv')
parser <- add_argument(parser, arg="--outdir",  type="character", help = "output folder", default =  "./output")
parser <- add_argument(parser, arg="--combine",  type="boolean", help = "combine graphics into a single file", default = TRUE)
parser <- add_argument(parser, arg="--k_per_page",  type="double", help = "number of graphics per page", default = 10)
parser <- add_argument(parser, arg="--individual",  type="boolean", help = "make individual graphics one per cluster", default = FALSE)
parser <- add_argument(parser, arg="--width",  type="double", help = "graphic width (inches)", default = 4)
parser <- add_argument(parser, arg="--height",  type="double", help = "graphic height (inches)", default = 2)
parser <- add_argument(parser, arg="--font_size",  type="double", help = "font size of gene names", default = 4)
parser <- add_argument(parser, arg="--axis_font_size",  type="double", help = "font size of logo axis text)", default = 6)
args = parse_args(parser)
print(args)
df = readr::read_csv(args$input)
chains = validate_input(df)

# SELECT APPROPRIATE FUNCTION, IF WE HAVE BOTH ALPHA AND BETA, OR ONLY SINGLE CHAIN
if (chains$a == T & chains$a == T){
  f = paired_chain_diagram
}else{
  f = single_chain_diagram
}

# GENE COLORS ARE SELECTED BASED ON IMGT V GENE
gene_colors3 = get_gene_colors()

# THE INPUT MUST HAVE NUMERIC CLUSTER IDs
cluster_ids = sort(unique(df$cluster_id))

# IF INDIVIDUAL IS TRUE INDIVIDUAL GRAPHICS WILL BE GENERATED FOR EACH CLUSTER
if (args$individual == TRUE){
  for (i in cluster_ids){
    cluster_data = df[df$cluster_id == i,]
    pdf_name = file.path(args$outdir, paste0("cluster_",i,'.pdf'))
    # WRITE EACH INDIVIDUAL GRAPHIC TO ITS OWN FILE
    print(paste0('generating ', pdf_name))
    pdf(pdf_name, width = args$width, height = args$height)
    f(external_data=cluster_data, 
      gene_colors = gene_colors3, 
      my_title = as.character(i),
      font_size = args$font_size, 
      axis_font_size = args$axis_font_size)
    dev.off()
  }
}


# IF COMBINE IS TRUE, 10 LOGOS PER PAGE WILL BE GENERATED, 
if (args$combine == TRUE){
  # HOW MANY CLUSTERS ARE THEER  
  n <- length(cluster_ids)
  # K is SET TO 10, BECAUSE THERE IS CURRENTLY NOT A WAY TO CHANGE THIS PROGRAMATICALLY
  if (args$k_per_page > 10){
    print("k must be 10 or less")
    k=10
  }else{
    k = args$k_per_page
  }
  
  
  cluster_per_page = split(cluster_ids, rep(1:ceiling(n/k), each=k)[1:n])
  
  counter = 0
  for (x in cluster_per_page){
    print(x)
    cluster_start = min(x)
    cluster_end   = max(x)
    
    mm = df %>% 
      filter(cluster_id %in% x)%>%
      group_by(cluster_id) %>%
      group_split()  
    
    ns = paste0("cluster ", x)
    multiple_motifs = purrr::map2(mm,ns, ~f(external_data = .x, 
                                            gene_colors = gene_colors3, 
                                            my_title = .y),
                                            font_size = args$font_size, 
                                            axis_font_size = args$axis_font_size)
    
    num_motifs = length(multiple_motifs)
    
    width = args$width
    height = args$height*k
    
    fig_name = file.path(args$outdir, paste0("cluster_", cluster_start, "_to_", cluster_end, "_w" , width, "_x_h", height, "___" , n ,"___", "motifs.pdf"))
    print(fig_name)
    
    # WE NEED TO FILL EMPTY ELEMENTS WITH THE LIST WITH SOMTHING BLANK
    if (num_motifs < k){
      for (j in  (num_motifs+1):k){
        multiple_motifs[[j]] = ggplot(data.frame())
      }
    }
    
    

    if (k == 10){
      pdf(fig_name, width = width, height = height )
      gridExtra::grid.arrange(multiple_motifs[[1]],
                              multiple_motifs[[2]],
                              multiple_motifs[[3]],
                              multiple_motifs[[4]],
                              multiple_motifs[[5]],
                              multiple_motifs[[6]],
                              multiple_motifs[[7]],
                              multiple_motifs[[8]],
                              multiple_motifs[[9]],
                              multiple_motifs[[10]],
                              ncol = 1, top = paste0("cluster_", cluster_start, "_to_", cluster_end))
      dev.off()
    }
    if (k == 9){
      pdf(fig_name, width = width, height = height )
      gridExtra::grid.arrange(multiple_motifs[[1]],
                              multiple_motifs[[2]],
                              multiple_motifs[[3]],
                              multiple_motifs[[4]],
                              multiple_motifs[[5]],
                              multiple_motifs[[6]],
                              multiple_motifs[[7]],
                              multiple_motifs[[8]],
                              multiple_motifs[[9]],
                              ncol = 1, top = paste0("cluster_", cluster_start, "_to_", cluster_end))
      dev.off()
    }
    if (k == 8){
      pdf(fig_name, width = width, height = height )
      gridExtra::grid.arrange(multiple_motifs[[1]],
                              multiple_motifs[[2]],
                              multiple_motifs[[3]],
                              multiple_motifs[[4]],
                              multiple_motifs[[5]],
                              multiple_motifs[[6]],
                              multiple_motifs[[7]],
                              multiple_motifs[[8]],
                              ncol = 1, top = paste0("cluster_", cluster_start, "_to_", cluster_end))
      dev.off()
    }
    if (k == 7){
      pdf(fig_name, width = width, height = height )
      gridExtra::grid.arrange(multiple_motifs[[1]],
                              multiple_motifs[[2]],
                              multiple_motifs[[3]],
                              multiple_motifs[[4]],
                              multiple_motifs[[5]],
                              multiple_motifs[[6]],
                              multiple_motifs[[7]],
                              ncol = 1, top = paste0("cluster_", cluster_start, "_to_", cluster_end))
      dev.off()
    }
    if (k == 6){
      pdf(fig_name, width = width, height = height )
      gridExtra::grid.arrange(multiple_motifs[[1]],
                              multiple_motifs[[2]],
                              multiple_motifs[[3]],
                              multiple_motifs[[4]],
                              multiple_motifs[[5]],
                              multiple_motifs[[6]],
                              ncol = 1, top = paste0("cluster_", cluster_start, "_to_", cluster_end))
      dev.off()
    }
    if (k == 5){
      pdf(fig_name, width = width, height = height )
      gridExtra::grid.arrange(multiple_motifs[[1]],
                              multiple_motifs[[2]],
                              multiple_motifs[[3]],
                              multiple_motifs[[4]],
                              multiple_motifs[[5]],
                              ncol = 1, top = paste0("cluster_", cluster_start, "_to_", cluster_end))
      dev.off()
    }
    if (k == 4){
      pdf(fig_name, width = width, height = height )
      gridExtra::grid.arrange(multiple_motifs[[1]],
                              multiple_motifs[[2]],
                              multiple_motifs[[3]],
                              multiple_motifs[[4]],
                              ncol = 1, top = paste0("cluster_", cluster_start, "_to_", cluster_end))
      dev.off()
    }
    if (k == 3){
      pdf(fig_name, width = width, height = height )
      gridExtra::grid.arrange(multiple_motifs[[1]],
                              multiple_motifs[[2]],
                              multiple_motifs[[3]],
                              ncol = 1, top = paste0("cluster_", cluster_start, "_to_", cluster_end))
      dev.off()
    }
    if (k == 2){
      pdf(fig_name, width = width, height = height )
      gridExtra::grid.arrange(multiple_motifs[[1]],
                              multiple_motifs[[2]],
                              ncol = 1, top = paste0("cluster_", cluster_start, "_to_", cluster_end))
      dev.off()
    }
    if (k == 1){
      pdf(fig_name, width = width, height = height )
      gridExtra::grid.arrange(multiple_motifs[[1]],
                              ncol = 1, top = paste0("cluster_", cluster_start, "_to_", cluster_end))
      dev.off()
    }
    
  }
  
}


#Rscript myprogram.R arg1 arg2 arg3
#args <- commandArgs(trailingOnly = TRUE)
#print(args)

