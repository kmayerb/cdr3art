#install.packages('ggseqlogo')
require(ggseqlogo)
require(ggsankey)
require(ggplot2)
require(dplyr)
require(purrr)
require(gridExtra)
require(readr)

#' cdr3art -- uses a data.frame input with the following columns 
#' index - c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "X", "-")
#' x1...xn are the frequency of adjusted frequencies
#' "cluster_id"  integer allows one to pile many motifs into a single file
#' "type" -- "raw", "subtracted", "gene_usage"
#' "chain" -- "a" or "b"
#' "v_a_gene"  (optional if chain is "b") - used in geneusage
#' "j_a_gene"  (optional if chain is "b") - used in geneusage
#' "v_b_gene"  (optional if chain is "a") - used in geneusage
#' "j_b_gene"  (optional if chain is "a") - used in geneusage
#' "x0"
#' "x1"
#' "x2"
#' "x3"
#' "x4"
#' "x5"
#' "x6"
#' "x7"
#' "x8"
#' "x9"
#' "x10"
#' "x11"
#' "x12"
#' "x13"
#' "x14"

validate_input <- function(df){
  # Make sure input has index, type, and chain columns
  # return whether this should be alpha, beta, or paired motif
  stopifnot('index' %in% names(df))
  stopifnot('type' %in% names(df))
  stopifnot('cluster_id' %in% names(df))
  stopifnot(all(unique(df$type) %in% c('raw','subtracted','gene_usage')))
  stopifnot('index' %in% names(df))
  
  stopifnot('chain' %in% names(df))
  stopifnot(all(unique(df$chain) %in% c('a','b', NA)))
  
  if ('a' %in% df$chain){a = TRUE}else(a = FALSE)
  if ('b' %in% df$chain){b = TRUE}else{b = FALSE}
  
  return(list('a' = a, 'b' = b))
}



# Using filename <f> load a gene table that can be used for the color map
get_gene_colors <- function(f= 'data/gene_table.csv'){
  genes = readr::read_csv(f)
  # select only *01 version
  genes = genes %>% 
    dplyr::mutate(primary = stringr::str_detect(string = id, pattern ="\\*01")) %>% 
    dplyr::filter(primary)
  # define a function that will add rainbow colors, which will be applied to each 
  # grouped subset of a dataframe, allowing the same pallete to be applied to TRAV, TRAJ, TRBV, TRBJ
  add_colors <- function(df){
    df$color = rainbow(dim(df)[1])
    return(df)
  }
  # <genes_2> data.frame with color information
  genes_2 = genes %>% 
    dplyr::group_by(chain, region) %>% 
    dplyr::group_split() %>% 
    purrr::map(., ~add_colors(.x) ) %>% 
    do.call(rbind, .)
  
  # <gene_colors> named vector will be passed to ggplot later one
  gene_colors2 = genes_2$color
  names(gene_colors2) = genes_2$id
  # Generate a legend
  legend_plot = ggplot(genes_2, aes(y = forcats::fct_rev(forcats::fct_reorder(id, family_int )), x = 1, col = id)) + geom_point() + 
    scale_color_manual(values = gene_colors2) + 
    facet_wrap(region~chain, scale = "free", ncol = 4)+ 
    theme_classic()  +
    ylab("")+ 
    xlab("")+
    theme(axis.text.x = element_blank())+
    theme(axis.text.y = element_text(size = 6, face= "bold"))+
    theme(legend.position = "None")
  
  # Assign color to abreviation names
  abv_names = gsub(names(gene_colors2), pattern = "TR", replacement = "") %>% 
    gsub(., pattern = "\\*01",replacement = "") 
  
  gene_colors3 = gene_colors2
  names(gene_colors3) = abv_names
  
  genes_2$abv = abv_names
  legend_plot3 = ggplot(genes_2, aes(y = forcats::fct_rev(forcats::fct_reorder(abv, family_int)), 
                      x = 1, 
                      col = abv )) + 
    geom_point() + 
    scale_color_manual(values = gene_colors3) + 
    facet_wrap(region~chain, scale = "free", ncol = 4)+ 
    theme_classic()  +
    ylab("")+ 
    xlab("")+
    theme(axis.text.x = element_blank())+
    theme(axis.text.y = element_text(size = 6, face= "bold"))+
    theme(legend.position = "None")
  
  return(gene_colors3)
}


#' single_chain_diagram
#'
#' @param external_data 
#' @param gene_colors 
#' @param my_title 
#' @param chain 
#' @param font_size 
#' @param axis_font_size 
#'
#' @return ggplot 
#' @export
#'
#' @examples
#' gene_colors3 = get_gene_colors()
#' data = readr::read_csv("data/motif_graphic_instructions.csv")
#' data %>% filter(cluster_id == 0) %>% 
#'  single_chain_diagram(., 
#'                 gene_colors = gene_colors3,
#'                 chain = "b",
#'                 font_size = 4,
#'                 my_title = "CLUSTER" )
single_chain_diagram <- function(external_data, 
                            gene_colors,
                            my_title = "", 
                            font_size = 4, 
                            axis_font_size = 6){
  
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  theme_extras = theme(legend.position = "none") + 
    theme(axis.line = element_line()) + 
    theme(axis.text = element_text(size= axis_font_size , face = "bold")) + 
    theme(axis.title.y = element_text(size=axis_font_size , face = "bold"))
  
  mr = external_data %>% filter(type == "raw", chain == chain)
  m  = external_data %>% filter(type == "subtracted", chain == chain) %>%
    mutate()
  
  gu = external_data %>% 
    filter(type == "gene_usage", chain == chain) %>%
    select(V = v_b_gene, 
           J=j_b_gene) %>%
    mutate(V_abr = gsub(V, pattern = "TR", replacement = "") %>% gsub(., pattern = "\\*01", replacement = ""),
             J_abr = gsub(J, pattern = "TR", replacement = "") %>% gsub(., pattern = "\\*01", replacement = ""))

  df <- gu %>%
    make_long(V_abr,J_abr)
  
  df$node_f = forcats::fct_rev(forcats::fct_infreq(df$node))
  
  index_b = external_data %>% filter(type == "raw", chain == chain) %>% pull(index)
  
  gg_tr = ggplot(df, aes(x = x, 
                           next_x = next_x, 
                           node = node_f, 
                           next_node = next_node,
                           fill = node_f)) +
    geom_sankey(flow.fill = "gray") +
    scale_fill_manual(values=gene_colors) +
    geom_sankey_text(aes(label= node), size = font_size)+
    theme_sankey(base_size = 16) + 
    theme(axis.text.x = element_blank())+ # element_text(angle = 90, size= 6)) + 
    theme(legend.position = "none")+
    xlab("")
  
  mb = external_data %>% 
    filter(type == "subtracted", chain == chain) %>%
    select(contains("x")) %>%
    select(-index) %>%
    select_if(~sum(!is.na(.)) > 0) %>%
    as.matrix()
  
  mbr = external_data %>% 
    filter(type == "raw", chain == chain) %>%
    select(contains("x")) %>%
    select(-index) %>%
    select_if(~sum(!is.na(.)) > 0) %>%
    as.matrix()
  mbr[mbr<0] = 0

  dimnames(mb) = list(index_b)
  dimnames(mbr) = list(index_b)
  
  ggb = ggseqlogo(mb, method='custom', seq_type='aa')   + ylab('Bits') + scale_y_continuous(labels=scaleFUN) + theme_extras 
  ggbr = ggseqlogo(mbr, method='probability', seq_type='aa') + ylab('Probability') + scale_y_continuous(labels=scaleFUN) + theme_extras 
  
    
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  layout_matrix = matrix(c(1,1, 2,2,2,2, 
                           1,1, 3,3,3,3), ncol = 6, byrow = T)
  chain_fig = gridExtra::grid.arrange(
    gg_tr + theme(plot.margin=unit(c(0,0,0,0), "cm")),
    ggb    + theme(plot.margin=unit(c(0,0,0,0), "cm")),
    ggbr   + theme(plot.margin=unit(c(0,0,0,0), "cm")),
    ncol=2, 
    layout_matrix = layout_matrix,
    top = my_title)
  
  return(chain_fig)
}



#' paired_chain_diagram
#'
#' The six part plot : 
#' (1) TRAVJ Sankey, 
#' (2) CDR3A (background subtracted motif)
#' (3) TABVJ Sankey 
#' (4) CDR3B (background_subtracted motif)
#' (5) CDRA (raw motif)
#' (6) CDRA ()
#' @param external_data - must contain columns x1..xN for allignment positions 
#' and <type> column (raw,subtracted,gene_usage)
#' and <chain> column (a,b)
#' and <cluster_id>
#' this single dataframe allows data to be simply ported out of Python tcrdist3 routine
#'
#' @return six_part_plot_fig a collection of grobs from gridExtra::grid.arrange
#' @export 
#'
#' @examples
paired_chain_diagram <- function(external_data, 
                          gene_colors = gene_colors3, 
                          abbr = T,
                          my_title = "",
                          font_size = 3, 
                          axis_font_size = 6){
  mr_a  = external_data %>% filter(type == "raw", chain == "a") %>% select(contains('x'))
  m_a   = external_data %>% filter(type =="subtracted", chain == "a") %>% select(contains('x'))
  mr_b  = external_data %>% filter(type == "raw", chain == "b") %>% select(contains('x'))
  m_b   = external_data %>% filter(type =="subtracted", chain == "b") %>% select(contains('x'))
  gu_a  = external_data %>% filter(type =="gene_usage") %>% select(TRAV = v_a_gene, TRAJ =j_a_gene)
  gu_b  = external_data %>% filter(type =="gene_usage") %>% select(TRBV = v_b_gene, TRBJ=j_b_gene)
  
  
  gu_a = gu_a %>% mutate(TRAV_abr = gsub(TRAV, pattern = "TR", replacement = "") %>% gsub(., pattern = "\\*01", replacement = ""),
                         TRAJ_abr = gsub(TRAJ, pattern = "TR", replacement = "") %>% gsub(., pattern = "\\*01", replacement = "")) 
  gu_b = gu_b %>% mutate(TRBV_abr = gsub(TRBV, pattern = "TR", replacement = "") %>% gsub(., pattern = "\\*01", replacement = ""),
                         TRBJ_abr = gsub(TRBJ, pattern = "TR", replacement = "") %>% gsub(., pattern = "\\*01", replacement = "")) 
  
  
  
  
  if (abbr){
    dfa <- gu_a%>%
      make_long(TRAV_abr,TRAJ_abr)
    dfb <- gu_b%>%
      make_long(TRBV_abr,TRBJ_abr)
  }else{
    dfa <- gu_a%>%
      make_long(TRAV,TRAJ)
    dfb <- gu_b%>%
      make_long(TRBV,TRBJ)
  }
  
  dfa$node_f = forcats::fct_rev(forcats::fct_infreq(dfa$node))
  
  # Make Sankey For TRAV-TRAJ
  gg_tra = ggplot(dfa, aes(x = x, 
                           next_x = next_x, 
                           node = node_f, #forcats::fct_infreq(node), 
                           next_node = next_node,
                           fill = node_f)) +
    geom_sankey(flow.fill = "gray") +
    scale_fill_manual(values=gene_colors) +
    geom_sankey_text(aes(label= node), size = font_size)+
    theme_sankey(base_size = 16) + 
    theme(axis.text.x = element_blank()) +   
    theme(legend.position = "none")+
    xlab("")
  
  
  dfb$node_f = forcats::fct_rev(forcats::fct_infreq(dfb$node))
  
  gg_trb = ggplot(dfb, aes(x = x, 
                           next_x = next_x, 
                           node = node_f, 
                           next_node = next_node,
                           fill = node_f)) +
    geom_sankey(flow.fill = "gray") +
    scale_fill_manual(values=gene_colors) +
    geom_sankey_text(aes(label= node), size = font_size)+
    theme_sankey(base_size = 16) + 
    theme(axis.text.x = element_blank())+ # element_text(angle = 90, size= 6)) + 
    theme(legend.position = "none")+
    xlab("") #+
  
  
  # Create a custom matrix 
  set.seed(123)
  
  index_a = external_data %>% filter(type == "raw", chain == "a") %>% pull(index)
  index_b = external_data %>% filter(type == "raw", chain == "a") %>% pull(index)
  all(index_a == index_b)
  
  mar = external_data %>% filter(type == "raw", chain == "a") %>%
    select(contains("x")) %>%
    select(-index) %>%
    select_if(~sum(!is.na(.)) > 0) %>%
    as.matrix()
  dimnames(mar) = list(index_a)

  
  mbr = external_data %>% filter(type == "raw", chain == "b") %>%
    select(contains("x")) %>%
    select(-index) %>%
    select_if(~sum(!is.na(.)) > 0) %>%
    as.matrix()
  dimnames(mbr) = list(index_b)
  
  ma = external_data %>% filter(type == "subtracted", chain == "a") %>%
    select(contains("x")) %>%
    select(-index) %>%
    select_if(~sum(!is.na(.)) > 0) %>%
    as.matrix()
  dimnames(ma) = list(index_a)
  
  mb = external_data %>% filter(type == "subtracted", chain == "b") %>%
    select(contains("x")) %>%
    select(-index) %>%
    select_if(~sum(!is.na(.)) > 0) %>%
    as.matrix()
  dimnames(mb) = list(index_b)
  
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  
  theme_extras = theme(legend.position = "none") + 
    theme(axis.line = element_line()) + 
    theme(axis.text = element_text(size= axis_font_size, face = "bold")) + 
    theme(axis.title.y = element_text(size= axis_font_size, face = "bold"))
  

  gga = ggseqlogo(ma, method='custom', seq_type='aa')   + ylab('Bits') + scale_y_continuous(labels=scaleFUN) + theme_extras 
  mar[mar<0] = 0
  
  ggar = ggseqlogo(mar, method='probability', seq_type='aa') + ylab('Probability') + scale_y_continuous(labels=scaleFUN) + theme_extras 
  
  ggb = ggseqlogo(mb, method='custom', seq_type='aa')   + ylab('Bits') + scale_y_continuous(labels=scaleFUN) + theme_extras 
  mbr[mbr<0] = 0
  
  ggbr = ggseqlogo(mbr, method='probability', seq_type='aa') + ylab('Probability') + scale_y_continuous(labels=scaleFUN) + theme_extras 
  
  
  layout_matrix = matrix(c(1,1, 2,2,2,2, 3,3, 4,4,4,4,
                           1,1, 5,5,5,5, 3,3, 6,6,6,6), 
                         ncol = 12, byrow = T)
  six_part_plot_fig = gridExtra::grid.arrange(
    gg_tra + theme(plot.margin=unit(c(0,0,0,0), "cm")),
    gga    + theme(plot.margin=unit(c(0,0,0,0), "cm")),
    gg_trb + theme(plot.margin=unit(c(0,0,0,0), "cm")),
    ggb    + theme(plot.margin=unit(c(0,0,0,0), "cm")),
    ggar   + theme(plot.margin=unit(c(0,0,0,0), "cm")), 
    ggbr   + theme(plot.margin=unit(c(0,0,0,0), "cm")),
    ncol= 4, layout_matrix =layout_matrix, top = my_title)
  return(six_part_plot_fig )
}


