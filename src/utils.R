# Note, many of the functions below are implemented in the microbiomer-package.
# See here: https://github.com/wsteenhu/microbiomer

# Save/load function for 'heavy' objects
sl <- function(name, ..., overwrite = FALSE, dir_path = here::here("results", "RData", subdir_name)) {
  # Possibility to add name as name or literal character string
  name <- as.character(substitute(name))
  assign(name, 
    if(file.exists(glue::glue("{dir_path}/{name}.Rds")) && !overwrite) {
     readRDS(glue::glue("{dir_path}/{name}.Rds"))
    }
    else { 
     dir.create(dir_path, showWarnings = F, recursive = T)
     saveRDS(..., file=glue::glue("{dir_path}/{name}.Rds"))
     readRDS(glue::glue("{dir_path}/{name}.Rds"))
    }, envir=.GlobalEnv)
}

# Format printing numbers with comma separating thousands
f <- function(x, ...) {
  format(x = x, big.mark=",", ...)
}

# Round
r <- function(x, ...) {
  round(x = x, digits = 3, ...)
}

# Keep trailing zero's
# Significant figures: # http://stackoverflow.com/questions/3245862/
s <- function(vec, digits = 3, format = "f", flag = ""){
  return(gsub("\\.$", "", 
              formatC(vec, 
                      digits = digits, 
                      # use "fg" for significant digits
                      format = format, 
                      flag = flag)))
}


# Conversion functions

to_RA <- function(ps) {
  phyloseq::transform_sample_counts(ps, function(OTU) OTU / sum(OTU))
}

pres_abund_filter <- function(ps, pres = 2, abund = 0.001, verbose = TRUE) { # Subramanian filter
  is_raw <- max(phyloseq::otu_table(ps)) > 1 # detect is ps has raw reads
  
  if(is_raw) { # if raw reads; convert to RA for filtering
    ps_raw <- ps
    ps <- ps_raw %>% to_RA()
  }
  
  ps_filt <- phyloseq::filter_taxa(ps, function(x) sum(x > abund) >= pres, TRUE)
  
  if(verbose) {
    message(glue::glue("A total of {phyloseq::ntaxa(ps_filt)} ASVs were found to be present at or above a level of confident detection ({(abund * 100)}% relative abundance) in at least {pres} samples (n = {phyloseq::ntaxa(ps) - phyloseq::ntaxa(ps_filt)} ASVs excluded)."))
  }
  
  if(!is_raw) {
    return(ps_filt)
  } else {
    return(phyloseq::prune_taxa(phyloseq::taxa_names(ps_filt), ps_raw))
  }
}

mean_abund_filter <- function(ps, mean_abund = 0.001) {
  filter_taxa(ps, function(x) mean(x) > mean_abund, TRUE) } 

asinsqrt <- function(otu_table) {
  asin(sqrt(otu_table))
}

get_topn <- function(ps, n = 15, residuals = TRUE) {
  otu_tab <- phyloseq::otu_table(ps)
  otu_tab_n <- otu_tab[order(rowMeans(otu_tab), decreasing = TRUE)[1:n], ]
  otu_tab_res <- otu_tab[-order(rowMeans(otu_tab), decreasing = TRUE)[1:n], ]
  
  if(residuals) {
    otu_tab_n <- otu_tab_n %>%
      t %>%
      data.frame(residuals = colSums(otu_tab_res), check.names = FALSE) %>%
      t
  }
  return(phyloseq::phyloseq(phyloseq::otu_table(otu_tab_n, taxa_are_rows = TRUE),
                            phyloseq::sample_data(ps)))
}

prep_bar <- function(ps, n, residuals = TRUE) {
  excl_cols <- c("sample_id2", colnames(phyloseq::sample_data(ps)))
  
  df_topn <- ps %>%
    get_topn(n = n, residuals = residuals) %>%
    ps_to_df(sample_name = "sample_id2") %>%
    tidyr::pivot_longer(-dplyr::all_of(excl_cols), names_to = "OTU", values_to = "value") %>%
    dplyr::mutate(OTU = format_OTU(.data$OTU) %>% forcats::fct_inorder() %>% forcats::fct_rev()) %>%
    dplyr::arrange(.data$sample_id2) %>%
    dplyr::mutate(sample_id2 = forcats::fct_inorder(.data$sample_id2))
  
  return(df_topn)
}

otu_tab_to_df <- function(ps, sample_name = "sample_id2") {
  otu_tab <- phyloseq::otu_table(ps)
  df <- otu_tab %>%
    t %>%
    data.frame() %>%
    tibble::rownames_to_column(sample_name)
  return(df)
}

meta_to_df <- function(ps, sample_name = "sample_id2") {
  meta <- phyloseq::sample_data(ps)
  df <- meta %>%
    data.frame() %>%
    tibble::rownames_to_column(sample_name)
  return(df)
}

ps_to_df <- function(ps, sample_name = "sample_id2") {
  df_meta <- meta_to_df(ps, sample_name)
  df_otu <- otu_tab_to_df(ps, sample_name)
  
  dplyr::left_join(df_meta, df_otu, by = sample_name) %>% tibble::as_tibble()
}

log10_px <- function(ps, pseudocount = 1) {
  ps %>% transform_sample_counts(., function(x) log10(x + pseudocount)) }

format_OTU <- function(OTU_names, short = F, parse = F) {
  
  OTU_names_tb <- tibble::tibble(OTU_names)
  
  OTU_names_num <- OTU_names_tb %>%
    dplyr::mutate(num = purrr::map_chr(OTU_names, ~stringr::str_extract(., "(?<=_)[0-9]+$")))
  
  OTU_names_tb_format <- OTU_names_num %>%
    dplyr::filter(!is.na(.data$num) & !duplicated(OTU_names)) %>%
    dplyr::mutate(
      name1 = unlist(purrr::map(stringr::str_split(OTU_names, "__|_"), ~ head(.x, -1) %>% glue::glue_collapse(., " "))), #everything except last
      name2 = unlist(purrr::map(stringr::str_split(OTU_names, "__|_"), ~ head(.x, 1))), #last
      long_name = as.character(glue::glue("*{name1}* ({num})")),
      short_name = as.character(glue::glue("*{name2}* ({num})")),
      parse_name = as.character(glue::glue("italic('{name2}')~({num})")))
  
  OTU_names_nonum <- OTU_names_num %>%
    dplyr::filter(is.na(.data$num) & !duplicated(OTU_names))
  
  OTU_names_final <- dplyr::left_join(
    OTU_names_tb,
    dplyr::bind_rows(OTU_names_tb_format, OTU_names_nonum) %>%
      dplyr::mutate(dplyr::across(dplyr::ends_with("_name"), ~dplyr::if_else(is.na(.data$num), stringr::str_to_title(OTU_names), .x)))
    , by = "OTU_names")
  
  if(!short & !parse) {
    return(OTU_names_final %>% dplyr::pull(.data$long_name))
  } else if(parse) {
    return(OTU_names_final %>% dplyr::pull(.data$parse_name))
  } else {
    return(OTU_names_final %>% dplyr::pull(.data$short_name)) }
}

# Plotting functions

prep_bar <- function(ps, n, residuals = TRUE) {
  excl_cols <- c("sample_id2", colnames(phyloseq::sample_data(ps)))
  
  df_topn <- ps %>%
    get_topn(n = n, residuals = residuals) %>%
    ps_to_df(sample_name = "sample_id2") %>%
    tidyr::pivot_longer(-dplyr::all_of(excl_cols), names_to = "OTU", values_to = "value") %>%
    dplyr::mutate(OTU = format_OTU(OTU) %>% forcats::fct_inorder() %>% forcats::fct_rev()) %>%
    dplyr::arrange(sample_id2) %>%
    dplyr::mutate(sample_id2 = forcats::fct_inorder(sample_id2))
  
  return(df_topn)
}

make_color_scheme <- function(name, n) {
  max_n <-  RColorBrewer::brewer.pal.info %>%
    tibble::rownames_to_column("name_pal") %>%
    dplyr::filter(.data$name_pal == name) %>%
    dplyr::pull(.data$maxcolors)
  
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(
    dplyr::case_when(n < 3 ~ 3, n < max_n ~ n, TRUE ~ as.numeric(max_n)), name))(n)
}

create_bar <- function(ps = NULL, df_topn = NULL, id = "sample_id2", y = "value", n = 15, ncol_legend = 3, name_legend = "OTU", RA = TRUE, colour = "white") {
  
  # accepts either ps (running prep_bar) or a dataframe already prepared with prep_bar
  if(any(class(ps)=="phyloseq")) { df_topn <- prep_bar(ps = ps, n = n) }
  
  bar <- df_topn %>%
    ggplot2::ggplot(ggplot2::aes_string(x = id, y = y, fill = "OTU")) +
    ggplot2::geom_bar(stat="identity", colour=colour) +
    ggplot2::scale_fill_manual(name = name_legend, values = c("grey90", rev(make_color_scheme("Paired", n)))) +
    ggplot2::labs(x = "" , y = "Number of reads") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   legend.position = "bottom", legend.text = ggtext::element_markdown(),
                   panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = ncol_legend))
  
  if(RA) {
    bar <- bar +
      ggplot2::scale_y_continuous(expand = c(0.01, 0.01), labels = scales::percent) +
      ggplot2::ylab("Relative abundance")
  }
  return(bar)
  #https://bookdown.org/rdpeng/RProgDA/non-standard-evaluation.html
  
  #TODO: add y = value to plot mean values 
}

create_dendro_gg <- function(hc, hang_height=0.05) {
  
  library(dendextend)
  library(ggdendro)
  
  dendro_data_bl <- hc %>% as.dendrogram %>% hang.dendrogram(hang_height=hang_height) %>% dendro_data
  dendro_data_bl$segments$yend[dendro_data_bl$segments$yend<0] <- 0
  
  dendro_data_gr <- hc %>% as.dendrogram %>% dendro_data
  hc_order <- dendro_data_gr$labels$label
  
  plot <- ggplot(segment(dendro_data_gr)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend),colour="grey75") + 
    geom_segment(data=segment(dendro_data_bl),aes(x=x, y=y, xend=xend, yend=yend)) +
    scale_x_continuous(expand = rep(1/length(hc_order)/2, 2)) + 
    scale_y_continuous(expand=c(0,0.02)) +
    #theme(plot.margin=unit(c(0,0,0,0),"lines")) +
    theme_void()
  return(list(plot=plot, hc_order=hc_order)) 
}

meta_bar <- function(data, var, name = NULL, color_scale = NULL, ...) {
  data <- data %>% mutate(index=1:nrow(.))
  
  p <- ggplot(data, aes_string(x = "index", y = 1, fill = var)) + 
    geom_tile() + 
    scale_y_continuous(expand = c(0,0), breaks = 1, labels = name) + 
    scale_x_continuous(expand = c(0,0)) + 
    theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), legend.key.size = unit(0.2, "in"), legend.spacing = unit(0.2, "line"), plot.margin=unit(c(0, 0, 0, 0), "lines"))
  
  if(!is.null(color_scale)) { 
    p <- p + color_scale
  } else {
    n_col <- length(levels(meta_hm_order[[var]])) %>% as.numeric()
    p <- p + scale_fill_manual(values = c(make_color_scheme("BrBG", n_col -1), "grey70"))
  }
  return(p)
}

#Hclust functions adapted from Kadi & Justyna

vegan_otu <-  function(physeq){   #source: https://rdrr.io/github/taowenmicro/ggClusterNet/src/R/utls.R
  OTU <-  otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <-  t(OTU)
  }
  return(as(OTU,"matrix"))
}

stacked.bar=function(ps,n_otus,hc_order=T,method="average",cols=NULL,var1=NULL,var2=NULL, plot=TRUE,by_genus=F, ra=F){
  
  if(by_genus==F){
    if(ra==F){
      otu_RA <- ps %>%
        transform_sample_counts(., function(x) x/sum(x)) 
    }else{otu_RA=ps}
    
    otu_RA= otu_RA%>%
      otu_table
  }else{
    otu_RA= t(genus_table(ps=ps))
  }
  
  otu_RA_ord <- otu_RA[order(rowMeans(otu_RA), decreasing=T), ]
  
  if(hc_order) {
    bc <- vegdist(t(otu_RA_ord), method="bray")
    hc <- hclust(bc, method=method) 
    otu_RA_ord <- otu_RA_ord[, hc$order]
  } 
  
  stack.df <- t(otu_RA_ord[1:n_otus, ]) %>% 
    data.frame(., Residuals=1-rowSums(.), check.names = F) %>%
    rownames_to_column("sample_id2") %>%
    gather(key=OTU, value=RA, -sample_id2) %>% 
    mutate(OTU = fct_rev(fct_inorder(OTU)),
           sample_id2=fct_inorder(sample_id2))
  
  stack.df_meta <- left_join(stack.df, sample_data(ps) %>%
                               data.frame(check.names=F) %>%
                               rownames_to_column("sample_id2") %>%
                               #select(sample_id2, {{var1}},{{var2}}) %>%
                               mutate(sample_id2=as.factor(sample_id2)), by="sample_id2") %>%
    mutate(OTU = format_OTU(OTU) %>% forcats::fct_inorder() %>% forcats::fct_rev())
  
  set.seed(19)
  
  if(is.null(cols)){
    cols <- c("grey90", rev(make_color_scheme("Paired", n=n_otus)))} else{cols=cols}
  #
  
  bar1 <- stack.df_meta %>%
    ggplot(aes(x = sample_id2, y = RA, fill = OTU)) + 
    geom_bar(stat = "identity", position = "stack", colour = "white") + 
    scale_fill_manual(name = "ASV", values = cols) +
    theme(legend.key.size =  unit(0.2, "in"), legend.spacing = unit(0.2, "line"), plot.margin=unit(c(0,0,0,0),"lines"), 
          axis.title.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "bottom",
          legend.text = ggtext::element_markdown()) +
    guides(fill = guide_legend(ncol = 4, reverse = T)) + 
    scale_x_discrete(name="Sample") + 
    scale_y_continuous(labels=scales::percent) +
    ylab("Relative abundance") 

  #if(!is.null(var)){
  #  return(bar1+var_grid(~var, scales = "free_x", space="free_x") )
  #}else{ 
  
  if(plot){return(bar1)
  }else{return(stack.df_meta)
  }}


meta_plot <- function(data, var, name = NULL, fillScale = NULL, ...) {
  data <- data %>% mutate(index=1:nrow(.))
  p <- ggplot(data, aes_string(x="index", y=1, fill=var)) + 
    geom_tile() + 
    scale_y_continuous(expand=c(0,0), breaks=1, labels=name) + 
    scale_x_continuous(expand=c(0,0)) + 
    theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), legend.key.size =  unit(0.2, "in"), legend.spacing = unit(0.2, "line"), plot.margin=unit(c(0,0,0,0),"lines"))
  if(!is.null(fillScale)) p <- p + fillScale
  return(p)
}

#Heatmap cluster plot

cluster_plot_hm_sampletype2_RA =  function(ps,k=NULL, n_otus,min_size=3,clust.var=NULL, plot=T,method="average",by_genus=F,clust_colours=F){
  #1.Create distance metric and dendogram
  if(by_genus==F){ 
    bc=vegdist(vegan_otu(ps),method = "bray")
  }else{
    bc=vegdist(genus_table(ps),method = "bray")
  }
  hc= hclust(bc, method={{method}})
  dendro_full <- create_dendro_gg(hc)
  hc_order <- dendro_full$hc_order
  dendro <- dendro_full$plot
  #2.Create heatmap
  
  hm_data <- stacked.bar(ps, n_otus = n_otus, hc_order = T, method = method, plot=F, ra=T) %>%
    mutate(RA=ifelse(RA==0, NA, RA))
  
  hm <- hm_data %>%
    mutate(sample_id2=factor(sample_id2, levels=hc_order)) %>%
    ggplot(aes(x=sample_id2, y=OTU, fill=RA)) +
    geom_tile(colour="white") +
    scale_fill_gradientn(name="Relative abundance", colors=brewer.pal(9, "Reds"), na.value = "grey90") +
    theme_bw() + 
    theme(axis.text.x = element_blank(),
          axis.text.y = ggtext::element_markdown(size=12),
          axis.title =element_text(size=12), legend.position = "bottom") +
    xlab("Sample") + ylab("ASV")
  
  #3.Add metadata
  if(is.null(clust.var)){
    meta_clust <- hm_data %>%
      dplyr:: slice(match(hc_order, sample_id2)) %>%
      mutate(clust=as.character(cutree(hc, k)[hc$order])) #k is the number of clusters
  }else{
    meta_clust <- hm_data %>%
      dplyr:: slice(match(hc_order, sample_id2))  
    meta_clust$clust=meta_clust[,{{clust.var}}]   #rename cluster variable if provided
    
    #meta_clust <-hm_data
    #meta_clust$clust=meta_clust[,{{clust.var}}]
  }
  
  if(is.null(clust.var)){
    min_size = {{min_size}} # recode small clusters to NA
    meta_clust_NA <- meta_clust %>% 
      group_by(clust) %>%
      mutate(n=n()) %>% ungroup %>%
      mutate(clust=ifelse(n<=min_size, NA, clust) %>% as.factor() %>% fct_infreq())
  }else{
    meta_clust_NA <- meta_clust
  }
  
  cols_clust <- brewer.pal(length(levels(meta_clust_NA$clust)),"Paired")
  clust_plot <- meta_clust_NA %>%
    meta_plot("clust", "Cluster", fillScale=scale_fill_manual(name="Cluster", 
                                      values=cols_clust, na.value="grey"))+
    theme(axis.text.y = element_text(size=12), legend.position = "none")
  
  #add the sample type metadata
  sample_types <- c("MO", "MS", "MN", "BM", "IS", "IN", "CS", "blank")
  sample_type_colors <- c("firebrick2", "darkorange2", "gold2", "springgreen3","dodgerblue2", "purple2","violet","gray40")
  names(sample_type_colors) <- sample_types
  sample_type_colors <- sample_type_colors[sample_types]
  
  visits <- c("V1","V3","V4","V5","V6","blank")
  visit_colors <- c("deeppink3","yellow2","palegreen3","skyblue2","mediumpurple","gray40")
  names(visit_colors) <- visits
  visit_colors <- visit_colors[visits]
  early_late_colors <- c("honeydew1","palegreen","springgreen","springgreen3","seagreen4", "gray40")
  names(early_late_colors) <- visits
  early_late_colors <- early_late_colors[visits]
  sampletype_visits <- c("V3.MO","V4.MO","V5.MO","V6.MO", "V3.MS","V4.MS","V5.MS","V6.MS",
                         "V3.MN","V4.MN","V5.MN","V6.MN","V3.BM","V4.BM","V5.BM","V6.BM",
                         "V3.IS","V4.IS","V5.IS","V6.IS","V3.IN","V4.IN","V5.IN","V6.IN", "V6.CS")
  sampletype_visit_colors <- c("firebrick1","firebrick2","firebrick3","firebrick4", 
                               "orange1","orange2","darkorange2","darkorange3","yellow1","yellow2","yellow3","yellow4",
                               "palegreen","green","green3","green4","deepskyblue","deepskyblue3","dodgerblue3","blue3",
                               "orchid1","mediumorchid2","magenta3","purple3","lightpink1")
  names(sampletype_visit_colors) <- sampletype_visits
  sampletype_visit_colors <- sampletype_visit_colors[sampletype_visits]
  
  sampletype <- meta_clust %>%
    meta_plot("sample_type", "Niche", fillScale=scale_fill_manual(name="", values=sample_type_colors)) +
    theme(axis.text.y = element_text(size=12), legend.position = "right", legend.title=element_blank()) + guides(fill = guide_legend(ncol = 6))
  visitnumber <- meta_clust %>%
    meta_plot("visit", "Timepoint", fillScale=scale_fill_manual(name="", values=visit_colors)) +
    theme(axis.text.y = element_text(size=12), legend.position = "right", legend.title=element_blank()) + guides(fill = guide_legend(ncol = 6))
  earlylate <- meta_clust %>%
    meta_plot("visit", "Timepoint", fillScale=scale_fill_manual(name="", values=early_late_colors)) +
    theme(axis.text.y = element_text(size=12), legend.position = "right", legend.title=element_blank()) + guides(fill = guide_legend(ncol = 5))
  niche.timepoint <- meta_clust %>%
    meta_plot("visit.type", "Timepoint.Niche", fillScale=scale_fill_manual(name="", values=sampletype_visit_colors)) +
    theme(axis.text.y = element_text(size=12), legend.position = "right", legend.title=element_blank()) + guides(fill = guide_legend(ncol = 5))
  inoculatedTF <- meta_clust %>%
    meta_plot("inoculated", "Inoculated") +
    theme(axis.text.y = element_text(size=12), legend.position = "right", legend.title=element_blank()) + 
    scale_fill_manual(values = c("TRUE" = "deepskyblue", "FALSE" = "lightpink")) + guides(fill = guide_legend(ncol = 6))
  siblingTF <- meta_clust %>%
    meta_plot("sib_5yr", "Sibling") +
    theme(axis.text.y = element_text(size=12), legend.position = "right", legend.title=element_blank()) + 
    scale_fill_manual(values = c("TRUE" = "goldenrod2", "FALSE" = "dimgrey")) + guides(fill = guide_legend(ncol = 6))
  
  #5. Combine
  comb_plot <- cowplot::plot_grid(dendro, clust_plot, sampletype, earlylate, inoculatedTF, hm , ncol=1, rel_heights = c(1,0.15,0.15, 0.15,0.15,2), align="v", axis="lr")
  
  if(plot){print(comb_plot)}else{return(meta_clust_NA)}
}
