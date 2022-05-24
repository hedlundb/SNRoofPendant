## plot_cladogram_modified.R
#  Cale Seymour
#  University of Nevada Las Vegas
#  Updated August 2020
#  Description: Modification to the plot_cladogram function from the microbiomeMarker package.
#  Allows the user to:
#  1) specify to ignore ranks (i.e., don't plot LEFSE results at the species and genus level)
#  2) Exclude "unclassified" taxa when plotting.
#  The purpose of both of these is to reduce the noisiness of the LEFSE trees, especially in larger datasets.
#  Usage: (in R) source('plot_cladogram_modified.R').
#  New params: <character length 1:[number of ranks], within c('p','c','o','f','g','s')> "exclude_rank"
#              <logical length 1> "drop_unclassifieds"
#  Usage example:
#   # clado = plot_cladogram_mod(lef, color = col, exclude_rank = 's', drop_unclassifieds = TRUE, clade_label_level = 5)
## Please cite Friel et al. (2022) and the R package microbiomeMarker if you use this script verbatim.

library('tidytree')
library('ape')
library('ggplot2')
library('ggtree')

## Source this file to get access to the modified version of the function.

plot_cladogram_mod = function (mm, color, exclude_rank = NULL, drop_unclassifieds = FALSE, branch_size = 0.2, alpha = 0.2, node_size_scale = 1, 
                               node_size_offset = 1, clade_label_level = 1, annotation_shape = 22, 
                               annotation_shape_size = 5, group_legend_param = list(), marker_legend_param = list()) 
{
    ps <- phyloseq(mm@otu_table, mm@tax_table)
    
    ## Modifications start here. Drop tips associated with the Species column of
    ## the phyloseq object.
    
    treedata = get_treedata_phyloseq(ps)
    
    if (!is.null(exclude_rank))
    {
        ## Identify offending nodes of the node class which we want to drop.
        offending = treedata@data$node[(treedata@data$node_class %in% exclude_rank)]
        
        newtree = ape::drop.tip(treedata@phylo, offending, trim.internal=FALSE, collapse.singles=FALSE)
        newdata = treedata@data[!treedata@data$node %in% offending,]
        newdata$node = match(newdata$node_label,c(newtree$tip.label, newtree$node.label))
        newdata$node_class = factor(newdata$node_class, levels = levels(newdata$node_class)[!levels(newdata$node_class) == 's'])
        
    } else {
        newtree = treedata@phylo
        newdata = treedata@data
    }
    
    if (drop_unclassifieds == TRUE)
    {
        offending = grep('[Uu]nclassified', newtree$tip.label)
        while(length(offending) > 0)
        {
            newtree = ape::drop.tip(newtree, offending, trim.internal=FALSE, collapse.singles=FALSE)
            newdata = newdata[!newdata$node %in% offending,]
            newdata$node = match(newdata$node_label,c(newtree$tip.label, newtree$node.label))
            offending = grep('[Uu]nclassified', newtree$tip.label)
        }
    }
    
    newtreedata = tidytree::treedata(phylo = newtree, data = newdata)
    tree = generate_taxa_tree(newtreedata, size = branch_size)
    
    
    ## Screw with the annotations a little bit.
    excludedNodes = setdiff(treedata@data$node_label, newtreedata@data$node_label)
    
    annotation = generate_cladogram_annotation(mm@marker_table, 
                                               color = color)
    annotation = annotation[!annotation$node %in% excludedNodes,]
    ## End modifications.
    
    annotation_info <- dplyr::left_join(annotation, tree$data, 
                                        by = c(node = "label")) %>% mutate(label = .data$node, 
                                                                           id = .data$node.y, level = as.numeric(.data$node_class))
    hilight_para <- dplyr::transmute(annotation_info, node = .data$id, 
                                     fill = .data$color, alpha = alpha, extend = get_offset(.data$level))
    hilights_g <- purrr::pmap(hilight_para, geom_hilight)
    tree <- purrr::reduce(hilights_g, `+`, .init = tree)
    hilights_df <- dplyr::distinct(annotation_info, .data$enrich_group, 
                                   .data$color) %>% arrange(.data$enrich_group)
    hilights_df$x <- 0
    hilights_df$y <- 1
    group_legend_param <- c(group_legend_param, list(title = NULL, 
                                                     order = 1, override.aes = list(fill = hilights_df$color)))
    group_lgd <- do.call(guide_legend, group_legend_param)
    tree <- tree + geom_rect(aes_(xmin = ~x, xmax = ~x, ymax = ~y, 
                                  ymin = ~y, fill = ~enrich_group), data = hilights_df, 
                             inherit.aes = FALSE) + guides(fill = group_lgd)
    nodes_colors <- rep("white", nrow(tree$data))
    nodes_colors[annotation_info$id] <- annotation_info$color
    node_size <- node_size_scale * log(tree$data$abd) + node_size_offset
    tree$data$node_size <- node_size
    tree <- tree + geom_point2(aes(size = I(node_size)), fill = nodes_colors, 
                               shape = 21)
    clade_label <- dplyr::transmute(annotation_info, node = .data$id, 
                                    offset = get_offset(.data$level) - 0.4, angle = purrr::map_dbl(.data$id, 
                                                                                                   get_angle, tree = tree) + 90, label = .data$label, 
                                    fontsize = 1.5 + sqrt(.data$level), barsize = 0, hjust = 0.5, 
                                    level = .data$level) %>% dplyr::arrange(desc(.data$level))
    ind <- clade_label$level < clade_label_level
    short_label <- get_short_label_id(clade_label, clade_label_level)
    clade_label_para <- mutate(clade_label, label = c(.data$label[!ind], 
                                                      short_label), level = NULL)
    clade_label_g <- purrr::pmap(clade_label_para, geom_cladelabel)
    tree <- purrr::reduce(clade_label_g, `+`, .init = tree)
    guide_label <- clade_label[ind, ] %>% mutate(label2 = paste0(short_label, 
                                                                 ": ", .data$label), color = annotation_info$color[match(.data$label, 
                                                                                                                         annotation_info$label)])
    marker_legend_param <- c(marker_legend_param, list(p = tree, 
                                                       color = guide_label$color, label = guide_label$label2, 
                                                       shape = annotation_shape, size = annotation_shape_size))
    p <- do.call(set_marker_annotation, marker_legend_param) + 
        theme(legend.position = "right", legend.title = element_blank())
    p
}

## Set the function within the microbiomeMarker namespace.
environment(plot_cladogram_mod) = asNamespace('microbiomeMarker')
