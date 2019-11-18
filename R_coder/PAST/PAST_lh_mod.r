
library(dplyr)
library(ggplot2)
library(qvalue)
library(utils)
library(doSNOW)

#' Load GWAS data
#'
#' @param association_file The association file
#' @param effects_file  The effects file
#' @param association_columns The names of the columns in your association data
#'   for Trait, Marker, Chromosome, Site, F, p, and marker_Rsquared
#' @param effects_columns The names of the columns in your effects data for
#'   Trait, Marker, Chromosome, Site, and effect
#' @return The association data and the effects data merged into a dataframe
#'   with one row for each SNP
#' @export
#' @importFrom rlang .data
#' @import utils
#' @import dplyr
#' @examples
#' demo_association_file = system.file("extdata", "association.txt.xz",
#'   package = "PAST", mustWork = TRUE)
#' demo_effects_file = system.file("extdata", "effects.txt.xz",
#'   package = "PAST", mustWork = TRUE)
#' gwas_data <- load_GWAS_data(demo_association_file, demo_effects_file)
load_GWAS_data <- function(association_file,
                           effects_file,
                           association_columns = c("Trait",
                                                   "Marker",
                                                   "Locus",
                                                   "Site",
                                                   "p",
                                                   "marker_R2"),
                           effects_columns = c("Trait",
                                               "Marker",
                                               "Locus",
                                               "Site",
                                               "Effect")) {

  stats <- read.table(association_file, header = TRUE, sep = "\t") %>%
    dplyr::mutate(.data,
                  Trait = !!as.name(association_columns[1]),
                  Marker_original = !!as.name(association_columns[2]),
                  Chr = !!as.name(association_columns[3]),
                  Pos = !!as.name(association_columns[4]),
                  Marker = paste0(Chr, "_", Pos),
                  p = !!as.name(association_columns[5]),
                  marker_R2 = !!as.name(association_columns[6])) %>%
    dplyr::select(.data$Marker,
                  .data$Marker_original,
                  .data$Trait,
                  .data$Chr,
                  .data$Pos,
                  .data$p,
                  .data$marker_R2)

  effects <- read.table(effects_file, header = TRUE, sep = "\t") %>%
    dplyr::mutate(.data,
                  Trait = !!as.name(effects_columns[1]),
                  Marker_original = !!as.name(effects_columns[2]),
                  Chr = !!as.name(effects_columns[3]),
                  Pos = !!as.name(effects_columns[4]),
                  Effect = !!as.name(effects_columns[5])) %>%
    dplyr::select(.data$Marker_original,
                  .data$Trait,
                  .data$Chr,
                  .data$Pos,
                  .data$Effect)

  # Delete all markers in effects and stats with more or less alleles than 2
  non_biallelic <- effects %>%
    dplyr::group_by(.data$Marker_original) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::filter(count != 2)
  effects <-
    effects %>% dplyr::filter(!(.data$Marker_original %in% non_biallelic$Marker_original))
  stats <-
    stats %>% dplyr::filter(!(.data$Marker_original %in% non_biallelic$Marker_original))

  # Remove all NaN data to prevent math with NaN
  stats <- stats %>% dplyr::filter(.data$marker_R2 != "NaN")

  # Split effects into even and odd rows and
  # recombine into a single row without duplicate columns
  odd_effects <- effects[seq(1, nrow(effects), by = 2), ]
  even_effects <- effects[seq(2, nrow(effects), by = 2), ]
  effects <- merge(odd_effects, even_effects, by = "Marker_original")
  effects <- dplyr::mutate(
    effects,
    Trait = effects$Trait.x,
    Trait.x = NULL,
    Trait.y = NULL
  )

  # Merge stats and effects and return
  all_data <- merge(stats, effects, by = "Marker_original") %>%
    dplyr::mutate(
      Trait = .data$"Trait.x",
      Trait.x = NULL,
      Trait.y = NULL
    ) %>%
    dplyr::select(
      .data$Marker,
      .data$Marker_original,
      .data$Chr,
      .data$Pos,
      .data$p,
      .data$marker_R2,
      .data$Effect.x,
      .data$Effect.y
    )
  all_data
}




#' Load Linkage Disequilibrium
#'
#' @param LD_file The file containing linkage disequilibrium data
#' @param LD_columns The names of the columns in your linkage disequilibrium
#'   data for the chromosome of the first SNP, the position of the first SNP,
#'   the site of the first SNP, the chromosome of the second SNP, the position
#'   of the second SNP, the site of the second SNP, the distance between the
#'   two SNPs, and the R.2
#' @importFrom rlang .data
#' @importFrom stats complete.cases
#' @import dplyr
#' @return The linkage disequilibrium data in a list containing
#'   dataframes for each chromosome.
#' @export
#'
#' @examples
#' demo_LD_file = system.file("extdata","LD.txt.xz",
#'   package = "PAST", mustWork = TRUE)
#' LD <- load_LD(demo_LD_file)
load_LD <- function(LD_file,
                    LD_columns = c("Locus1",
                                   "Position1",
                                   "Site1",
                                   "Position2",
                                   "Site2",
                                   "Dist_bp",
                                   "R.2")) {

  LD <- read.table(LD_file, header = TRUE) %>%
    dplyr::mutate(Locus = as.character(!!as.name(LD_columns[1])),
                  Position1 = !!as.name(LD_columns[2]),
                  Site1 = !!as.name(LD_columns[3]),
                  Position2 = !!as.name(LD_columns[4]),
                  Site2 = !!as.name(LD_columns[5]),
                  Dist_bp = !!as.name(LD_columns[6]),
                  R.2 = as.numeric(!!as.name(LD_columns[7])),
                  Dist_bp = ifelse(.data$Dist_bp == "N/A",
                                   NA,
                                   .data$Dist_bp)) %>%
    dplyr::select(.data$Locus,
                  .data$Position1,
                  .data$Site1,
                  .data$Position2,
                  .data$Site2,
                  .data$Dist_bp,
                  .data$R.2)
  LD <- LD[complete.cases(LD), ]
  split(LD, f = LD$Locus)
}










#' Determine Linkage
#'
#' @param chunk A chunk of data to be processed
#' @param r_squared_cutoff The R^2 value to check against
#' @importFrom rlang .data
#' @import dplyr
#' @return Either the first unlinked SNP or a set of linked SNPs
determine_linkage <- function(chunk, r_squared_cutoff) {
  return_chunk = NULL
  for (i in seq_along(chunk)){
    
    #print(i)  
    data <- chunk[[i]]
    data <- data %>% dplyr::arrange(.data$Dist_bp)
    linked <- data %>% dplyr::filter(.data$R.2 >= r_squared_cutoff)
    unlinked <- data %>% dplyr::filter(.data$R.2 < r_squared_cutoff)
    
    if (nrow(unlinked) == nrow(data)) {
      return_chunk <- rbind(return_chunk, unlinked[1, ])
    } else {
      return_chunk <- rbind(return_chunk,
                            linked %>% filter(.data$R.2 >= r_squared_cutoff))
    }
  }
  return_chunk
}

#' Find representative SNP for a chunk of SNPs
#'
#' @param chunk A chunk of data to parse
#' @param r_squared_cutoff The R^2 value to check against when counting SNPs
#' @importFrom rlang .data
#' @import dplyr
#' @return A single SNP representing the whole chunk
find_representative_SNP <- function(chunk, r_squared_cutoff) {
pb()
  chunk <- chunk %>%
    dplyr::mutate(linked_snp_count = nrow(chunk),
                  linked_snp_count = ifelse(.data$R.2 < r_squared_cutoff,
                                            0,
                                            .data$linked_snp_count))
  # count the number of negative/positive effects in each chunk
  negative <- sum(chunk$SNP2_effect < 0)
  positive <- sum(chunk$SNP2_effect > 0)
  
  # Find SNP with largest negative or positive effect
  if (positive > negative) {
    # sort in descending order
    chunk <- chunk %>%
      dplyr::arrange(desc(.data$SNP2_effect), .data$Dist_bp)
    
    # get top row
    chunk <- chunk[1, ]
    
  } else if (negative > positive) {
    
    # sort in ascending order
    chunk <- chunk %>%
      dplyr::arrange(.data$SNP2_effect, .data$Dist_bp)
    
    # get top row
    chunk <- chunk[1, ]
    
  } else if (negative == positive) {
    if (chunk$SNP1_effect[[1]] > 0) {
      
      # sort in descending order
      chunk <- chunk %>%
        dplyr::arrange(desc(.data$SNP2_effect))
      
      # get top row
      chunk <- chunk[1, ] %>%
        dplyr::mutate(linked_snp_count = -1)
    } else{
      # sort in ascending order
      chunk <- chunk %>%
        dplyr::arrange(.data$SNP2_effect)
      
      # get top row
      chunk <- chunk[1, ] %>%
        dplyr::mutate(linked_snp_count = -1)
    }
  }
  chunk
  
}

#' Assign SNPs in a chunk to genes
#'
#' @param gff The GFF data for the chromosome being parsed
#' @param chunk The dataframe containing SNP data
#' @param window The search window around the SNPs
#' @importFrom rlang .data
#' @import dplyr
#' @return tagSNPs labeled with gene names
assign_chunk <- function(gff, chunk, window) {
  
  # set up data.frame of conditions to select tagSNP
  conditions <- chunk %>%
    dplyr::mutate(unlinked = .data$linked_snp_count == 0,
                  single_linked = .data$linked_snp_count == 1,
                  many_linked = .data$linked_snp_count > 1,
                  problem_linked = .data$linked_snp_count == -1,
                  effects_same_sign = .data$SNP1_effect * .data$SNP2_effect > 0,
                  largest_effect = ifelse(abs(.data$SNP1_effect) >
                                            abs(.data$SNP2_effect),
                                          "SNP1",
                                          "SNP2"),
                  effects_same_sign_use = ifelse(.data$SNP1_effect == .data$SNP2_effect,
                                                 "SNP2",
                                                 .data$largest_effect),
                  pvalue_equal = ifelse(.data$SNP1_pvalue == .data$SNP2_pvalue,
                                        TRUE,
                                        FALSE),
                  pvalue_not_equal_use = ifelse(.data$SNP1_pvalue > .data$SNP2_pvalue,
                                                "SNP2",
                                                "SNP1"),
                  single_pvalue_use = ifelse(.data$pvalue_equal == FALSE,
                                             .data$pvalue_not_equal_use,
                                             "PROBLEM"),
                  single_use = ifelse(.data$effects_same_sign == TRUE,
                                      .data$effects_same_sign_use,
                                      .data$pvalue_not_equal_use),
                  distance_use = ifelse(.data$Dist_bp < window,
                                        .data$effects_same_sign_use,
                                        "SNP2"),
                  many_use = ifelse(.data$effects_same_sign == TRUE,
                                    .data$distance_use,
                                    "SNP2"),
                  problem_use = ifelse(.data$effects_same_sign == TRUE,
                                       .data$effects_same_sign_use,
                                       "PROBLEM"),
                  chromosome = .data$Locus,
                  position = NA,
                  effect = NA,
                  p.value = NA,
                  distance = .data$Dist_bp) %>%
    dplyr::select(.data$unlinked,
                  .data$single_linked,
                  .data$many_linked,
                  .data$problem_linked,
                  .data$single_use,
                  .data$many_use,
                  .data$problem_use,
                  .data$Position1,
                  .data$Position2,
                  .data$SNP1_pvalue,
                  .data$SNP1_effect,
                  .data$SNP2_pvalue,
                  .data$SNP2_effect,
                  .data$chromosome,
                  .data$position,
                  .data$effect,
                  .data$p.value,
                  .data$distance,
                  .data$linked_snp_count,
                  .data$Marker_original
    )
  
  # handle unlinked SNPs
  conditions <- rbind(conditions %>% dplyr::filter(.data$unlinked == FALSE),
                      conditions %>% dplyr::filter(.data$unlinked == TRUE) %>%
                        dplyr::mutate(position = .data$Position1,
                                      effect = .data$SNP1_effect,
                                      p.value = .data$SNP1_pvalue)
  ) %>%
    dplyr::select(-.data$unlinked)
  
  # handle single linked SNPs
  conditions <- rbind(conditions %>%
                        dplyr::filter(.data$single_linked == FALSE),
                      conditions %>%
                        dplyr::filter(.data$single_linked == TRUE) %>%
                        dplyr::mutate(position = ifelse(.data$single_use == "SNP1",
                                                        .data$Position1,
                                                        .data$Position2),
                                      effect = ifelse(.data$single_use == "SNP1",
                                                      .data$SNP1_effect,
                                                      .data$SNP2_effect),
                                      p.value = ifelse(.data$single_use == "SNP1",
                                                       .data$SNP1_pvalue,
                                                       .data$SNP2_pvalue))
  ) %>%
    dplyr::select(-.data$single_linked,
                  -.data$single_use)
  
  # handle multiply-linked SNPs
  conditions <- rbind(conditions %>%
                        dplyr::filter(.data$many_linked == FALSE),
                      conditions %>%
                        dplyr::filter(.data$many_linked == TRUE) %>%
                        dplyr::mutate(position = ifelse(.data$many_use == "SNP1",
                                                        .data$Position1,
                                                        .data$Position2),
                                      effect = ifelse(.data$many_use == "SNP1",
                                                      .data$SNP1_effect,
                                                      .data$SNP2_effect),
                                      p.value = ifelse(.data$many_use == "SNP1",
                                                       .data$SNP1_pvalue,
                                                       .data$SNP2_pvalue))
  ) %>%
    dplyr::select(-.data$many_linked,
                  -.data$many_use)
  
  # handle problem SNPs
  tagSNPs <- rbind(conditions %>%
                     dplyr::filter(.data$problem_linked == FALSE),
                   conditions %>%
                     dplyr::filter(.data$problem_linked == TRUE) %>%
                     dplyr::mutate(position = ifelse(.data$problem_use == "SNP1",
                                                     .data$Position1,
                                                     .data$Position2),
                                   effect = ifelse(.data$problem_use == "SNP1",
                                                   .data$SNP1_effect,
                                                   .data$SNP2_effect),
                                   p.value = ifelse(.data$problem_use == "SNP1",
                                                    .data$SNP1_pvalue,
                                                    .data$SNP2_pvalue))
  ) %>%
    dplyr::mutate(marker = .data$Marker_original) %>%
    dplyr::select(.data$chromosome,
                  .data$position,
                  .data$marker,
                  .data$effect,
                  .data$p.value,
                  .data$distance,
                  .data$linked_snp_count) %>%
    dplyr::mutate(window_start = .data$position - window,
                  window_end = .data$position + window) %>%
    dplyr::arrange(.data$position)
  
  # assigned_genes = tagSNPs %>%
  #   rowwise() %>%
  #   mutate(name = ifelse(length(which((.data$window_start >= gff$start &
  #                                 .data$window_end <= gff$end) |
  #                                (gff$start >= .data$window_start &
  #                                   gff$end <= .data$window_end) |
  #                                (.data$window_start >= gff$start &
  #                                   .data$window_start <= gff$end) |
  #                                (.data$window_end >= gff$start &
  #                                   .data$window_end <= gff$end)) > 0),
  #                        gff[which((.data$window_start >= gff$start &
  #                                     .data$window_end <= gff$end) |
  #                                    (gff$start >= .data$window_start &
  #                                       gff$end <= .data$window_end) |
  #                                    (.data$window_start >= gff$start &
  #                                       .data$window_start <= gff$end) |
  #                                    (.data$window_end >= gff$start &
  #                                       .data$window_end <= gff$end)), 4],
  #                        NA))
  #
  # as.data.frame(assigned_genes %>%
  #                 dplyr::select(chromosome,
  #                               position,
  #                               effect,
  #                               p.value,
  #                               linked_snp_count,
  #                               name)) %>%
  #   dplyr::filter(is.na(.data$name) == FALSE)
  
  # assign genes
  inner_join(tagSNPs, gff, by = c("chromosome" = "seqid")) %>%
    dplyr::filter(
      (
        .data$window_start >= .data$start &
          .data$window_end <= .data$end
      ) |
        (
          .data$start >= .data$window_start & .data$end <= .data$window_end
        ) |
        (
          .data$window_start >= .data$start &
            .data$window_start <= .data$end
        ) |
        (
          .data$window_end >= .data$start &
            .data$window_end <= .data$end
        )
    ) %>%
    dplyr::select(-c(
      .data$distance,
      .data$window_start,
      .data$window_end,
      .data$start,
      .data$end
    )) %>%
    dplyr::mutate(name = .data$Name,
                  Name = NULL)
}

#' Find the SNP-gene assignment that represents SNPs assigned to a gene
#'
#' @param chunk A chunk of gene assignments
#' @importFrom rlang .data
#' @import dplyr
#' @return A single SNP-gene assignment representing all SNPS assigned to the
#'   same gene
#' to a gene
find_representative_SNP_gene_pairing <- function(chunk) {
  neg_genes <- chunk %>%
    dplyr::filter(.data$effect < 0) %>%
    dplyr::arrange(.data$effect, .data$p.value) %>%
    dplyr::mutate(linked_snp_count = ifelse(.data$linked_snp_count <= 0,
                                            1,
                                            .data$linked_snp_count))
  pos_genes <- chunk %>%
    dplyr::filter(.data$effect > 0) %>%
    dplyr::arrange(desc(.data$effect), .data$p.value) %>%
    dplyr::mutate(linked_snp_count = ifelse(.data$linked_snp_count <= 0,
                                            1,
                                            .data$linked_snp_count))
  negative <- sum(neg_genes$linked_snp_count)
  positive <- sum(pos_genes$linked_snp_count)
  
  if (positive > negative) {
    chunk <- chunk %>%
      dplyr::arrange(desc(.data$effect), desc(.data$p.value))
    representative <- chunk[1, ]
  } else if (negative > positive) {
    chunk <- chunk %>%
      dplyr::arrange(.data$effect, .data$p.value)
    representative <- chunk[1, ]
  } else if (positive == negative) {
    pos_max <- pos_genes[1, ]
    neg_max <- neg_genes[1, ]
    
    if (pos_max$effect > abs(neg_max$effect)) {
      representative <- pos_max
    }
    else{
      representative <- neg_max
    }
  }
  representative %>% dplyr::mutate(linked_snp_count = negative + positive)
}

#' Assign SNPs to genes
#'
#' @param gwas_data Merged association and effects data from merge_data()
#' @param LD Linkage disequilibrium data from parse_LD()
#' @param gff_file The path to a GFF file
#' @param window The search window for genes around the SNP
#' @param r_squared_cutoff The R^2 value used to determine SNP significance
#' @param num_cores The number of cores to use in parallelizing PAST
#' @importFrom rlang .data
#' @importFrom rtracklayer readGFF
#' @import dplyr
#' @import foreach
#' @import doParallel
#' @return A dataframe of genes from the SNP data
#' @export
#'
#' @examples
#' example("load_GWAS_data")
#' example("load_LD")
#' demo_genes_file = system.file("extdata", "genes.gff",
#'   package = "PAST", mustWork = TRUE)
#' genes <-assign_SNPs_to_genes(gwas_data, LD, demo_genes_file, 1000, 0.8, 2)

###
assign_SNPs_to_genes <-
  function(gwas_data,
           LD,
           gff_file,
           window,
           r_squared_cutoff,
           num_cores) {
    
    full_gff <- rtracklayer::readGFF(gff_file) %>%
      filter(.data$type == "gene") %>%
      select(.data$seqid,
             .data$start,
             .data$end,
             .data$Name)
    chromosomes <- as.character.factor(
      full_gff %>%
        dplyr::select(.data$seqid) %>%
        dplyr::arrange(.data$seqid) %>%
        unique() %>%
        dplyr::pull(.data$seqid)
    )
    
    cl <- parallel::makeCluster(num_cores,type = "FORK")
    registerDoSNOW(cl)
    print(paste("Register Cluster of ",num_cores," cores",sep="") )
    
    all_genes <- NULL
    
    # UP/DOWNSTREAM LOOP
    for (n in seq_along(c(1, 2))) {
      
      # BEGIN PROCESSING BY CHROMOSOMES LOOP
      for (chromosome in names(LD)) {
        print(paste("procrssing chr",chromosome,sep=""))
        
        if (chromosome %in% chromosomes) {
          temp_data <-
            LD[[chromosome]] %>%
            dplyr::mutate(Marker1 = paste0(.data$Locus,
                                           "_",
                                           .data$Position1),
                          Marker2 = paste0(.data$Locus,
                                           "_",
                                           .data$Position2)) %>%
            dplyr::arrange(.data$Position1)
          
          if (n == 2) {
            temp_data <-
              temp_data %>%
              dplyr::mutate(
                temp_p2 = .data$Position1,
                temp_s2 = .data$Site1,
                temp_Marker2 = .data$Marker1,
                Position1 = .data$Position2,
                Site1 = .data$Site2,
                Marker1 = .data$Marker2
              ) %>%
              dplyr::mutate(
                Site2 = .data$temp_s2,
                Position2 = .data$temp_p2,
                Marker2 = .data$temp_Marker2,
                temp_s2 = NULL,
                temp_p2 = NULL,
                temp_Marker2 = NULL
              ) %>%
              dplyr::arrange(.data$Site1)
          }
          split <- 1000
          temp_data_list <- split(temp_data, temp_data$Position1)
          temp_data_list_split <- split(temp_data_list,
                                        ceiling(seq_along(temp_data_list)
                                                /split))
          
          rm(temp_data)
          rm(temp_data_list)  
          print("DL")    
          t1<-proc.time()      
          temp_data <-
            foreach(
              data = temp_data_list_split,
              .combine = rbind,
              .packages = "dplyr",
              .export =  c("determine_linkage")
            ) %dopar% {
              determine_linkage(data, r_squared_cutoff)
              
            }
          rm(temp_data_list_split)
          proc.time()-t1
          print("DL done")
          gwas_data_for_chromosome <-
            dplyr::filter(gwas_data, .data$Chr == as.integer(chromosome))
          
          print("look up p-value and effect data for SNP1")
          temp_data <-
            merge(temp_data, gwas_data_for_chromosome,
                  by.x = "Marker1",
                  by.y = "Marker") %>%
            dplyr::mutate(SNP1_pvalue = .data$p,
                          SNP1_effect = .data$Effect.x) %>%
            dplyr::select(
              .data$Locus,
              .data$Position1,
              .data$Position2,
              .data$Site1,
              .data$Site2,
              .data$Dist_bp,
              .data$R.2,
              .data$Marker1,
              .data$SNP1_pvalue,
              .data$SNP1_effect,
              .data$Marker2,
              .data$Marker_original
            )
          
          print("look up p-value and effect data for SNP2")
          temp_data <-
            merge(temp_data, gwas_data_for_chromosome,
                  by.x = "Marker2",
                  by.y = "Marker") %>%
            dplyr::mutate(SNP2_pvalue = .data$p,
                          SNP2_effect = .data$Effect.x,
                          Marker_original = Marker_original.x) %>%
            dplyr::select(
              .data$Locus,
              .data$Position1,
              .data$Position2,
              .data$Site1,
              .data$Site2,
              .data$Dist_bp,
              .data$R.2,
              .data$Marker1,
              .data$SNP1_pvalue,
              .data$SNP1_effect,
              .data$Marker2,
              .data$SNP2_pvalue,
              .data$SNP2_effect,
              .data$Marker_original
            ) %>%
            dplyr::arrange(.data$Position1)
          
          
          singles = temp_data %>% dplyr::group_by(.data$Marker1) %>%
            dplyr::summarise(count = n()) %>%
            dplyr::filter(count == 1)
          
          singles <- temp_data %>% filter(.data$Marker1 %in% singles$Marker1) %>%
            mutate(linked_snp_count = 1)
          blocks <- temp_data %>% filter(!(.data$Marker1 %in% singles$Marker1))
          index <- c(0, cumsum(abs(diff(blocks$Site2)) > 1))
          rm(temp_data)
          #print(ls())
           save(index,blocks,singles,file=paste("blocks_chr",chromosome,"_",n,".RData",sep=""))
          print("Now going hard")
          
          temp_data_list <- split(blocks, paste(blocks$Position1, index))
          print(length(temp_data_list))
          rm(blocks)
          temp_data<-NULL #
          temp_data_list<-head(temp_data_list,5000) #
          ####doing in sequencial
          temp_data <-
            foreach(
              data = temp_data_list,
              .combine = rbind,
              .packages = "dplyr"

            ) %do% {
           
              find_representative_SNP(data, r_squared_cutoff)

              
            }

####doing in doSNOW 
#pb <- txtProgressBar(max = length(temp_data_list), style = 3)
#progress <- function(n) setTxtProgressBar(pb, n)
#opts <- list(progress = progress)
#temp_data <- foreach(d=iter(temp_data_list, by = "row"),
#                    .combine = rbind, 
#                        .packages = "dplyr",
#                        .export =  c("find_representative_SNP"),
#                  .options.snow = opts) %dopar%
#{
# find_representative_SNP(d, r_squared_cutoff)
#}
#close(pb)



          rm(temp_data_list)
          temp_data <- rbind(temp_data, singles)
          rm(singles)
          split <- 1000
          snp_list <-
            split(temp_data,
                  rep(seq_len(split),
                      length.out = nrow(temp_data),
                      each = ceiling(nrow(temp_data) / split)
                  ))
          rm(temp_data)
           gc()
          # subset gff to only handle this chromosome
          gff <- dplyr::filter(full_gff, .data$seqid == chromosome) %>%
            mutate(seqid = as.character.factor(.data$seqid))
          print("assign_chunk")
          print(length(snp_list))
          # get genes in parallel
          chromosome_genes <-
            foreach(
              
              snp_chunk = snp_list,
              .combine = rbind,
              .packages = "dplyr",
              .export =  c("assign_chunk")
            ) %dopar% {
              assign_chunk(gff, snp_chunk, window)
            }
          
          all_genes <- rbind(all_genes, chromosome_genes)
          if(n==2){
          write.csv(all_genes,file=paste("all_chr_",chromosome,".csv",sep=""),row.names = FALSE)
          }
          rm(snp_list)
        }
      }
    }
      print("ass_gene")
    
    group_by_gene <- split(all_genes, f = all_genes$name)
pb <- txtProgressBar(max = length(group_by_gene), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
representative_genes <-
      foreach(
        chunk = group_by_gene,
        .combine = rbind,
        .packages = c("dplyr", "PAST"),
		.options.snow = opts
      ) %dopar% {
        find_representative_SNP_gene_pairing(chunk)
      }
       close(pb)
       write.csv(representative_genes,file=paste("repr_genes_chr",chr,".csv",sep=""),row.names = FALSE)
       save(representative_genes,file=paste("repr_genes_chr",chr,".RData",sep=""))
  gc()
  
    stopCluster(cl)
    
    representative_genes
  }
##############################


get_factors <- function(gene_ranks) {
  factors <- vector("logical", length(gene_ranks))
  factors[1] <- gene_ranks[1] - 1
  for (i in 2:length(gene_ranks)) {
    factors[i] <- gene_ranks[i] - gene_ranks[i - 1] - 1 + factors[i - 1]
  }
  factors
}

get_pmisses <- function(temp_data, genes_in_pathway, factors) {
  pmisses <- factors / (nrow(temp_data) - nrow(genes_in_pathway))
  pmisses
}

get_phits <- function(gene_effects) {
  NR <- sum(abs(gene_effects))
  absolute <- abs(gene_effects)
  phits <- vector("logical", length(gene_effects))
  phits[1] <- (absolute[1] / NR)
  for (i in 2:length(gene_effects)) {
    phits[i] <- absolute[i] / NR + phits[i - 1]
  }
  phits
}

get_running_enrichment_score <- function(phits, pmisses) {
  running_enrichment_score <- phits - pmisses
  running_enrichment_score
}

find_max <- function(running_enrichment_score) {
  running_enrichment_score <- sort(running_enrichment_score, decreasing = TRUE)
  max <- running_enrichment_score[[1]]
  max
}

# This isn't really global, but R CMD check
# doesn't understand how to deal with variable
# created by foreach loops
globalVariables("pathway")

#' Find Pathway Significance
#'
#' @param genes Genes from assign_SNPs_to_genes()
#' @param pathways_file A file containing the pathway IDs, their names, and the
#'   genes in the pathway
#' @param gene_number_cutoff A cut-off for the minimum number of genes in a
#'   pathway
#' @param mode increasing/decreasing
#' @param sample_size How many times to sample the effects data during random
#'   sampling
#' @param num_cores The number of cores to use in parallelizing PAST
#' @importFrom rlang .data
#' @importFrom stats pnorm
#' @importFrom stats sd
#' @import dplyr
#' @import foreach
#' @import doParallel
#' @import iterators
#' @import parallel
#' @import qvalue
#' @return Rugplots data
#' @export
#'
#' @examples
#' example("assign_SNPs_to_genes")
#' demo_pathways_file = system.file("extdata", "pathways.txt.xz",
#'   package = "PAST", mustWork = TRUE)
#' rugplots_data <- find_pathway_significance(genes, demo_pathways_file, 5,
#'   "increasing", 1000, 2)
find_pathway_significance <-
  function(genes,
           pathways_file,
           gene_number_cutoff = 5,
           mode,
           sample_size = 1000,
           num_cores) {
    # load pathways
    pathways <-
      read.table(pathways_file,
                 sep = "\t",
                 header = TRUE,
                 quote = "")
    
    # sample to create 1000 random distributions
    effects <- genes %>% dplyr::select(.data$name, .data$effect)
    effects <-
      cbind(effects, vapply(seq_len(sample_size),
                            function(i) sample(effects$effect),
                            FUN.VALUE = double(nrow(effects))))
    
    pathways_unique <- unique(select(pathways, .data$pathway_id))
    pathways_unique[] <- lapply(pathways_unique, as.character)
    print("makeCluster")
    # process sample columns
    cl <- parallel::makeCluster(num_cores, outfile = "debug.txt")
    #registerDoParallel(cl)
    
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = 1000, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    column_observations <-
      foreach(
        i = iter(2:(sample_size + 2)),
        .combine = cbind,
        .packages = c("dplyr", "PAST", "foreach", "iterators","doSNOW"),
        .export=c("get_factors","get_pmisses","get_phits","get_running_enrichment_score","find_max"),
        .options.snow = opts
      ) %dopar% {
       
        temp_data <-
          data.frame(matrix("NA", ncol = 2, nrow = nrow(effects)))
        colnames(temp_data) <- c("gene", "effect")
        temp_data$effect <- effects[, i]
        temp_data$gene <- effects[, 1]
        
        if (mode == "decreasing") {
          temp_data <- temp_data %>%
            dplyr::arrange(.data$effect) %>%
            dplyr::mutate(rank = row_number(),
                          effect = abs(.data$effect))
        } else if (mode == "increasing") {
          temp_data <- temp_data %>%
            dplyr::arrange(desc(.data$effect)) %>%
            dplyr::mutate(rank = row_number())
        } else {
          stop("Incorrect mode.")
        }
        
        column_observation <- data.frame()
        
        foreach(pathway = iter(pathways_unique$pathway_id, by = "row")) %do% {
          genes_in_pathway <-
            dplyr::filter(pathways, pathways$pathway_id == pathway) %>%
            mutate(gene_id = as.character.factor(gene_id))
          
          ## get ranks and effects and sort by rank
          genes_in_pathway <-
            merge(genes_in_pathway,
                  temp_data,
                  by.x = "gene_id",
                  by.y = "gene") %>%
            dplyr::arrange(.data$rank) %>% unique()
          
          # check cutoff
          if (nrow(genes_in_pathway) >= gene_number_cutoff) {
            
            # get factors using rank
            factors <- get_factors(genes_in_pathway$rank)
            
            # get pmisses
            pmisses <- get_pmisses(temp_data, genes_in_pathway, factors)
            
            # get phits
            phits <- get_phits(genes_in_pathway$effect)
            
            # get phits-pmisses
            running_enrichment_score <- get_running_enrichment_score(phits, pmisses)
            find_max(running_enrichment_score)
            # store max phit_pmisses
            column_observation <-
              rbind(column_observation, find_max(running_enrichment_score))
          } else {
            column_observation <- rbind(column_observation, NA)
          }
        }
        column_observation
      }
    
    close(pb)
    #print("hear")
    stopCluster(cl)
    
    pathways_unique <- cbind(pathways_unique, column_observations)
    colnames(pathways_unique) <-
      c("Pathway", "ES_Observed", seq_len(sample_size))
    colnames(pathways_unique)[3:(sample_size + 2)] <-
      paste0("ES", colnames(pathways_unique)[3:(sample_size + 2)])
    pathways_unique <- pathways_unique %>%
      dplyr::filter(!is.na(.data$ES_Observed))
    pathways_unique <- pathways_unique %>%
      dplyr::mutate(permutation_mean =
                      apply(pathways_unique[, 3:(sample_size + 2)], 1, mean))
    pathways_unique <- pathways_unique %>%
      dplyr::mutate(permutation_standard_deviation =
                      apply(pathways_unique[, 3:(sample_size + 2)], 1, sd))
    pathways_unique <- pathways_unique %>%
      dplyr::mutate(
        NES_Observed = (.data$ES_Observed - .data$permutation_mean) / .data$permutation_standard_deviation
      )
    for (i in seq_along(sample_size)) {
      pathways_unique[i + 2] <-
        (pathways_unique[i + 2] - pathways_unique[sample_size + 3]) / pathways_unique[sample_size + 4]
    }
    pathways_unique <- pathways_unique %>%
      dplyr::mutate(pvalue = 1-pnorm(.data$NES_Observed))
    pathways_unique <-dplyr::arrange(pathways_unique, .data$pvalue)
    pathways_unique <- pathways_unique %>%
      dplyr::select(.data$Pathway,
                    .data$ES_Observed,
                    .data$NES_Observed,
                    .data$pvalue) %>%
      dplyr::mutate(qvalue = qvalue(
        pathways_unique$pvalue,
        lambda = 0,
        fdr.level = 0.05
      )$qvalues)
    
    pathways_significant <- pathways_unique %>%
      dplyr::select(.data$Pathway,
                    .data$NES_Observed,
                    .data$pvalue,
                    .data$qvalue) %>%
      dplyr::mutate(pathway_number = row_number(), NESrank = NULL)
    
    temp_data <-
      data.frame(matrix("NA", ncol = 2, nrow = nrow(effects)))
    colnames(temp_data) <- c("Gene", "Effect")
    temp_data$Effect <- effects[, 2]
    temp_data$Gene <- effects[, 1]
    colnames(temp_data) <- c("Gene", "Effect")
    if (mode == "decreasing") {
      temp_data <- temp_data %>%
        dplyr::arrange(.data$Effect) %>%
        dplyr::mutate(rank = row_number(),
                      effect = abs(.data$Effect))
    } else if (mode == "increasing") {
      temp_data <- temp_data %>%
        dplyr::arrange(desc(.data$Effect)) %>%
        dplyr::mutate(rank = row_number())
    } else {
      stop("Incorrect mode.")
    }
    
    rugplots_data <- NULL
    rugplots_data <-
      foreach(
        pathway = iter(pathways_significant$Pathway, by = "row"),
        .combine = rbind
      ) %do% {
        genes_in_pathway <-
          dplyr::filter(pathways, pathways$pathway_id == pathway)
        
        ## get ranks and effects and sort by rank
        genes_in_pathway <-
          merge(genes_in_pathway,
                temp_data,
                by.x = "gene_id",
                by.y = "Gene") %>%
          dplyr::arrange(rank)
        
        # get factors using rank
        genes_in_pathway$factors <-
          get_factors(genes_in_pathway$rank)
        
        # get pmisses
        genes_in_pathway$pmisses <-
          get_pmisses(temp_data, genes_in_pathway, genes_in_pathway$factors)
        
        # get phits
        genes_in_pathway$phits <- get_phits(genes_in_pathway$Effect)
        
        # get phits-pmisses
        genes_in_pathway$running_enrichment_score <-
          get_running_enrichment_score(genes_in_pathway$phits, genes_in_pathway$pmisses)
        
        # append rows with NESrank
        genes_in_pathway <-
          merge(genes_in_pathway,
                pathways_significant,
                by.x = "pathway_id",
                by.y = "Pathway")
        genes_in_pathway
      }
    rugplots_data <- merge(
      rugplots_data %>%
        dplyr::select(
          .data$pathway_id,
          .data$pathway_number,
          .data$gene_id,
          .data$rank,
          .data$running_enrichment_score
        ),
      pathways %>%
        dplyr::select(.data$pathway_id, .data$pathway_name) %>%
        unique(),
      by = "pathway_id"
    )
    rugplots_data <- merge(
      rugplots_data,
      dplyr::select(
        pathways_unique,
        .data$Pathway,
        .data$pvalue,
        .data$qvalue
      ),
      by.x = "pathway_id",
      by.y = "Pathway"
    ) %>%
      dplyr::arrange(.data$pathway_number)
    rugplots_data
  }
plot_pathways <-
  function(rugplots_data,
           filter_type,
           filter_parameter,
           mode,
           output_directory) {

    rugplots_data <- rugplots_data %>%
      dplyr::arrange(.data$pathway_number)
    write.table(
      rugplots_data,
      file = paste0(output_directory, "/", mode, ".full.txt"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )

    if (filter_type == "rank") {
      rugplots_data <- rugplots_data %>%
        dplyr::filter(.data$pathway_number <= filter_parameter)
    } else if (filter_type == "pvalue") {
      rugplots_data <- rugplots_data %>%
        dplyr::filter(.data$pvalue <= filter_parameter)
    } else {
      print("Incorrect filtering type. Filtering at p-value <= 0.05")
      rugplots_data <- rugplots_data %>%
        dplyr::filter(.data$pvalue <= 0.05)
    }

    write.table(
      rugplots_data,
      file = paste0(output_directory, "/", mode, ".filtered.txt"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    rugplots_split <- split(rugplots_data, rugplots_data$pathway_number)

    for (rank in names(rugplots_split)) {
      temp_data <- rugplots_split[[rank]]
      title <-
        paste0(unique(as.character(temp_data$pathway_id)), " - ",
               unique(as.character(temp_data$pathway_name)))
      intercept <- temp_data %>%
        dplyr::arrange(desc(.data$running_enrichment_score)) %>%
        dplyr::select(.data$rank)
      intercept <- intercept[, 1][1]
      rugplot <-
        ggplot(temp_data,
               aes(x = temp_data$rank,
                   y = temp_data$running_enrichment_score)) +
        geom_line(stat = "identity") +
        geom_rug(sides = "t", position = "jitter") +
        geom_vline(xintercept = intercept,
                   color = "black",
                   linetype = "longdash") +
        ggtitle(title) +
        labs(x = "Gene Rank", y = "Running Enrichment Score") +
        scale_x_continuous(breaks = c(0, 5000, 10000, 15000, 20000, 25000)) +
        theme(
          axis.text = element_text (color = "black"),
          panel.background = element_rect (color = "black", fill = "pink")
        )
      ggsave(paste0(
        output_directory,
        "/",
        mode,
        ".",
        unique(as.character(temp_data$pathway_id)),
        ".png"
      ),
      rugplot)
    }
  }


