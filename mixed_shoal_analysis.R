#### Libraries ####
rm(list=ls())

library(tidyverse)
library(magrittr)
library(rentrez)
library(msa)
library(ape)
library(strataG)
library(adegenet)
library(starmie)
library(janitor)
library(posterior)
library(limma)
library(ggforce)

#
#### Read in Data ####
full_cope_data <- read_csv('Raw Data.csv')

#### Morphology Assignment ####
full_cope_data <- full_cope_data %>%
  mutate(morph_species = case_when(is.na(AIC_pores) ~ 'unknown',
                                   AIC_pores == 1 ~ 'cpers',
                                   AIC_pores == 2 ~ 'chya',
                                   TRUE ~ 'error'))

#### CO1 Processing ####
## Download Baldwin et al. Sequences
NCBI_sequence_data <- entrez_search(db = "nucleotide",
                                    term = '(Coryphopterus personatus[Organism] OR Coryphopterus hyalinus[Organism]) AND COI[Gene] AND 
                                    Reconciling Genetic Lineages with Species in Western Atlantic Coryphopterus (Teleostei: Gobiidae)',
                                    retmax=9999)$ids %>%
  entrez_fetch(db = 'nucleotide', id = ., rettype = 'fasta') %>%
  str_split(pattern = '>') %>% 
  unlist %>%
  str_subset('mitochondrial') %>%
  str_remove_all('\n') %>%
  str_remove('mitochondrial') %>%
  tibble(COI_sequence = .) %>%
  separate(COI_sequence, into = c('info', 'COI_sequence'), sep = ';') %>%
  mutate(COI_sequence = str_trim(COI_sequence)) %>%
  separate(info, into = c('genbankID_nuc','info'), sep = ' ', extra = 'merge') %>%
  mutate(species = case_when(
    str_detect(info, 'personatus') ~ 'Cpers',
    str_detect(info, 'lipernes') ~ 'Clip',
    str_detect(info, 'hyalinus') ~ 'Chya'
  )) %>%
  filter(str_detect(species, 'Clip', negate = TRUE)) %>%
  filter(str_detect(info, 'UNVERIFIED', negate = TRUE)) %>%
  mutate(ID = str_c(species,1:n(), sep='_')) %>%
  dplyr::select(ID, COI_sequence, genbankID_nuc)

## Align with our sequences
aligned_sequence <- NCBI_sequence_data %>%
  bind_rows(full_cope_data, .id = 'source') %>%
  mutate(source = case_when(source == 1 ~ 'NCBI',
                            TRUE ~ 'Collection')) %>%
  filter(!is.na(COI_sequence)) %>%
  dplyr::select(source, ID, COI_sequence) %>%
  mutate(aligned = msa(COI_sequence, type='dna', verbose = TRUE, order = 'input') %>%
           as.character)

## Trim to core region 
start_end_bases <- aligned_sequence %>%
  filter(source != 'NCBI') %>%
  select(-source) %>%
  mutate(aligned_length = map_int(aligned, str_length)) %>%
  mutate(align_trim = map_chr(aligned, str_remove_all, pattern = c('^-+|-+$'))) %>%
  mutate(trimmed_length = map_int(align_trim, str_length)) %>%
  mutate(start_dash = map2_int(aligned, aligned_length, function(x, y) y - str_length(str_remove_all(x, '^-+'))),
         end_dash = map2_int(aligned, aligned_length, function(x, y) y - str_length(str_remove_all(x, '-+$')))) %>%
  
  filter(start_dash < 60,
         end_dash < 100) %>%
  summarise(start = max(start_dash), end = max(end_dash))

trimmed_sequence <- aligned_sequence %>%
  mutate(aligned = map_chr(aligned, ~str_sub(.x, start = start_end_bases$start + 1L, end = str_length(.x) - start_end_bases$end))) %>%
  mutate(trim_len = str_length(aligned)) %>%
  mutate(start_dash = map2_int(aligned, trim_len, function(x, y) y - str_length(str_remove_all(x, '^-+'))),
         end_dash = map2_int(aligned, trim_len, function(x, y) y - str_length(str_remove_all(x, '-+$')))) %>%
  filter(start_dash == 0,
         end_dash == 0) %>%
  select(source, ID, aligned)

## Find diagnostic SNPs
## Based on only the NCBI fastas
per_position_base <- trimmed_sequence %>%
  mutate(bases = str_split(aligned, '')) %>%
  unnest(bases) %>%
  dplyr::select(-aligned) %>%
  group_by(ID) %>%
  mutate(position = 1:n()) %>%
  ungroup 

diagnostic_snps <- per_position_base %>%
  filter(source == 'NCBI') %>%
  dplyr::select(ID, bases, position) %>%
  separate(ID, into = c('species'), sep = '_', remove = FALSE, extra = 'drop') %>%
  group_by(species, position) %>%
  mutate(number_bases = n_distinct(bases)) %>%
  filter(number_bases == 1) %>%
  ungroup %>%
  group_by(position) %>% 
  mutate(number_bases = n_distinct(bases)) %>%
  filter(number_bases > 1) %>%
  group_by(species, position, bases) %>%
  summarise(.groups = 'drop') %>%
  
  ungroup %>%
  filter(bases != '-', bases != 'N') %>%
  group_by(position) %>%
  filter(n() == 2) %>%
  ungroup %>%
  spread(species, bases)

## Use diagnostic SNPs to ID fish ##
## Prior is entirely weighted to either 0 (COHY) or 1 (COPE)
a<-0.5; b<-0.5 #Jeffreys Prior - expect to be either one or the other

identification_probability <- per_position_base %>%
  filter(source != 'NCBI') %>%
  inner_join(diagnostic_snps, ., by = 'position') %>%
  filter(bases != '-') %>%
  dplyr::select(-source) %>%
  
  mutate(match = case_when(bases == Chya ~ 'chya',
                           bases == Cpers ~ 'cpers',
                           TRUE ~ 'neither')) %>%
  
  group_by(ID) %>%
  summarise(n_chya = sum(match == 'chya'),
            n_cpers = sum(match == 'cpers'),
            total = n()) %>%
  gather(co1_species, number, -ID, -total) %>%
  mutate(co1_species = str_remove(co1_species, 'n_')) %>%
  
  mutate(prob_id = (number+a)/(total+a+b),
         id_lwr.95 = qbeta(0.025, number+a, total-number+b),
         id_lwr.50 = qbeta(0.25, number+a, total-number+b),
         id_upr.50 = qbeta(0.75, number+a, total-number+b),
         id_upr.95 = qbeta(0.975, number+a, total-number+b)) %>%
  dplyr::select(-number, -total) %>%
  group_by(ID) %>%
  filter(prob_id == max(prob_id)) %>%
  mutate(co1_species = case_when(id_lwr.95 > 0.5 ~ co1_species, 
                                 TRUE ~ 'unknown')) %>%
  ungroup

identification_probability %>%
  ggplot(aes(x = ID, y = prob_id, colour = co1_species)) +
  geom_linerange(aes(ymin = id_lwr.95, ymax = id_upr.95), lty = 1, size = 1) +
  geom_linerange(aes(ymin = id_lwr.50, ymax = id_upr.50), size = 3) +
  #geom_point(size = 6) +
  scale_y_continuous(limits = c(0,1)) +
  scale_color_manual(values = c('cpers'='#F8766D', 'chya'='#00BFC4', 'unknown'='black')) +
  theme_classic() +
  ylab("Identification Probability") +
  xlab('Individual') +
  theme(axis.title.y = element_text(colour = 'black', size = 36),
        axis.text.y = element_text(colour = 'black', size = 30),
        axis.title.x = element_text(colour = 'black', size = 36),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "transparent", colour = NA), 
        panel.border = element_rect(fill = NA, colour = 'black'), 
        plot.background = element_rect(size = 1, linetype = 'solid', colour = 'black', fill = "transparent"),
        legend.position='none')

## Set CO1 ID
full_cope_data <- full_cope_data %>%
  left_join(identification_probability %>%
              dplyr::select(ID, co1_species), by = 'ID') %>%
  mutate(co1_species = if_else(is.na(co1_species), 'unknown', co1_species))

## Make Haplotype network 
haplotype_data <- trimmed_sequence %>%
  dplyr::select(source, ID, aligned) %>%
  full_join(full_cope_data %>%
              dplyr::select(ID, co1_species, AIC_pores),
            by = 'ID') %>% 
  filter(!is.na(aligned)) %>% 
  mutate(species = case_when(is.na(co1_species) ~ str_remove(ID, '_.*'),
                             TRUE ~ str_to_sentence(co1_species))) %>%
  dplyr::select(-co1_species) %>%
  
  mutate(source = if_else(source == 'NCBI', source, str_c('z',source,sep='')),
         label_1 = str_c(species, source, sep = '-'),
         label_2 = case_when((is.na(AIC_pores) & source == 'NCBI') ~ species,
                             AIC_pores == 1 ~ 'one',
                             AIC_pores == 2 ~ 'two',
                             TRUE ~ 'unknown'))

write_popart <- function(x, file){
  popart_data_header <- "\nBEGIN TRAITS;\n  Dimensions NTRAITS=5;\n  Format labels=yes missing=? separator=Comma;\n  TraitLabels Chya Cpers one two unknown;\n  Matrix\n"
  
  tmp <- x %$%
    str_to_lower(aligned) %>%
    str_split('', simplify = TRUE) %>%
    set_rownames(x$ID) %>%
    as.DNAbin %T>%
    write.nexus.data(file)
  
  popart_data <- x %>%
    select(ID, label_2) %>%
    mutate(number = 1) %>%
    spread(label_2, number, fill = 0) %>%
    arrange(ID) %>%
    mutate(for_popart = str_c(ID, str_c(Chya, Cpers, one, two, unknown, sep = ','), sep = ' '))
  
  str_c(colnames(popart_data)[2:6], collapse = ' ')
  
  write_lines(x = popart_data_header, file, append = TRUE)
  write_lines(x = popart_data$for_popart, file, append = TRUE)
  
}

write_popart(haplotype_data, "haplotypes.nex")
## Go to popart to make the hap network ##

## Calculate nucleotide divergence
the_meta <- haplotype_data %>%
  filter(source != 'NCBI') %>%
  select(-source, -AIC_pores, -label_1, -label_2) %>%
  select(-aligned) %>%
  as.data.frame %>%
  set_rownames(.$ID)

the_seq <- haplotype_data %>%
  filter(source != 'NCBI') %>%
  select(-source, -AIC_pores, -label_1, -label_2) %>%
  mutate(aligned = str_to_lower(aligned),
         aligned = map(aligned, str_split, pattern = ''),
         aligned = map(aligned, unlist)) %$%
  aligned %>%
  set_names(the_meta$ID)

df2gtypes(the_meta[,c('ID', 'species', 'ID')], ploidy = 1, sequences = the_seq) %>% nucleotideDivergence


#### Microsat Processing ####
genind2structure <- function(obj, file="", pops=FALSE){
  #https://github.com/lvclark/R_genetics_conv/blob/master/genind2structure.R
  
  if(!"genind" %in% class(obj)){
    warning("Function was designed for genind objects.")
  }
  
  # get the max ploidy of the dataset
  pl <- max(obj@ploidy)
  # get the number of individuals
  S <- adegenet::nInd(obj)
  # column of individual names to write; set up data.frame
  tab <- data.frame(ind=rep(indNames(obj), each=pl))
  # column of pop ids to write
  if(pops){
    popnums <- 1:adegenet::nPop(obj)
    names(popnums) <- as.character(unique(adegenet::pop(obj)))
    popcol <- rep(popnums[as.character(adegenet::pop(obj))], each=pl)
    tab <- cbind(tab, data.frame(pop=popcol))
  }
  loci <- adegenet::locNames(obj) 
  # add columns for genotypes
  tab <- cbind(tab, matrix(-9, nrow=dim(tab)[1], ncol=adegenet::nLoc(obj),
                           dimnames=list(NULL,loci)))
  
  # begin going through loci
  for(L in loci){
    thesegen <- obj@tab[,grep(paste("^", L, "\\.", sep=""), 
                              dimnames(obj@tab)[[2]]), 
                        drop = FALSE] # genotypes by locus
    al <- 1:dim(thesegen)[2] # numbered alleles
    for(s in 1:S){
      if(all(!is.na(thesegen[s,]))){
        tabrows <- (1:dim(tab)[1])[tab[[1]] == indNames(obj)[s]] # index of rows in output to write to
        tabrows <- tabrows[1:sum(thesegen[s,])] # subset if this is lower ploidy than max ploidy
        tab[tabrows,L] <- rep(al, times = thesegen[s,])
      }
    }
  }
  
  # export table
  write.table(tab, file=file, sep="\t", quote=FALSE, row.names=FALSE)
  
  read_lines(file) %>%
    str_remove('ind') %>%
    write_lines(file)
}

runStructure2 <- function (path_to_structure, input_file, main_params, extra_params, 
                           out_prefix, n_K, n_replicates, n_cores) 
{
  #modified from https://github.com/sa-lee/starmie
  stopifnot(file.exists(path_to_structure))
  stopifnot(file.exists(input_file))
  stopifnot(file.exists(main_params))
  stopifnot(file.exists(extra_params))
  if (!is.integer(n_K) & n_K < 1) {
    stop("Assumed populations must be greater than 1")
  }
  if (!is.integer(n_replicates) & n_replicates < 1) {
    stop("Number of replicates must be greater than 1")
  }
  K <- 1L:n_K
  replicates <- 1L:n_replicates
  out_files <- outer(replicates, K, function(x, y) paste0(out_prefix, 
                                                          stringr::str_pad(x, width = 2, pad = 0), "_K_", stringr::str_pad(y, 
                                                                                                                           width = 2, pad = 0), ".out"))
  log_files <- gsub("out", "log", out_files)
  run_structure_single <- function(out_file, log_file) {
    k <- as.integer(stringr::str_extract(out_file, "[0-9]{2}\\b"))
    seed <- round(runif(1) * 1e+08)
    cmd <- paste(path_to_structure, "-K", k, "-i", input_file, 
                 "-m", main_params, "-e", extra_params, "-D", seed, 
                 "-o", out_file, ">", log_file)
    system(cmd)
    return(cmd)
  }
  message(paste("Running STRUCTURE on", n_cores, "core with", 
                n_K, "populations with", n_replicates, "replicates."))
  if (requireNamespace("parallel", quietly = TRUE)) {
    if (n_cores > parallel::detectCores()) {
      stop("Number of cores greater than number available on machine.")
    }
    cmds_run <- parallel::mcmapply(run_structure_single, 
                                   out_files, log_files, mc.cores = n_cores, mc.set.seed = TRUE)
    message(paste("Commands run\n", paste(cmds_run, collapse = "\n")))
  }
  else {
    if (n_cores > 1L) {
      stop("parallel package needed to run multiple calls. Please install it.", 
           call. = FALSE)
    }
    else {
      cmds_run <- mcmapply(run_structure_single, out_files, 
                           log_files)
      message(paste("Commands run\n", paste(cmds_run, 
                                            collapse = "\n")))
    }
  }
}

loadStructure2 <- function (filename, logfile = NULL) 
{
  #modified from https://github.com/sa-lee/starmie
  if (!is.na(filename) & !is.character(filename)) 
    stop("filename must be a string.")
  if (!is.null(logfile)) {
    if (!is.character(logfile) & !is.na(logfile)) 
      stop("logfile must be a string.")
  }
  structure_obj <- struct()
  s_f <- readr::read_lines(filename)
  s_f <- s_f[s_f != ""]
  prop_membership <- grep("^Proportion of membership.*", s_f)
  locprior <- TRUE
  if (length(prop_membership) == 0) {
    prop_membership <- grep("^Overall proportion of membership*", 
                            s_f)
    locprior <- FALSE
  }
  run_lines <- s_f[(grep("Run parameters:", s_f) + 1):(prop_membership - 
                                                         2)]
  run_lines <- str_trim(run_lines)
  run_lines <- str_split_fixed(run_lines, " ", n = 2)
  run_params <- data.frame(Parameter = run_lines[, 2], Value = run_lines[, 
                                                                         1], stringsAsFactors = FALSE)
  pops <- as.numeric(run_params[run_params$Parameter == "populations assumed", 
                                2])
  mem_lines <- s_f[(prop_membership + 3):(grep("^Allele-freq", 
                                               s_f) - 2)]
  mem_lines <- str_trim(mem_lines)
  if (!locprior) {
    mem_lines <- str_split_fixed(mem_lines, "\\s+", n = pops)
    mem_lines <- t(mem_lines)
    class(mem_lines) <- "numeric"
    mem_df <- data.frame(Cluster = mem_lines[, 1], Proportion = mem_lines[, 
                                                                          2], stringsAsFactors = FALSE)
  }
  else {
    mem_lines <- str_split_fixed(mem_lines, "\\s+", n = pops + 
                                   2)
    mem_lines[, 1] <- str_replace(mem_lines[, 1], ":", "")
    mem_df <- data.matrix(data.frame(mem_lines[-1, ], stringsAsFactors = FALSE))
    colnames(mem_df) <- mem_lines[1, ]
  }
  alle_lines <- s_f[(which(grepl("^Allele-freq", s_f)) + 4):which(grepl("^Average distances.*", 
                                                                        s_f)) - 1]
  alle_lines <- str_trim(alle_lines)
  alle_lines <- str_split_fixed(alle_lines, "\\s+", n = pops + 
                                  1)
  alle_freqs <- alle_lines[, 2:ncol(alle_lines)]
  suppressWarnings(class(alle_freqs) <- "numeric")
  avg_dist_lines <- s_f[(which(grepl("^Average distances.*", 
                                     s_f)) + 1):(which(grepl("^Estimated Ln.*", s_f)) - 2)]
  avg_dist_lines <- str_trim(avg_dist_lines)
  avg_dist_lines <- str_split_fixed(avg_dist_lines, "  : ", 
                                    n = 2)
  avg_dist_lines[, 1] <- str_replace(avg_dist_lines[, 1], 
                                     "cluster +", "")
  class(avg_dist_lines) <- "numeric"
  avg_dist_df <- data.frame(Cluster = avg_dist_lines[, 1], 
                            Avg.dist = avg_dist_lines[, 2])
  fit_lines <- s_f[(which(grepl("^Estimated Ln.*", s_f))):(which(grepl("^Mean value of Fst_1 .*", 
                                                                       s_f)) - 1)]
  fit_lines <- str_trim(fit_lines)
  fit_lines <- str_split_fixed(fit_lines, "= ", n = 2)
  fit_lines[, 1] <- str_trim(fit_lines[, 1])
  fit_stats_df <- data.frame(Statistic = fit_lines[, 1], Value = as.numeric(fit_lines[, 
                                                                                      2]))
  fst_lines <- s_f[(grep("^Mean value of Fst_1 .*", s_f)):(grep("^Inferred ancestry of.*", 
                                                                s_f) - 1)]
  fst_lines <- str_trim(fst_lines)
  fst_lines <- str_split_fixed(fst_lines, "= ", n = 2)
  fst_lines[, 1] <- str_trim(fst_lines[, 1])
  fst_lines[, 1] <- str_replace(fst_lines[, 1], "Mean value of Fst_", 
                                "")
  fst_df <- data.frame(Fst.Group = as.numeric(fst_lines[, 
                                                        1]), Value = as.numeric(fst_lines[, 2]))
  ances_lines <- s_f[(grep("^Inferred ancestry of.*", s_f) + 
                        1):(grep("^Estimated Allele Frequencies .*", s_f) - 
                              1)]
  ances_lines <- str_trim(ances_lines)
  header <- gsub("\\(|\\)|:", "", ances_lines[1], " ")
  header <- str_split(header, " ")[[1]]
  ances_lines <- ances_lines[-1]
  missing_proportions <- as.numeric(gsub("[\\(\\)]", "", regmatches(ances_lines, 
                                                                    regexpr("\\(.*?\\)", ances_lines))))
  sample_label <- str_trim(gsub("\\(", "", regmatches(ances_lines, 
                                                      regexpr(".*\\(", ances_lines))))
  sample_label <- str_split_fixed(sample_label, "\\s+", n = 2)[, 
                                                               -1]
  if (!locprior) {
    sample_summary <- data.frame(sample_label, missing_proportions)
    split_n <- 2
  }
  else {
    population_assignment <- as.integer(str_trim(gsub("\\)|:", 
                                                      "", regmatches(ances_lines, regexpr("\\).*?:", ances_lines)))))
    sample_summary <- data.frame(sample_label, missing_proportions, 
                                 population_assignment)
    split_n <- 3
  }
  ancest_matrix <- gsub(":  ", "", regmatches(ances_lines, 
                                              regexpr(":(.*)", ances_lines)))
  ancest_matrix <- str_split_fixed(ancest_matrix, "\\s+", 
                                   n = pops)
  class(ancest_matrix) <- "numeric"
  ancest_df <- data.frame(sample_summary, ancest_matrix, stringsAsFactors = FALSE)
  colnames(ancest_df)[1:split_n] <- header[1:split_n]
  colnames(ancest_df)[(split_n + 1):ncol(ancest_df)] <- paste("Cluster", 
                                                              seq(1, ncol(ancest_matrix)))
  clust_allel_lines <- s_f[(grep("^First column gives.*", 
                                 s_f) + 1):(grep("^Values of parameters used.*", s_f) - 
                                              1)]
  pos <- grep("^Locus .*", clust_allel_lines)
  clust_allel_lines <- gsub("[()%]", "", clust_allel_lines)
  clust_allel_lines <- unname(split(clust_allel_lines, cumsum(seq_along(clust_allel_lines) %in% 
                                                                pos)))
  clust_allele_list <- purrr::map(clust_allel_lines[1:2], 
                                  function(x) {
                                    list(Locus = as.numeric(str_split(x[[1]], "\\s+")[[1]][2]), 
                                         AlleleNumber = as.numeric(str_split(x[[2]], 
                                                                             "\\s+")[[1]][1]), MissingDataPercentage = as.numeric(str_split(x[[3]], 
                                                                                                                                            "\\s+")[[1]][1]), FreqMatrix = apply(str_split_fixed(str_trim(x[4:length(x)]), 
                                                                                                                                                                                                 "\\s+", n = pops + 2), 2, as.numeric))
                                  })
  structure_obj$K = pops
  structure_obj$run_params = run_params
  structure_obj$mem_df = mem_df
  structure_obj$allele_freqs = alle_freqs
  structure_obj$avg_dist_df = avg_dist_df
  structure_obj$fit_stats_df = fit_stats_df
  structure_obj$fst_df = fst_df
  structure_obj$ancest_df = ancest_df
  structure_obj$clust_allele_list = clust_allele_list
  if (!is.null(logfile)) {
    
    l_f <- readr::read_lines(logfile)
    l_f <- l_f[l_f != ""]
    
    column_names <- l_f %>%
      str_subset('Rep#') %>%
      str_trim %>%
      unique %>%
      str_subset('Est') %>%
      str_split('  +') %>%
      unlist
    
    full_log <- l_f %>%
      str_subset('^[0-9]+:') %>%
      str_subset('COPE', negate = TRUE) %>%
      str_trim %>%
      tibble(raw_lines = .) %>%
      separate(raw_lines, into = c(column_names), sep = '[: ]+', fill = 'right') %>%
      clean_names %>%
      mutate(ln_like = if_else(ln_like == '--', NA_character_, ln_like)) %>%
      mutate_all(as.numeric)  
    
    
    structure_obj$burn_df = full_log %>%
      filter(is.na(ln_like)) %>%
      as_draws_df() %>%
      mutate(.iteration = rep_number) %>%
      select(-rep_number) 
    
    structure_obj$nonburn_df = full_log %>%
      filter(!is.na(ln_like)) %>%
      as_draws_df() %>%
      mutate(.iteration = rep_number) %>%
      select(-rep_number)
  }
  structure_obj
}


microsats_genind <- full_cope_data %>%
  select(ID, starts_with('COPE'), starts_with('CPER')) %>%
  rowwise %>%
  filter(sum(!is.na(c_across(-ID))) != 0) %>%
  ungroup %>%
  select(-COPE10, -CPER52) %$%
  df2genind(.[,c(-1)], sep='/', ind.names=ID, NA.char = NA, type='codom') %T>%
  genind2structure('Structure_Run/microsats.structure')

# Takes a very long time to run
# runStructure2(path_to_structure = '~/Programs/structure', 
#               input_file = 'Structure_Run/microsats.structure', 
#               main_params = 'Structure_Run/mainparams',
#               extra_params = 'Structure_Run/extraparams',
#               out_prefix = 'Structure_Run/Results/',
#               n_K = 15, 
#               n_replicates = 10, 
#               n_core = 40)

structure_results <- tibble(output_file = list.files('Structure_Run/Results/', pattern = ".out_f$", full.names = TRUE)) %>%
  mutate(log_file = str_replace(output_file, ".out_f$", '.log'),
         replicate = str_extract(output_file, '/[0-9]+'),
         replicate = str_remove(replicate, '/'),
         replicate = as.integer(replicate),
         
         k = str_extract(output_file, 'K_[0-9]+'),
         k = str_remove(k, 'K_'),
         k = as.integer(k)) %>%
  mutate(structure_results = map2(output_file, log_file, loadStructure2)) 

#Confirm chain convergence - trace plots & rhat
structure_diagnostics <- structure_results %>%
  mutate(chains = map(structure_results, ~.x$nonburn_df)) %>%
  rowwise %>%
  mutate(chains = list(mutate(chains, .chain = replicate))) %>%
  group_by(k) %>%
  summarise(chains = list(do.call(rbind, chains)), .groups = 'rowwise') %>%
  mutate(diagnostics = list(select(chains, ln_like, starts_with('.')) %>% summary(rhat, ess_bulk, ess_tail)))

structure_diagnostics %>%
  select(-chains) %>%
  unnest(diagnostics) %>%
  filter(variable == 'ln_like')

structure_diagnostics %>%
  select(-diagnostics) %>%
  unnest(chains) %>%
  ggplot(aes(x = .iteration, y = ln_like, colour = as.character(.chain))) +
  geom_path() +
  facet_wrap(~k, labeller = label_both) +
  labs(colour = 'Chain', 
       x = 'Iteration',
       y = 'ln(likelihood)') +
  theme_classic()

## Choose Best K
structList(structure_results$structure_results) %>%
  bestK(method = 'evanno')

structList(structure_results$structure_results) %>%
  bestK(method = 'structure')

## Structure Assignments 
structure_assignments <- structure_results %>%
  filter(k == 2, replicate == 1) %>%
  pull(structure_results) %>%
  extract2(1) %>%
  plotBar(facet = FALSE, plot = FALSE) %>%
  as_tibble 

struct_plotting <- structure_assignments %>%
  group_by(Label) %>%
  mutate(plurality_clust = Cluster[value == max(value)]) %>%
  ungroup %>%
  
  spread(Cluster, value) %>%
  arrange(plurality_clust, desc(`Cluster 1`)) %>%
  mutate(Label = as.character(Label),
         Label = factor(Label, levels = .$Label)) %>%
  gather(Cluster, value, -Label, -plurality_clust)


the_box <- struct_plotting %>%
  mutate(Label_number = as.integer(Label)) %>%
  filter(value > 0.1 & value < 0.9) %>%
  filter(Label_number == min(Label_number) | Label_number == max(Label_number)) %>%
  select(Label) %>%
  distinct %>% 
  mutate(val = c("A", "B")) %>%
  pivot_wider(names_from = 'val',
              values_from = 'Label') %>%
  mutate(C = 0.1, D = 0.9)

struct_plotting %>%
  ggplot(aes(x = Label, y = value, fill = Cluster)) +
  geom_bar(stat = "identity", width=1) +
  geom_rect(data = the_box, aes(xmin = A, xmax = B, ymin = C, ymax = D), inherit.aes = FALSE, fill = NA, colour = 'white', size = 0.5) +
  scale_y_continuous(expand = c(0,0)) +
  # scale_fill_brewer(palette = 'Greys') +
  scale_fill_manual(values = c("#9f9f9fff", "#333333")) +
  ylab('Assignment Probability') +
  xlab("Individual") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10, colour = 'black'))
# ggsave("../Figures/Structure Plot.svg", width = 7, height = 7)

structure_classification <- structure_assignments %>%
  dplyr::rename(ID = Label, structure_cluster = Cluster, structure_assignment_prob = value) %>%
  group_by(ID) %>%
  filter(structure_assignment_prob == max(structure_assignment_prob)) %>%
  ungroup %>%
  mutate_if(is.factor, as.character)

full_cope_data %>%
  select(ID, morph_species, co1_species) %>%
  left_join(structure_classification, by = 'ID') %>%
  group_by(co1_species, structure_cluster) %>%
  summarise(n = n(), .groups = 'drop') %>%
  filter(co1_species != 'unknown', !is.na(structure_cluster))
#Cluster 1 is chya, cluster 2 is cpers - may flip based on RNG

full_cope_data <- full_cope_data %>%
  left_join(structure_classification, by = 'ID') %>%
  mutate(microsat_species = case_when(structure_cluster == 'Cluster 1' ~ 'chya',
                                      structure_cluster == 'Cluster 2' ~ 'cpers',
                                      TRUE ~ 'unknown')) 

#### Joint species ID ####
full_cope_data <- full_cope_data %>%
  mutate(joint_species = case_when(co1_species == microsat_species & microsat_species == morph_species ~ microsat_species,
                                   co1_species == microsat_species & morph_species == 'unknown' ~ microsat_species))

#### Calculate Shoal composition ####
## Flat prior
a <- 1; b<-1 
shoal_composition <- full_cope_data %>%
  group_by(Site, Shoal, X, Y) %>%
  summarise(cpers = sum(joint_species == 'cpers', na.rm = TRUE), 
            chya = sum(joint_species == 'chya', na.rm = TRUE), 
            total = cpers + chya,
            .groups = 'drop') %>%
  filter(total > 0) %>%
  ungroup %>%
  mutate(Shoal = LETTERS[as.integer(as.factor(Shoal))]) %>%
  
  mutate(prop_cpers = (cpers+a)/(total+a+b),
         cpers_lwr.95 = qbeta(0.025, cpers+a, total-cpers+b),
         cpers_lwr.50 = qbeta(0.25, cpers+a, total-cpers+b),
         cpers_upr.50 = qbeta(0.75, cpers+a, total-cpers+b),
         cpers_upr.95 = qbeta(0.975, cpers+a, total-cpers+b))

a <- 1; b<-1
percent_cpers_overall <- full_cope_data %>%
  group_by(Site) %>%
  summarise(n_shoal = n_distinct(Shoal), 
            cpers = sum(joint_species == 'cpers', na.rm = TRUE),
            chya = sum(joint_species == 'chya', na.rm = TRUE),
            total = cpers + chya) %>%
  mutate(site_cpers = (cpers+a)/(total+a+b),
         site_cpers_lwr.95 = qbeta(0.025, cpers+a, total-cpers+b),
         site_cpers_lwr.50 = qbeta(0.25, cpers+a, total-cpers+b),
         site_cpers_upr.50 = qbeta(0.75, cpers+a, total-cpers+b),
         site_cpers_upr.95 = qbeta(0.975, cpers+a, total-cpers+b)) %>%
  
  mutate(xmin = 'A\n25', xmax = 'M\n35')


a_shoal <- 1; b_shoal<-1; a_site <- 1; b_site<-1 # 
difference_from_site <- shoal_composition %>%
  dplyr::select(Site:total) %>%
  left_join(percent_cpers_overall %>%
              dplyr::select(Site, cpers, total) %>%
              dplyr::rename(site_cpers = cpers, site_total = total)) %>%
  mutate(mean_diff_site = (a_shoal+cpers)/(a_shoal+b_shoal+total)-(a_site+site_cpers)/(a_site+b_site+site_total),
         sd_diff_site = sqrt(((a_shoal+cpers)*(b_shoal+total-cpers))/((a_shoal+b_shoal+total)^2*(a_shoal+b_shoal+total+1)) + ((a_site+site_cpers)*(b_site+site_total-site_cpers))/((a_site+b_site+site_total)^2*(a_site+b_site+site_total+1))),
         lwr_95_site_diff = qnorm(0.025,mean_diff_site, sd_diff_site),
         upr_95_site_diff = qnorm(0.975,mean_diff_site, sd_diff_site),
         different = case_when(lwr_95_site_diff < 0 & upr_95_site_diff < 0 ~ TRUE,
                               lwr_95_site_diff > 0 & upr_95_site_diff > 0 ~ TRUE,
                               TRUE ~ FALSE))
difference_from_site %>%
  filter(Shoal %in% c('D', 'H')) %>%
  select(Shoal, lwr_95_site_diff, upr_95_site_diff)

shoal_composition %>%
  filter(Shoal %in% c('D', 'H'))


shoal_composition %>%
  ungroup %>%
  left_join(difference_from_site %>% 
              ungroup %>%
              dplyr::select(Site, Shoal, different)) %>%
  left_join(percent_cpers_overall %>%
              dplyr::select(-n_shoal:-total)) %>%
  mutate(majority = case_when(prop_cpers > site_cpers & different == TRUE ~ 'Cope',
                              prop_cpers < site_cpers & different == TRUE ~ 'Cohy',
                              TRUE ~ 'mix')) %>% 
  mutate(Shoal = str_c(Shoal, total, sep = '\n')) %>%
  
  
  ggplot(aes(x = Shoal, y = prop_cpers)) + #colour = majority
  
  geom_rect(data = percent_cpers_overall, aes(x = xmin, xmin = xmin, xmax = xmax, y = site_cpers, ymin = site_cpers_lwr.95, ymax = site_cpers_upr.95), colour = NA, fill = 'black', alpha = 0.3) +
  geom_rect(data = percent_cpers_overall, aes(x = xmin, xmin = xmin, xmax = xmax, y = site_cpers, ymin = site_cpers_lwr.50, ymax = site_cpers_upr.50), colour = NA, fill = 'black', alpha = 0.3) +
  geom_segment(data = percent_cpers_overall, aes(x = xmin, xend = xmax, y = site_cpers, yend = site_cpers), colour = 'black', size = 2, lty  = 4) +
  
  geom_linerange(aes(ymin = cpers_lwr.95, ymax = cpers_upr.95), lty = 1, size = 1) +
  geom_linerange(aes(ymin = cpers_lwr.50, ymax = cpers_upr.50), size = 3) +
  #geom_point(size = 6) +
  scale_y_continuous(limits = c(0,1), labels = scales::percent) +
  geom_point(data = tibble(prop_cpers = 0.8, Shoal = c('D\n25', 'H\n9')), shape = '*', size = 12) +
  # scale_color_manual('Shoal Composition', values = c('Cope'='grey50', 'Cohy'='grey75', 'mix'='black'), 
  #                    labels = c('C. hyalinus', 'C. personatus', 'Mixture')) +
  theme_classic() +
  ylab(expression(paste("Percent ", italic("C. personatus")))) +
  guides(color = guide_legend(title.position="top", title.hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 24),
        axis.text.y = element_text(colour = 'black', size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = 'black', size = 16),
        axis.ticks.x = element_blank(), 
        panel.background = element_rect(fill = "transparent", colour = NA), 
        panel.border = element_rect(fill = NA, colour = 'black'), 
        plot.background = element_blank(),
        legend.position='bottom',
        legend.title = element_text(colour = 'black', size = 20),
        legend.text = element_text(colour = 'black', size = 18),
        legend.background = element_blank(), 
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        
        strip.placement = "outside",
        strip.text = element_text(colour = 'black', size = 18))
# ggsave("../Figures/Shoal Composition_grey.svg", bg = "transparent", width = 7.5, height = 7.5)

#### Summary info for paper ####
full_cope_data %>%
  dplyr::select(morph_species, co1_species, microsat_species, joint_species) %>%
  gather(method, species) %>%
  group_by(method, species) %>%
  summarise(n = n()) %>%
  filter(!is.na(species),
         species != 'unknown') %>%
  mutate(prop = n/sum(n)) %T>%
  print %>%
  group_by(species) %>%
  summarise(mean(prop))


full_cope_data %>%
  dplyr::select(ID, morph_species, co1_species, microsat_species, structure_assignment_prob) %>%
  filter((morph_species != co1_species & morph_species != 'unknown' & co1_species != 'unknown') | 
           (microsat_species != co1_species & microsat_species != 'unknown' & co1_species != 'unknown') |
           (microsat_species != morph_species & microsat_species != 'unknown' & morph_species != 'unknown')) %>%
  mutate_at(vars(morph_species, co1_species, microsat_species), ~str_replace(., 'chya','C. hyalinus')) %>%
  mutate_at(vars(morph_species, co1_species, microsat_species), ~str_replace(., 'cpers','C. personatus')) %>%
  mutate_at(vars(morph_species, co1_species, microsat_species), ~str_replace(., 'unknown','-')) 



#### Venn diagram of sampling type overlap ####
df.venn <- data.frame(x = c(0, 0.866, -0.866),
                      y = c(1, -0.5, -0.5),
                      labels = factor(c('Morphological', 'Mitochondrial', 'Microsatellite'), 
                                      levels = c('Morphological', 'Mitochondrial', 'Microsatellite')))

venn.data <- full_cope_data %>%
  dplyr::select(ID, morph_species, co1_species, microsat_species) %>%
  dplyr::rename(morphology = morph_species, 
                co1 = co1_species, 
                microsat = microsat_species) %>%
  mutate_at(vars(morphology, co1, microsat), ~if_else(. == 'unknown', NA_character_, .)) %>%
  mutate_at(vars(morphology, co1, microsat), list(~!is.na(.))) %>%
  dplyr::select(-ID) %>%
  vennCounts() %>% 
  as.matrix %>%
  extract(-1,) %>% #the comma is correct
  as_tibble() %>%
  mutate(x = c(-1.2, 1.2, 0, 0, -0.8, 0.8, 0),
         y = c(-0.6, -0.6, -1, 1.2, 0.5, 0.5, 0))

label_data <- venn.data %>%
  filter(Counts %in% c(27, 8, 50)) %>%
  mutate(category = c('Microsatellite', 'Mitochondrial', 'Morphological'),
         y = 2 * c(-1.5, -1.5, 1.5),
         x = if_else(category == 'Morphological', x, 1.1 * x))


label_data <- df.venn %>%
  mutate(y = if_else(y < 0, y - 1.75, y + 1.75),
         x = if_else(y < 0, 1.5 * x, x))

ggplot(df.venn) +
  geom_circle(aes(x0 = x, y0 = y, r = 1.5), alpha = 0.5, size = 1, colour = 'black') +
  coord_fixed() +
  
  labs(fill = "Method of Identification") +
  annotate("text", x = venn.data$x, y = venn.data$y, label = venn.data$Counts, size = 18) +
  annotate("text", x = label_data$x, y = label_data$y, label = label_data$labels, size = 14) +
  theme_void() +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5)) +
  theme(legend.position = 'bottom',
        legend.title = element_text(colour = 'black', size = 20),
        legend.text = element_text(colour = 'black', size = 18),
        panel.background = element_rect(fill = "transparent", colour = NA), 
        panel.border = element_rect(fill = NA, colour = NA), 
        plot.background = element_blank())


#### Further Microsatellite Analysis ####
library(poppr)
library(pegas)
library(hierfstat)

just_certain_id <- full_cope_data %>%
  select(ID, starts_with('COPE'), starts_with('CPER'), morph_species, co1_species, microsat_species, 
         structure_assignment_prob, joint_species) %>%
  filter(co1_species != 'unknown',
         # structure_assignment_prob > 0.9,
         !is.na(joint_species)) %>%
  dplyr::select(ID, joint_species, COPE5:CPER188) %$%
  df2genind(.[,c(-1:-2)], sep='/', ind.names=ID, pop=joint_species, NA.char =NA, type='codom') 

#### Basic Summary
seppop(just_certain_id) %>%
  map(summary)

summary_stats <- basic.stats(just_certain_id, diploid = TRUE)

just_certain_id %>%
  seppop() %>%
  map(hw.test, B = 10000) %>%
  map(~as_tibble(.) %>%
        mutate(p_adj = p.adjust(Pr.exact, 'holm'))) %>%
  bind_rows(.id = 'species') %>%
  mutate(sig = p_adj < 0.05) %T>%
  print %>%
  filter(sig)


private_alleles(just_certain_id, count.alleles = FALSE, report = 'data.frame') %>%
  as_tibble() %>%
  separate(allele, into = c('locus', 'allele')) %>%
  group_by(population, locus) %>%
  summarise(count = sum(count), .groups = 'drop') %>%
  pivot_wider(names_from = 'population',
              values_from = 'count')



#### Inbreeding
summary_stats$Fis

tukey.nonadditivity.test <- function(the.aov) {
  the.model <- the.aov$model
  y <- the.model[,1]; f1 <- the.model[,2]; f2 <- the.model[,3]
  lm.1 <- lm(y~f1+f2)
  interact.term <- fitted(lm.1)^2
  lm.2 <- lm(y~f1+f2+interact.term)
  return(anova(lm.2)[3,])
}

library(broom)
library(emmeans)

fis_aov <- summary_stats$Fis %>%
  as_tibble(rownames = 'locus') %>%
  pivot_longer(cols = -locus,
               names_to = 'species',
               values_to = 'fis') %>%
  aov(fis ~ species + locus, data = .) 

tukey.nonadditivity.test(fis_aov)
diag.plots(fis_aov)
car::Anova(fis_aov, type = 2)
effectsize(car::Anova(fis_aov, type = 2), partial = TRUE)

emmeans(fis_aov, ~locus) %>% multcomp::cld(Letters = LETTERS) %>%
  as_tibble %>%
  mutate(.group = str_trim(.group),
         locus = fct_reorder(locus, emmean)) %>%
  ggplot(aes(x = locus, y = emmean, ymin = emmean - SE, ymax = emmean + SE)) +
  geom_pointrange() +
  geom_text(aes(y = 1.05 * (emmean + SE), label = .group))


#### Population Differentiation
fstat(just_certain_id, fstonly = TRUE)
gstat.randtest(just_certain_id, nsim = 9999)


just_certain_id %>%
  as.loci() %>%
  Fst

locus_specific_p <- just_certain_id %>%
  seploc() %>%
  map(gstat.randtest, nsim = 9999) %>%
  map_dbl(~.x$pvalue)

p.adjust(locus_specific_p, 'holm')