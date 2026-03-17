if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = c("tidyverse", "shape", "pafr"), character.only = TRUE)

get_species_name_from_ref_genome <- function(ref_name, position_in_name=2){
  #' Get the name of the species using the name of the reference genome
  #' @description
  #' This function splits the standardised name of the reference genome to get the name of the species.
  #' The standard name should look something like: `Jaera_albifrons_chromosomes.fasta`.
  #' @param ref_name This is a string with the name of the reference genome. It is like the example in the description.
  #' @param position_in_name `default=2`. This argument describes the position of the name to return in the name of the
  #' reference genome.
  #'
  #' @returns string. It returns the name of the species that was stored in the name of the reference genome.
  return(str_split_1(ref_name, "_")[position_in_name])
}

align_ref_genomes <- function(){
  #' Align reference genomes
  #' @description This way too long function is used to align reference genomes on one another.
  #' It iterates on all unique combination of reference genomes to algin them on one another. A part 
  #' of the function is taylored to my study as I am looking to align reference genomes in a specific
  #' order. The lines can be commented if it is not useful and these lines will be indicated in commentry.
  #' This function calls the minimap2 aligner using a bash script to run it in parallel for all the 
  #' combinations of reference genomes.
  #'
  #' @param None
  #'
  #' @returns None
  
  # Make the unique pairs of reference genomes
  combination_of_genomes <- combn(vector_ref_genomes, 2)
  # Iterate over these pairs of reference genomes
  for (col in 1:ncol(combination_of_genomes)){
    # Get the names of the reference genomes
    ref_genome_1 <- combination_of_genomes[1, col]
    ref_genome_2 <- combination_of_genomes[2, col]
    # Check that the firt genome is either praehirsuta or albifrons. Otherwise, change them
    # These lines can be commented as they are taylored to my study. You can comment from here ...
    if ((grepl("praehirsuta", ref_genome_2)) & (grepl("albifrons", ref_genome_1))){
      .temp_var. <- ref_genome_2
      ref_genome_2 <- ref_genome_1
      ref_genome_1 <- .temp_var.
      rm(.temp_var.)
    }
    # ... to here
    # Get the name of the output file
    name_output <- paste0(
      out_folder,
      get_species_name_from_ref_genome(ref_genome_2),
      "_aligned_on_",
      get_species_name_from_ref_genome(ref_genome_1),
      ".paf"
    )
    
    # Build the command line
    command_line  <- paste0("/shared/projects/sexisol/script/Basile/Cross_alignment_of_ref_genomes/Scripts/run_minimap2_for_reference_genomes.sh ",
                            folder_ref_genomes, ref_genome_1, " ", folder_ref_genomes, ref_genome_2, " ",
                            name_output)
    
    print(paste("command_line: sbatch", command_line))
    
    # Run minimap2 using the bash script on a cluster
    system2(command = "sbatch",
            args = paste0("/shared/projects/sexisol/script/Basile/Cross_alignment_of_ref_genomes/Scripts/run_minimap2_for_reference_genomes.sh ",
                          folder_ref_genomes, ref_genome_1, " ", folder_ref_genomes, ref_genome_2, " ",
                          name_output))
  }
}

are_there_paf_files_in_outfolder <- function(folder_to_check){
  #' Check if there are any paf files in a folder
  #' @description This function simply checks if there are any `.paf` files in the
  #' specified folder.
  #'
  #' @param folder_to_check This string is the path of the folder to check for `.paf files`
  #'
  #' @returns TRUE/FALSE The function returns TRUE if there is at least one `.paf` file in the
  #' specified folder
  
  number_paf_files_output <- length(list.files(folder_to_check, patt = ".paf"))
  return(number_paf_files_output > 0)
}

import_chromosome_length_from_dict <- function(path){
  #' Imports the lengths of the chromosomes using the ".dict" file from the
  #' reference genome indexing.
  #' @description This function takes as input the path to the ".dict" file
  #' obtained when the reference genome is indexed and keeps only the name
  #' of the chromosome and its length.
  #' @param path This string of characters is the path to the ".dict" file
  #' 
  #' @returns A table containing the name and length of each chromosomes.
  read.table(path, skip = 1) %>% 
    select(V2, V3) %>% 
    rename(Chromosome = V2,
           Length = V3) %>% 
    mutate(across(everything(), ~ str_split_fixed(., ":", 2)[, 2]),
           Length = as.numeric(Length),
           Chromosome = str_split_fixed(Chromosome, "\\.", 2)[, 1]) %>% 
    return()
}

import_paf <- function(path){
  #' Import a "paf" file and transform it into a tidy table
  #' @description This function uses the `read_paf` function from the "pafr"
  #' package to import a ".paf" file and transforms the list into a tidy table
  #' that can me modified or plotted with the tidy format.
  #' @param path This string of characters is the path to the "paf" file
  #' 
  #' @returns A table with the contents of the "paf" file
  read_paf(path) %>% 
    plyr::ldply() %>% 
    column_to_rownames(".id") %>% 
    t %>% 
    as_tibble() %>% 
    mutate(across(everything(), ~ suppressWarnings(convert_to_numeric_if_possible(.)))) %>% 
    return()
}

get_corresponding_chromosomes_species <- function(all_alignments){
  #' Get the correspondences between chromosomes of different species
  #' @description This function uses the alignments of the reference genomes got
  #' from the "paf" files to know which chromosomes are the closest, thus 
  #' which chromosome will be called "Chrom1" in the different species.
  #' To do this, we get the total length of aligned sequences per pair of 
  #' chromosomes and filter it to keep only the pair of chromosomes for which
  #' the total length of aligned sequence is the greatest.
  #' @param all_alignments This data frame contains all the alignments that we
  #' want to plot.
  #' 
  #' @returns A table with four columns: two columns with the names of the species
  #' the name of the chromosome in the first species and the name of the 
  #' corresponding chromosome in the second species
  all_alignments %>% 
    group_by(query, target, qname, tname) %>% 
    summarise(total_length = sum(alen), .groups = "drop_last") %>% 
    arrange(desc(total_length)) %>% 
    mutate(sp_pairs = paste(query, target, sep = "_")) %>% 
    dplyr::slice_head(n = 2) %>% 
    mutate(name = row_number()) %>% 
    pivot_wider(names_from = "name", values_from = c("tname", "total_length")) %>% 
    find_correspondences_for_rearrangements() %>%  
    select(-total_length) %>%
    ungroup %>% 
    return()
}

order_magnitude <- function(x){
  #' Get the order of magnitude of a number
  #' @description
    #' This function gets the order of magnitude of a number using the `floor`
    #' and `log10` functions.
  #' @param x This is a number.
  #' 
  #' @returns An integer
  
  return(floor(log10(abs(x))))
}

is_same_order_of_magnitude <- function(a, b){
  #' Checks if two numbers have the same order of magnitude
  #' @description
    #' This function computes the order of magnitude of the two provided numbers
    #' using the `order_magnitude` function provided above.
  #' @param a This is a number.
  #' @param b This is a number.
  #' 
  #' @returns Boolean: `TRUE` if the two numbers have the same order of
  #' magnitude and `FALSE` otherwise.
  
  order_magnitude_a <- order_magnitude(a)
  order_magnitude_b <- order_magnitude(b)
  
  return(order_magnitude_a == order_magnitude_b)
}

add_not_attributed_chromosomes <- function(information_rearrangements, chroms_not_attributed){
  #' Adds chromosomes that were not added in the `find_correspondences_for_rearrangements`
  #' function.
  #' @description
    #' This function selects the column of the input containing chromosomes that
    #' were not added to the list of correspondences between chromosomes and adds
    #' them to the list of correspondences with the chromosome that corresponds
    #' the most.
  #' @param information_rearrangements This data frame contains the top two 
  #' correspondences for each chromosome (two chromosomes that have the largest 
  #' cumulative matching sequence). 
  #' @param chroms_not_attributed This vector contains the names of the chromosomes
  #' that did not get attributed to another chromosome yet.
  #' 
  #' @returns A data frame containing the correspondence between the chromosomes 
  #' that were not attributed and their most corresponding chromosome in another
  #' species.
  
  # The name of the chromosome can be either in the query or target columns, so 
  # to select the right column, we need to know where it is located. To do this,
  where_are_unattributed_chroms <- information_rearrangements %>% 
    # we first try to match both columns to the chromosomes that were not attributed
    mutate(across(contains("name"), ~ . %in% chroms_not_attributed)) %>% 
    # we select only these columns
    select(contains("name")) %>% 
    # And we sum these columns
    colSums(na.rm = TRUE)
  # Once this is done, we can filter the name of the column containing the name 
  # of the chromosome that was not attributed
  colname <- names(where_are_unattributed_chroms)[where_are_unattributed_chroms >= 1]
  
  # Then, we need to know which columns to select depending on the column name
  if (colname %in% c("tname_1", "query")){
    columns_to_select <- c("query", "target", "qname", "tname_1", "total_length_1")
    length_column <- "total_length_1"
  }else{
    columns_to_select <- c("query", "target", "qname", "tname_2", "total_length_2")
    length_column <- "total_length_2"
  }
  
  # And finally, we filter and select the right columns to be returned.
  information_rearrangements %>% 
    filter(!!sym(colname) == chroms_not_attributed) %>% 
    select(all_of(columns_to_select)) %>% 
    rename(!!sym(str_split_fixed(colname, "_", 2)[, 1]) := !!sym(colname),
           total_length := !!sym(length_column)) %>% 
    return()
  
}

find_correspondences_for_rearrangements <- function(information_rearrangements){
  #' Get the names of the corresponding chromosomes when there are some complex 
  #' rearrangements and the name of the chromosome that matches most is not enough
  #' @description
    #' This function takes the top two chromosomes of species B that match the 
    #' chromosome in species A and counts the number of times these top two
    #' chromosomes of species B are the same. If the pair of chromosomes of species
    #' B is unique, the first chromosome of species B (the one with the largest 
    #' cumulative mapping sequence) is given as correspondence to the chromosome
    #' of species A. If the pair of chromosomes of species B is present more than
    #' once, it checks the order of magnitude of the cumulative mapping sequence
    #' for the top two chromosomes of species B. If the order of magnitude is 
    #' not the same, the first chromosome of species B (the one with the largest 
    #' cumulative mapping sequence) is given as correspondence. If the order of
    #' magnitude is the same, the function attributes the first chromosome of
    #' species A to the first chromosome of species B (the one with the largest 
    #' cumulative mapping sequence) and attributes the remaining chromosome of 
    #' species A to the remaining chromosome of species B.
  #' @param information_rearrangements This data frame contains the top two
  #' chromosomes of species B that matches each chromosome of species A.
  #' 
  #' @returns A data frame that contains all the names of chromosomes from
  #' species A and the name of the corresponding chromosome of species B.

  # Get the names of all the chromosomes in the target species
  tnames <- information_rearrangements %>%
    ungroup %>% 
    select(tname_1, tname_2) %>%
    pivot_longer(cols = everything(), names_to = "toto", values_to = "tname") %>%
    pull(tname) %>%
    unique
  
  # Get the names of all chromosomes in the query species
  qnames <- information_rearrangements %>% 
    pull(qname) %>% 
    unique
  
  # For each pair of chromosome from species B, count the number of pairs
  counts_correspondences_chromosomes <- information_rearrangements %>% 
    group_by(sp_pairs, tname_1, tname_2) %>% 
    summarize(count = n(),
              .groups = "drop_last") %>% 
    left_join(information_rearrangements,
              by = c("sp_pairs", "tname_1", "tname_2"))
  
  # For pairs of chromosomes of species B that only exist once, the first
  # chromosome of the pair is used as correspondence for the chromosome of
  # species A.
  obvious_correspondences <- counts_correspondences_chromosomes %>% 
    ungroup %>% 
    filter(count == 1) %>% 
    select(query, target, qname, tname_1, total_length_1) %>% 
    rename(tname = tname_1,
           total_length = total_length_1)
  
  # For the pairs of chromosomes from species B that are present more than once,
  # we compare the order of magnitude of the cumulative matching sequence on the
  # two chromosomes. We differentiate fusions from translocations.
  complex_rearrangements <- counts_correspondences_chromosomes %>% 
    ungroup() %>% 
    filter(count > 1) %>% 
    rowwise %>% 
    mutate(Rearrangement_type = ifelse(
      is_same_order_of_magnitude(total_length_1, total_length_2),
      "Large_translocation", "Fusion"
    )) 
  
  # For the fusions, the first chromosome of species B is taken as correspondence
  complex_fusions <- complex_rearrangements %>% 
    filter(Rearrangement_type == "Fusion") %>% 
    select(query, target, qname, tname_1, total_length_1) %>% 
    rename(tname = tname_1,
           total_length = total_length_1)
  
  # For the translocations, the first chromosome of species A is attributed to 
  # the first chromosome of species B and the second chromosome of species A is
  # attributed to the second chromosome of species B.
  complex_translocations <- complex_rearrangements %>% 
    filter(Rearrangement_type == "Large_translocation") %>%
    group_by(tname_1, tname_2) %>%
    mutate(tname = case_when(
      row_number() == 1 ~ tname_1,
      row_number() == 2 ~ tname_2
    ),
    total_length = case_when(
      row_number() == 1 ~ total_length_1,
      row_number() == 2 ~ total_length_2
    )) %>% 
    ungroup %>% 
    select(query, target, qname, tname, total_length)
  
  # Here, we combine the tables that were created
  all_rearrangements <- obvious_correspondences %>% 
    rbind(complex_fusions,
          complex_translocations)
  
  # In these manipulations, some chromosomes may not have been attributed, so
  # we add them here
  chroms_not_attributed <- which_not_in(c(qnames, tnames), c(all_rearrangements$qname, all_rearrangements$tname)) %>% 
    discard(is.na)
  
  missing_chromosomes <- information_rearrangements %>% 
    ungroup %>% 
    add_not_attributed_chromosomes(chroms_not_attributed)
  
  all_rearrangements %>% 
    rbind(missing_chromosomes) %>% 
    return()
}

get_chromosome_names_using_reference <- function(correspondences_chromosomes, reference_chromosomes, pair_of_species, reference_species){
  #' Get the name of the chromosomes using the names of one species as a reference
  #' @description This function allows to name the chromosomes of one species
  #' using the names given to the chromosomes of a reference species. We use
  #' correspondences between the chromosome of the reference species and the
  #' chromosome of the focal species to return the correct names.
  #' @param correspondences_chromosomes This data frame contains the
  #' correspondences between the chromosomes of the focal species and the
  #' reference species.
  #' @param reference_chromosomes This data frame contains three columns: the
  #' name of the reference species, the names of the chromosomes given when they
  #' were assembled and the name of the chromosomes that are given when doing the
  #' analysis.
  #' @param pair_of_species This string contains the name of the focal species
  #' and the name of the reference species separated by a "_".
  #' 
  #' @returns A table with three columns: the name of the focal species, the name
  #' of the chromosomes of this species when they were assembled and the new name
  #' given when doing the data analysis.
  
  # First, we separate the names of the species
  species <- str_split_fixed(pair_of_species, "_", 2) %>% 
    as.vector

  # And consider which columns to sample depending on if the reference species was
  # used to map or was mapped on in the minimap2 step.
  if (species[1] == reference_species){
    columns_to_return <- c("target", "tname")
  }else{
    columns_to_return <- c("query", "qname")
  }
  
  # Finally, we join the table with the reference species 
  correspondences_chromosomes %>% 
    filter(query == species[1],
           target == species[2]) %>% 
    left_join(reference_chromosomes %>% 
                select(-Species),
              by = join_by("qname" == "Chromosome")) %>% 
    left_join(reference_chromosomes %>% 
                select(-Species),
              by = join_by("tname" == "Chromosome"),
              relationship = "many-to-many") %>% 
    # Depending on where the join was done, there are some chromosome names in
    # two columns, so we fuse these two columns.
    mutate(Chromosome_name.x = case_when(
      !is.na(Chromosome_name.y) ~ Chromosome_name.y,
      TRUE ~ Chromosome_name.x
    )) %>% 
    select(-Chromosome_name.y) %>%
    rename(Chromosome_name = Chromosome_name.x) %>%
    select(all_of(columns_to_return), Chromosome_name) %>% 
    rename(Species = columns_to_return[1],
           Chromosome = columns_to_return[2]) %>% 
    return()
}

define_reference_for_chromosome_naming <- function(all_chromosomes, prefix = "Chrom"){
  #' Define a reference species and chromosome names for the chromosome naming
  #' @description This function takes the first species alphabetically in the
  #' table of all chromosomes from the different species, orders the chromosomes
  #' from the longest to the shortest and names them from "Chrom_1" to "Chrom_n",
  #' n being the number of chromosomes that are given for the species chosen as 
  #' reference.
  #' @param all_chromosomes This data frame contains at least three columns: 
  #' the name of the chromosomes that was given to them when they were 
  #' assembled, the length of the chromosomes and the species to which the 
  #' chromosomes belong.
  #' @param prefix This string is the beginning of the name of the chromosomes.
  #' Its default is "Chrom" to have "Chrom_1", ..., but it can be changed to 
  #' any string. Usual ones are "Chr", "Chrom" or "LG".
  #' 
  #' @returns This function returns a table containing three columns: The name 
  #' of the chromosome that was given when it was assembled for the species 
  #' chosen as reference, the new name of the chromosome and the species that was
  #' chosen as reference.
  
  # Define the species to use as reference
  reference_species <- all_chromosomes %>% 
    arrange(Species) %>% 
    dplyr::slice(1) %>% 
    pull(Species)
  
  all_chromosomes %>% 
    filter(Species == reference_species) %>% 
    arrange(desc(Length)) %>% 
    mutate(Chromosome_name = paste(prefix, 1:nrow(.), sep = "_")) %>% 
    select(Chromosome, Chromosome_name, Species) %>% 
    return()
}

are_all_in <- function(vec1, vec2){
  #' Check if all occurrences of the first vector are in the second one
  #' @description Checks if all the occurrences of the first 
  #' vector are also contained in the second one
    #' @param vec1 The vector to test
    #' @param vec2 The vector in which to test
  #' 
  #' @returns boolean TRUE if all the occurrences of the first vector are also in
  #' the second one and FALSE otherwise.
  
  return(all(vec1 %in% vec2))
}

get_new_reference_species <- function(sps, already_done_species){
  #' Get a new reference species in the ones that were already run
  #' @description
    #' This function checks if the first species given has already been run. If
    #' it is the case, it uses this as a new reference. Otherwise, it uses the 
    #' second species specified.
  #' @param sps Vector containing two species names to test
  #' @param already_done_species vector containing all the species for which the
  #' chromosome naming was already done
  #' 
  #' @returns the name of the new reference species to use
  if (sps[1] %in% already_done_species){
    new_ref_species <- sps[1]
  }else{
    new_ref_species <- sps[2]
  }
  return(new_ref_species)
}

any_string_in_vector_contains_pattern <- function(pattern, vec_to_check){
  #' Checks if any string in a vector contains a pattern
  #' @description
    #' Checks if any string in the vector contains the specified pattern
  #' @param pattern Pattern to look for in the vector
  #' @param vec_to_check Vector in which to look for the pattern
  #' 
  #' @returns Boolean: TRUE if the pattern is in any of the vector's string and
  #' FALSE otherwise
    return(any(grepl(pattern, vec_to_check)))
}

iterate_over_pairs_of_species_to_get_name_of_chromosomes <- function(pair_of_species_to_do,
                                                                     already_done_pairs_of_species,
                                                                     named_chromosomes_all_sps,
                                                                     correspondences_chromosomes,
                                                                     reference_chromosomes,
                                                                     reference_species = NULL){
  #' Go over all the pair of species given as parameter to name the chromosomes
  #' according to a reference species
  #' @description
    #' This function uses a for loop to iterate over all the combination of species
    #' that are given in input to name all the corresponding chromosomes between
    #' pairs of species. A reference can be specified, but is not necessarily
    #' needed.
    #' This function is called recursively, which is why some of its outputs are
    #' also given in input. If you run this function for the first time, the 
    #' arguments that are given in the output can be given as empty.
  #' @param pair_of_species_to_do Vector containing all the pairs of species to
  #' iterate over in the function. The pairs of species must be under the
  #' following format: "sp1_sp2"
  #' @param already_done_pairs_of_species Vector containing all the species for 
  #' which the chromosomes are already named. In this vector, each component is 
  #' the name of one species.
  #' @param named_chromosomes_all_sps Data frame containing the name that was 
  #' given to a chromosome when it was assembled, the name of the chromosome that
  #' is given during the analysis and the name of the species in which this 
  #' naming was done. If it is the first time that the function is run, this can
  #' be an empty data frame.
  #' @param correspondences_chromosomes Data frame containing four columns:
  #' the names of the query and target species and the names of their 
  #' corresponding chromosomes. This is the output of the function `get_corresponding_chromosomes_species`.
  #' @param reference_chromosomes Data frame containing three columns: the name 
  #' of the species used as reference, the name of the chromosomes for this reference
  #' species when it was assembled and the name of the chromosomes that were 
  #' given in the data analysis.
  #' @param reference_species Name of the reference species.
  #' 
  #' @returns list of three parameters:
  #'  1- named_chromosomes_all_sps: Data frame containing the name that was 
  #'    given to a chromosome when it was assembled, the name of the chromosome that
  #'    is given during the analysis and the name of the species in which this 
  #'    naming was done
  #'  2-already_done: vector containing the list of the species that were already done
  #'  3-to_do_later: vector containing the pairs of species to run in another call to the function

  # Initialise the loop
  to_do_later <- c()
  for (pair_of_species in pair_of_species_to_do){
    # Separate the pair to get independent names
    sps <- str_split_fixed(pair_of_species, "_", 2) %>% 
      as.vector
    # Check if the pair was already done
    if (are_all_in(sps, already_done_pairs_of_species)) next
    
    if (is.null(reference_species)){
      # If none of the species in the pair was done, we add it to the list to do later
      # as we can't compare any of the species to any other one.
      if (!any(sps %in% already_done_pairs_of_species)){
        to_do_later <- c(to_do_later, pair_of_species)
        next
      }
      # Otherwise, we use one of them as the new reference species
      reference_species <- get_new_reference_species(sps, already_done_pairs_of_species)
      reference_chromosomes <- named_chromosomes_all_sps %>% 
        filter(Species == reference_species)
    }else{
      # We add the reference chromosomes' names to the output table if it is not
      # already done
      if (reference_species %!in% unique(named_chromosomes_all_sps$Species)){
        named_chromosomes_all_sps <- add_table_to_df_in_iteration(named_chromosomes_all_sps,
                                                                  reference_chromosomes)
      }
      # Here, we check that the species used as reference is in the pair of species to 
      # run as we need it to name the chromosomes of the focal species.
      if (!grepl(reference_species, pair_of_species)){
        to_do_later <- c(to_do_later,
                         pair_of_species)
        next
      }
    }
    # We get the names of the chromosomes of the focal species.
    table_corresp_pair_species <- get_chromosome_names_using_reference(
      correspondences_chromosomes = correspondences_chromosomes,
      reference_chromosomes = reference_chromosomes,
      pair_of_species = pair_of_species,
      reference_species = reference_species
    )
    
    already_done_pairs_of_species <- c(already_done_pairs_of_species,
                                       sps) %>% unique
    named_chromosomes_all_sps <- add_table_to_df_in_iteration(named_chromosomes_all_sps,
                                                              table_corresp_pair_species)
  }
  return(list(
    "named_chromosomes_all_sps" = named_chromosomes_all_sps,
    "already_done" = already_done_pairs_of_species,
    "to_do_later" = to_do_later
  ))
}

contains_dot <- function(string){
  #' This function checks if there is a dot in the string
  #' @description
    #' This function uses the `grepl` function to check the string for a dot
  #' @param string String in which to look for a dot
  #' 
  #' @returns Boolean: `TRUE` if the pattern contains a dot and `FALSE` otherwise
  
  return(grepl("\\.", string))
}

contains_hyphen <- function(string){
  #' This function checks if there is a hyphen in the string
  #' @description
  #' This function uses the `grepl` function to check the string for a hyphen
  #' @param string String in which to look for a hyphen
  #' 
  #' @returns Boolean: `TRUE` if the pattern contains a hyphen and `FALSE` otherwise
  
  return(grepl("-", string))
}

get_names_of_fused_chromosomes <- function(string){
  #' Separates and distributes a prefix separated by a dot to two suffixes separated
  #' by hyphens
  #' @description
    #' This function splits the string using a hyphen as delimiter, splits the string
    #' using a dot as delimiter and returns the first part of the string (before
    #' the hyphen) and the second part of the string (after the hyphen) with the
    #' suffix (part after the dot).
  #' @param string The string to be split
  #' 
  #' @returns A list with two levels: `name_first_chrom` is the prefix (part
  #' before the dot) of the string, and `name_second_chrom` is the prefix with
  #' the second part of the string (after the hyphen).
  #' 
  #' @examples
    #' get_name_of_fused_chromosomes("LG1.a-b")
    #' ## This returns the following list
    #' # $name_first_chrom
    #' # [1]  "LG1.a"
    #' # $name_second_chrom
    #' # [2]  "LG1.b"
  
  portions <- str_split_fixed(string, "-", 2) %>% 
    as.vector
  
  prefix <- (str_split_fixed(portions[1], "\\.", 2) %>% 
    as.vector)[1]
  
  return(list(
    "name_first_chrom" = portions[1],
    "name_second_chrom" = paste(prefix, portions[2], sep = ".")
  ))
}

get_names_unfused_chromosomes <- function(string){
  #' Separates and distributes a prefix separated by a dot to a suffix and the
  #' following letter of the alphabet.
  #' @description
    #' This function splits the string using a dot as delimiter, identifies which
    #' letter constitutes the suffix and returns the string as well as the prefix
    #' pasted with the next letter of the alphabet.
  #' @param string The string to be split
  #' 
  #' @returns A list with two levels: `name_first_chrom` is the string given as
  #' input, and `name_second_chrom` is the prefix (part of the string before the 
  #' dot) with the next letter of the alphabet.
  #' 
  #' @examples
    #' get_name_of_fused_chromosomes("LG1.a")
    #' ## This returns the following list
    #' # $name_first_chrom
    #' # [1]  "LG1.a"
    #' # $name_second_chrom
    #' # [2]  "LG1.b"
  
  prefix <- str_split_fixed(string, "\\.", 2) %>% 
    as.vector
  letter_nb <- which(letters == prefix[2])
  
  return(list(
    "name_first_chrom" = string,
    "name_second_chrom" = paste(prefix[1], letters[letter_nb + 1], sep = ".")
  ))
}

rename_portions_of_chromosomes_fusion <- function(temp_names_for_sp_for_chr, all_alignments, chromosome_correspondences){
  #' Reorder and rename the portions of chromosomes that are in different orders
  #' in different species (especially the case for fusions/fissions of chromosomes).
  #' @description
    #' This function computes the mean mapping position of all the chromosomes of
    #' species B that map on one chromosome of species A in the case of fusions.
    #' It reorders the fragments of the fractionned chromosome for them to be in
    #' accordance with the chromosome of species A. For example, if we have LG1.a-b
    #' in species A that is separated into LG1.a and LG1.b in species B and that 
    #' the portion called LG1.a maps with a portion of LG1.a-b that is located
    #' after the portion of LG1.a-b that maps with LG1.b, this function will 
    #' swap the names of LG1.a and LG1.b so that the names make sense with the 
    #' mapping.
  #' @param temp_names_for_sp_for_chr This data frame contains the names of the 
  #' chromosomes of species A and B, their respective mapping positions and the
  #' names of species A and B.
  #' @param all_alignments This data frame contains all the alignments between the
  #' two species.
  #' @param chromosome_correspondences This data frame contains the correspondences
  #' between chromosomes of the two species. It is the output of the `get_corresponding_chromosomes_species`
  #' function.
  #' 
  #' @returns A data frame containing the new names of the chromosomes
  
  # First we compute the mean mapping position of all the chromosomes of species
  # B that map on one chromosome of species A 
  positions_to_where_chroms_map <- temp_names_for_sp_for_chr %>%
    left_join(correspondence_chroms, by = join_by("Chromosome" == "qname")) %>%
    left_join(all_correspondences, by = join_by("Chromosome" == "qname", "query", "target", "tname")) %>%
    group_by(Chromosome, tname, query, target) %>%
    summarize(mean_map = mean(tstart),
              .groups = "drop_last") %>% 
    ungroup %>% 
    rename(qname = Chromosome)

  # We get the names of the considered chromosomes (both the name given in the
  # analysis and the names given in the assembly)
  name_chrom_given <- temp_names_for_sp_for_chr %>% 
    pull(Chromosome_name) %>% 
    unique
  
  name_first_chrom_assembly <- positions_to_where_chroms_map %>% 
    filter(mean_map == min(mean_map)) %>% 
    pull(qname)
  
  name_second_chrom_assembly <- positions_to_where_chroms_map %>% 
    filter(mean_map == max(mean_map)) %>% 
    pull(qname)
  
  # Then, we consider how to rename the chromosomes given their original name
  if (name_first_chrom_assembly != name_second_chrom_assembly){
    # Here, we have a fusion of the reference chromosome into two chromosomes
    if (contains_dot(name_chrom_given)){
      # If we already know the fusion, we treat it as such
      if (contains_hyphen(name_chrom_given)){
        names_to_give_to_chromosomes <- get_names_of_fused_chromosomes(name_chrom_given)
      }else{
        names_to_give_to_chromosomes <- get_names_unfused_chromosomes(name_chrom_given)
      }
    }
  }
  # And give the new names to the chromosomes
  temp_names_for_sp_for_chr %>% 
    mutate(Chromosome_name = case_when(
      Chromosome == name_first_chrom_assembly ~ names_to_give_to_chromosomes$name_first_chrom,
      Chromosome == name_second_chrom_assembly ~ names_to_give_to_chromosomes$name_second_chrom
    )) %>% 
    return()
}

rename_chromosomes_after_calling <- function(temp_names_chroms, all_alignments, chromosome_correspondences){
  #'
  
  final_names_chromosomes <- data.frame()
  for (species in unique(temp_names_chroms$Species)){
    # print(species)
    temp_names_for_sp <- temp_names_chroms %>% 
      filter(Species == species)
    for (chromosome in unique(temp_names_for_sp$Chromosome_name)){
      # print(paste0("   ", chromosome))
      temp_names_for_sp_for_chr <- temp_names_for_sp %>% 
        filter(Chromosome_name == chromosome)
      
      if (nrow(temp_names_for_sp_for_chr) == 1){
        final_names_chromosomes <- add_table_to_df_in_iteration(final_names_chromosomes,
                                                                temp_names_for_sp_for_chr)
      }else{
        final_names_chromosomes <- final_names_chromosomes %>%
          add_table_to_df_in_iteration(temp_names_for_sp_for_chr %>%
                                         rename_portions_of_chromosomes_fusion(all_alignments, chromosome_correspondences))
      }
    }
    
    for (chromosome in unique(temp_names_for_sp$Chromosome)){
      # print(paste0("second    ", chromosome))
      temp_names_for_sp_for_chr <- temp_names_for_sp %>% 
        filter(Chromosome == chromosome)
      if (nrow(temp_names_for_sp_for_chr) == 1) next
      else{
        letters_in_rearrangement <- temp_names_for_sp_for_chr %>% 
          mutate(letter_in_chrom_name = str_split_fixed(Chromosome_name, "\\.", 2)[, 2]) %>% 
          arrange(letter_in_chrom_name) %>%
          pull(letter_in_chrom_name) %>% 
          paste(collapse = "-")
        
        
        final_names_chromosomes <- temp_names_for_sp_for_chr %>% 
          mutate(Chromosome_name = str_split_fixed(Chromosome_name, "\\.", 2)[, 1] %>% 
                   paste(letters_in_rearrangement, sep = ".")) %>% 
          unique %>% 
          rbind(final_names_chromosomes %>% 
                  filter(Chromosome != chromosome))
      }
    }
  }
  return(final_names_chromosomes)
  
}

get_chromosome_names_all_species <- function(all_chromosomes, correspondences_chromosomes, all_alignments, reference_chromosomes = NULL){
  #' Get the names of all the chromosomes for all the species that are specified
  #' @description
    #' This function uses a species as reference to name the chromosomes of all the species 
    #' in the same way. For example, if we have two species, A and B, and we
    #' know that chromosome 1 of species A corresponds to chromosome OZ1292456
    #' from species B, this function will allow to name the chromosome OZ1292456
    #' as Chromosome 1 in species B. This function automates this process for all
    #' the pairs of species that are aligned on each other.
  #' @param all_chromosomes This data frame contains at least the name of the
  #' chromosomes that were given when the chromosomes were assembled, and to which
  #' species this assembly corresponds.
  #' @param correspondences_chromosomes This data frame contains four columns: 
  #' the species of the query and target species and the names of the corresponding
  #' chromosomes from these two species. This table is the output of the `get_corresponding_chromosomes_species`
  #' function.
  #' @param reference_chromosomes This data fame contains the name of the species
  #' used as reference, the names of the chromosomes when assembled and the names
  #' of the chromosomes in the data analysis. The default value of this argument 
  #' is `NULL`. If this is the case, a default species will be chosen as reference.
  #' 
  #' @returns This function returns a table with at least three columns: the names
  #' of all species, the names of all chromosomes that were given when assembled
  #' and the names of the chromosomes given in the data analysis. Any supplementary
  #' column contained in the `all_chromosomes` table will also be found in the 
  #' output table.
  
  # Check if there is a reference species. If not, one will be appointed
  if (is.null(reference_chromosomes)){
    reference_chromosomes <- define_reference_for_chromosome_naming(all_chromosomes)
  }
  
  # Get the name of the reference species given as input or appointed
  reference_species <- reference_chromosomes %>% 
    pull(Species) %>%
    unique
  
  # Get all the pair of species that can be considered
  pairs_species <- correspondences_chromosomes %>% 
    mutate(pair_species = paste(query, target, sep = "_")) %>% 
    pull(pair_species) %>% 
    unique
  
  # Initialize variables in loop
  already_done_species <- c()
  named_chromosomes_all_sps <- data.frame()
  
  while (!is.null(pairs_species)){
    
    list_outputs_iteration <- iterate_over_pairs_of_species_to_get_name_of_chromosomes(
      pair_of_species_to_do = pairs_species,
      already_done_pairs_of_species = already_done_species,
      named_chromosomes_all_sps = named_chromosomes_all_sps,
      correspondences_chromosomes = correspondences_chromosomes,
      reference_chromosomes = reference_chromosomes,
      reference_species = reference_species
    )
    
    pairs_species <- list_outputs_iteration$to_do_later
    already_done_species <- list_outputs_iteration$already_done
    named_chromosomes_all_sps <- list_outputs_iteration$named_chromosomes_all_sps
    
    if ((!is.null(reference_species)) && (!any_string_in_vector_contains_pattern(reference_species, pairs_species))){
      reference_species <- NULL
      reference_chromosomes <- NULL
    }
  }
  
  named_chromosomes_all_sps %>% 
    rename_chromosomes_after_calling(all_alignments) %>% 
    left_join(all_chromosomes,
              by = c("Chromosome", "Species")) %>% 
    return()
}

reverse_chromosome <- function(all_data, name_chromosome_to_turn){
  #' Reverse the order of the positions along a specified chromosome
  #' @description This function determines the maximum values of a chromosome
  #' given its name and rearranges the positions of this chromosome by defining
  #' its last position as its first and its first position as its last.
  #' @param all_data This data frame contains all the coordinates of the aligned
  #' genomes on one another.
  #' @param name_chromosome_to_chose This string is the name that was given to the
  #' chromosome to reorder when it was assembled.
  #' 
  #' @returns A data frame similar to the input `all_data`, but with the
  #' positions of the focal chromosome reversed.
  
  # Get the maximum position of the chromosome
  max_qstart <- all_data %>% 
    filter(qname == name_chromosome_to_turn) %>% 
    pull(qstart) %>% 
    max(na.rm = TRUE)
  max_qend <- all_data %>% 
    filter(qname == name_chromosome_to_turn) %>% 
    pull(qend) %>% 
    max(na.rm = TRUE)
  max_tstart <- all_data %>% 
    filter(tname == name_chromosome_to_turn) %>% 
    pull(tstart) %>% 
    max(na.rm = TRUE)
  max_tend <- all_data %>% 
    filter(tname == name_chromosome_to_turn) %>% 
    pull(tend) %>% 
    max(na.rm = TRUE)
  
  
  all_data %>% 
    filter(qname == name_chromosome_to_turn | tname == name_chromosome_to_turn) %>% 
    mutate(across(c(qstart, tstart, qend, tend), ~ ., .names = "{col}_ori"),
           qend = ifelse(qname == name_chromosome_to_turn, max_qstart - qstart_ori, qend),
           qstart = ifelse(qname == name_chromosome_to_turn, max_qend - qend_ori, qstart),
           tend = ifelse(tname == name_chromosome_to_turn, max_tstart - tstart_ori, tend),
           tstart = ifelse(tname == name_chromosome_to_turn, max_tend - tend_ori, tstart)) %>% 
    select(-contains("ori")) %>% 
    rbind(all_data %>% 
            filter(qname != name_chromosome_to_turn & tname != name_chromosome_to_turn)) %>% 
    return()
}

is_chromosome_pair_reversed <- function(all_data, chromosome_1, chromosome_2){
  #' Checks if the pair of chromosome is reversed or not
  #' @description This function filters the data to keep the correspondences
  #' between the two specified chromosomes and runs a linear regression on the 
  #' corresponding positions in the two chromosomes to know if they are reversed
  #' @param all_data This data frame contains all the correspondences for all
  #' chromosomes.
  #' @param chromosome_1 This string is the name of the first chromosome
  #' @param chromosome_2 This string is the name of the second chromosome
  #' 
  #' @returns Boolean. TRUE if the slope of the linear regression is negative, i.e.
  #' if the chromosomes are reversed and FALSE if the chromosomes are well arranged.
  coefs_lm <- all_data %>% 
    filter(qname == chromosome_1, tname == chromosome_2) %>% 
    lm(data = ., qstart ~ tstart) %>% 
    coef() %>% 
    unname
  
  return(coefs_lm[2] < 0)
}

reverse_all_chromosomes <- function(all_alignments, correspondences_chromosomes){
  #' Reverse chromosomes that are not oriented in the same way in all the reference
  #' genomes
  #' @description
  #' This function goes through all the pairs of species and the pairs of
  #' chromosomes, checks if they are reversed compared to each other
  #' and reverses one of the two chromosomes if it is not already done.
  #' @param all_alignments This data frame is the output of the paf alignments 
  #' transformed into a table
  #' @param correspondences_chromosomes This data frame contains the correspondences
  #' of the chromosomes from one species to the other. It is the output of the 
  #' `get_corresponding_chromosomes_species` function.
  #' 
  #' @returns A similar data frame as the `all_alignments` input, but with 
  #' chromosomes aligned in the same direction in the different species.
  
  species_pairs <- all_alignments %>% 
    mutate(sp_pair = paste(query, target, sep = "_")) %>% 
    pull(sp_pair) %>% 
    unique
  
  already_reversed_chroms <- c()
  for (sp_pair in species_pairs){
    sps <- str_split_fixed(sp_pair, "_", 2) %>% 
      as.vector
    
    major_correspondences_sps_pair <- correspondences_chromosomes %>% 
      filter(query == sps[1], target == sps[2]) %>% 
      left_join(all_alignments,
                by = c("query", "target", "qname", "tname")) %>% 
      mutate(chrom_pair = paste(qname, tname, sep = "_"))
    for (chrom_pair in unique(major_correspondences_sps_pair$chrom_pair)){
      chrom_pair <- str_split_fixed(chrom_pair, "_", 2) %>% 
        as.vector()
      
      if (is_chromosome_pair_reversed(all_alignments, chrom_pair[1], chrom_pair[2])){
        chromosome_to_reverse <- which_not_in(chrom_pair, already_reversed_chroms)
        if (length(chromosome_to_reverse) > 1){
          chromosome_to_reverse <- chromosome_to_reverse[1]
        }
        all_alignments <- all_alignments %>% 
          reverse_chromosome(chromosome_to_reverse)
        already_reversed_chroms <- c(already_reversed_chroms,
                                     chromosome_to_reverse)
      }
    }
    
  }
  return(all_alignments)
}


compute_cumulative_positions_along_genome <- function(temp_data_to_plot, name_of_chromosomes){
  temp_data_to_plot %>% 
    left_join(name_of_chromosomes %>% 
                select(Chromosome, Species, Cumul_position_start_chromosome) %>% 
                unique,
              by = join_by("qname" == "Chromosome", "query" == "Species")) %>%
    mutate(qstart = qstart + Cumul_position_start_chromosome,
           qend = qend + Cumul_position_start_chromosome) %>% 
    select(-Cumul_position_start_chromosome) %>% 
    left_join(name_of_chromosomes %>% 
                select(Chromosome, Species, Cumul_position_start_chromosome),
              by = join_by("tname" == "Chromosome", "target" == "Species")) %>% 
    mutate(tstart = tstart + Cumul_position_start_chromosome,
           tend = tend + Cumul_position_start_chromosome) %>%  
    return()
}

rename_species_according_to_order <- function(temp_data_to_plot, order_species){
  
  temp_data_to_plot <- temp_data_to_plot %>% 
    mutate(nb_species = NA)
  
  for (sp in unique(temp_data_to_plot$Species)){
    number_of_species <- as.numeric(names(order_species)[which(sp == order_species)])
    temp_data_to_plot <- temp_data_to_plot %>% 
      mutate(nb_species = case_when(Species == sp ~ number_of_species,
                                    TRUE ~ nb_species))
  }
  return(temp_data_to_plot)
}

draw_ellipses_for_chromosome_delims <- function(names_of_chromosomes, order_species, reduction_length_proportion = 0.5){
  names_of_chroms <- names_of_chromosomes %>% 
    group_by(Species) %>% 
    arrange(Chromosome_name) %>% 
    mutate(Cumul_position_start_chromosome = cumsum(lag(Length + 1e7, default = 0)),
           Center_chrom = Cumul_position_start_chromosome + Length / 2,
           Length_ellipse = reduction_length_proportion * Length)
  
  
  data_ellipses <- data.frame()
  for (chromosome in unique(names_of_chroms$Chromosome_name)){
    for (species in unique(names_of_chroms$Species)){
      temp_df <- names_of_chroms %>% 
        filter(Chromosome_name == chromosome,
               Species == species) %>% 
        rename_species_according_to_order(order_species)
      if (nrow(temp_df) == 0) next
      temp_out_df <- getellipse(rx = temp_df$Length_ellipse, ry = 0.1, mid = c(temp_df$Center_chrom, temp_df$nb_species)) %>% 
        as_tibble() %>% 
        mutate(Chromosome_name = chromosome,
               Species = species,
               Group = paste0(Chromosome_name, Species)) %>% 
        rename(x = V1, y = V2)
      
      data_ellipses <- data_ellipses %>% 
        add_table_to_df_in_iteration(temp_out_df)
    }
  }
  return(data_ellipses)
}


make_synteny_plot <- function(alignments, name_of_chromosomes, 
                              only_chromosome = NULL,
                              small_alignment_threshold = 5e4,
                              order_species = c("praehirsuta", "albifrons", "ischiosetosa"),
                              file_name = NULL,
                              shape_chroms = "rectangle"){
  #' Traces a synteny plot on the provided data
  #' @description
    #' This function takes the alignments and the names of the chromosomes,
    #' computes the cumulative positions of the alignments along the genome
    #' and plots it. This function is mainly graphics
  #' @param alignments This table contains the alignments for the different species
  #' combinations. It is the equivalent of the `.paf` files, but transformed as tables
  #' @param name_of_chromosomes This table contains the names of the chromosomes
  #' for the different species. It is the output of the `get_chromosome_names_all_species`
  #' function
  #' @param only_chromosome This argument is a string or vector of strings that 
  #' contains the name of one chromosome if you want to look at it individually.
  #' It allows to filter the data to look only at this chromosome. The default
  #' value is `NULL`.
  #' @param small_alignment_threshold This number allows to filter small fragments
  #' that map. It defines the threshold of what the user specifies as "small". The 
  #' default value is `5e4`
  #' @param order_species This vector contains the names of the species that we
  #' want on the graph. It is important that the species be given in order from
  #' top to bottom of the graph. The default value is `c("praehirsuta", "albifrons", "ischiosetosa")`.
  #' However, this is to be changed with your dataset
  #' @param file_name This string contains the name of the file to be produced in
  #' your working directory. The default is `NULL`. In this case, it plots the graph
  #' in your plot viewer. As soon as a string is specified, it will save the plot
  #' without plotting it.
  #' @param shape_chroms This is a graphic parameter that allows to change the 
  #' shape that the chromosomes have according to your preferences. Its possible
  #' values are `"round" or "rectangle"`. It will draw chromosomes as ellipses or
  #' rectangles respectively.
  #' 
  #' @returns The synteny plot as a ggplot2 object
  
  # Prepare the order of the species along the y axis and prepare the y axis
  order_species <- rev(order_species)
  names(order_species) <- as.character(1:length(order_species))
  y_scale <- scale_y_discrete(labels = order_species)
  
  
  # Get the value of the gap to leave between chromosomes
  order_magnitude_length_genome <- name_of_chromosomes %>% 
    mutate(max_position_chromosome = Cumul_position_start_chromosome + Length) %>% 
    pull(max_position_chromosome) %>% 
    max(na.rm = TRUE) %>% 
    order_magnitude()
  
  gap_between_chroms <- 0.1 * 10 ^ order_magnitude_length_genome
  
  # Prepare the data for plotting
  data_to_plot <- alignments %>%
    # We filter the data to keep alignments larger than the specified threshold
    filter(alen > small_alignment_threshold) %>% 
    # Add the correspondence of the chromosomes for the query and target species
    left_join(name_of_chromosomes %>% 
                select(Chromosome, Chromosome_name),
              by = join_by("qname" == "Chromosome")) %>% 
    left_join(name_of_chromosomes %>% 
                select(Chromosome, Chromosome_name),
              by = join_by("tname" == "Chromosome"),
              suffix = c("_query", "_target"))
  
  # If we want to compute the data only for one chromosome, we filter the data
  if (!is.null(only_chromosome)){
    regex_filter_only_chrom <- paste0(only_chromosome, collapse = "|")
    data_to_plot <- data_to_plot %>% 
      filter(grepl(regex_filter_only_chrom, Chromosome_name_query),
             grepl(regex_filter_only_chrom, Chromosome_name_target))
    
    name_of_chromosomes_to_plot <- name_of_chromosomes %>% 
      filter(grepl(only_chromosome, Chromosome_name)) %>% 
      group_by(Species) %>% 
      mutate(Cumul_position_start_chromosome = cumsum(lag(Length + gap_between_chroms, default = 0)),
             Center_chrom = Cumul_position_start_chromosome + Length / 2) %>% 
      ungroup %>% 
      rename_species_according_to_order(order_species)
  }else{
    name_of_chromosomes_to_plot <- name_of_chromosomes %>% 
      group_by(Species) %>% 
      mutate(Cumul_position_start_chromosome = cumsum(lag(Length + gap_between_chroms, default = 0)),
             Center_chrom = Cumul_position_start_chromosome + Length / 2) %>% 
      ungroup %>%
      rename_species_according_to_order(order_species)
  }
  
  # Prepare the data for plotting
  data_to_plot <- data_to_plot %>%
    # First, compute the cumulative positions along the genome
    compute_cumulative_positions_along_genome(name_of_chromosomes_to_plot)
  
  # Then, we plot the data
  p <- data_to_plot %>%
    select(contains("Chromosome"), qstart, qend, tstart, tend, query, target) %>% 
    # Here, we prepare the data so it is plotable using polygons
    mutate(Position = row_number(),
           Chromosome_name_query = str_split_fixed(Chromosome_name_query, "\\.", 2)[, 1]) %>% 
    pivot_longer(cols = c(qstart, qend, tstart, tend), names_to = "name", values_to = "delims_polygon") %>% 
    pivot_longer(cols = c(query, target), names_to = "direction", values_to = "Species") %>% 
    filter((grepl("^q", name) & grepl("^q", direction)) | (grepl("^t", name) & grepl("^t", direction))) %>% 
    unique %>% 
    select(-Cumul_position_start_chromosome) %>%
    group_by(Position) %>% 
    # This also means that we have to rearrange the data in a special way
    arrange_for_polygons() %>% 
    rename_species_according_to_order(order_species) %>%
    # Start to plot
    ggplot() +
    geom_polygon(aes(x = delims_polygon, y = factor(nb_species), group = Position, fill = Chromosome_name_query),
                 alpha = 0.2)
  
  # Here, we add the chromosomes to the plot. You can choose if you prefere having them as ellipses or rectangles
  if (shape_chroms == "round"){
    data_ellipses <- name_of_chromosomes %>% 
      draw_ellipses_for_chromosome_delims(order_species) %>% 
      mutate(Chromosome_name = str_split_fixed(Chromosome_name, "\\.", 2)[, 1])
    
    p <- p +
      geom_polygon(data = data_ellipses, aes(x = x, y = y, group = Group, fill = Chromosome_name),
                   colour = "black", lwd = 1.2)
  }else{
    # On the graph, there will be some rectangles that will be used to delimitate the
    # chromosomes, but we have to build the data for these rectangles
    delimitations_chromosomes <- name_of_chromosomes_to_plot %>% 
      rename(start_chrom = Cumul_position_start_chromosome) %>% 
      mutate(End_chromosome = start_chrom + Length) %>% 
      rename_species_according_to_order(order_species) %>%
      mutate(Group_id = paste(Chromosome_name, Species, sep = "_"),
             Chromosome_name = str_split_fixed(Chromosome_name, "\\.", 2)[, 1]) 
    p <- p +
      geom_rect(data = delimitations_chromosomes,
                aes(xmin = start_chrom, xmax = End_chromosome,
                    ymin = nb_species - 0.05, ymax = nb_species + 0.05,
                    fill = Chromosome_name, group = Group_id),
                colour = "black", lwd = 1.2, linejoin = "round", lineend = "round")
  }
  
  # We finish the plot with the names of the chromosomes
  p <- p +
    geom_text(data = name_of_chromosomes_to_plot,
              aes(x = Center_chrom, y = factor(nb_species), label = Chromosome_name)) +
    labs(x = "Position along genome",
         y = "Species",
         fill = "Chromosome") +
    y_scale +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    guides(color = guide_legend(override.aes = list(alpha = 1)))
  
  
  # Choose if you want to save the figure or plot it
  if (is.null(file_name)){
    return(p)
  }else{
    ggsave(plot = p, file_name, scale = 3, width = 2000, height = 1500, units = "px")
    message(paste0('File saved as: "', file_name, '"'))
  }
}

arrange_for_polygons <- function(positions_to_plot){
  
  
  df_to_return <- data.frame()
  orientation <- "simple"
  for (species in unique(positions_to_plot$Species)){
    if (orientation == "simple"){
      temp_df <- positions_to_plot %>% 
        filter(Species == species) %>% 
        arrange(delims_polygon)
      orientation <- "opposite"
    }else{
      temp_df <- positions_to_plot %>% 
        filter(Species == species) %>% 
        arrange(desc(delims_polygon))
      orientation <- "simple"
    }
    df_to_return <- add_table_to_df_in_iteration(df_to_return,
                                                 temp_df)
    
    
  }
  
  return(df_to_return)
}
