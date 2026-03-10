
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
    summarise(total_length = sum(alen)) %>% 
    arrange(desc(total_length)) %>% 
    dplyr::slice(1) %>% 
    ungroup() %>% 
    select(-total_length) %>% 
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
              by = join_by("tname" == "Chromosome")) %>% 
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
      if (!any(sps %in% already_done_pair_of_species)){
        to_do_later <- c(to_do_later, pair_of_species)
        next
      }
      # Otherwise, we use one of them as the new reference species
      reference_species <- get_new_reference_species(sps, already_done_pairs_of_species)
      reference_chromosomes <- named_chromosomes_all_species %>% 
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


get_chromosome_names_all_species <- function(all_chromosomes, correspondences_chromosomes, reference_chromosomes = NULL){
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
  named_chromosomes_all_species <- data.frame()
  
  while (!is.null(pairs_species)){
    
    list_outputs_iteration <- iterate_over_pairs_of_species_to_get_name_of_chromosomes(
      pair_of_species_to_do = pairs_species,
      already_done_pairs_of_species = already_done_species,
      named_chromosomes_all_sps = named_chromosomes_all_species,
      correspondences_chromosomes = correspondences_chromosomes,
      reference_chromosomes = reference_chromosomes,
      reference_species = reference_species
    )
    
    pairs_species <- list_outputs_iteration$to_do_later
    already_done_species <- list_outputs_iteration$already_done
    named_chromosomes_all_species <- list_outputs_iteration$named_chromosomes_all_sps
    
    if ((!is.null(reference_species)) && (!any_string_in_vector_contains_pattern(reference_species, pairs_species))){
      reference_species <- NULL
      reference_chromosomes <- NULL
    }
  }
  
  named_chromosomes_all_species %>% 
    left_join(all_chromosomes,
              by = c("Chromosome", "Species")) %>% 
    return()
}
