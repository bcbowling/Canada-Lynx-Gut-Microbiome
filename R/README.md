This is a folder for all R code

## Clean Up 
Order for fastq file clean up and prep
1. install_packages.R (install any missing packages needed for clean up)
2. clean_fastq.R (clean up data and create sequence table)
3. assign_taxonomy.R (create taxonomy table)
4. phylo_tree.R (create phylogenetic tree)
5. create_phyloseq.R (combine data into phyloseq object for analysis)
6. remove_contaminants.r (remove contaminants from phyloseq object)

## Analysis (in Progress)
Steps one and two need to be done in order, the rest don't necessarily
1. install_analysis.r (install any missing packages needed for analysis)
2. taxonomy_filtering.R (filter out low prevalence sequences and group by genus)
3. diversity_analysis.r (calculate alpha and beta diversity metrics)
4. time_series_analysis.R (generate horizon plots with BiomeHorizon)
5. statistical_analyses.R (calculate statistics and significance)
    Question 1: Does overall community composition change through time?
    Step 1. Compute a community distance matrix: Bray–Curtis / Weighted UniFrac
    Step 2. (visualization in diversity_analysis.r)
    Step 3. PERMANOVA with ScatID as strata 
    Question 2: Does alpha diversity change through time?
    Step 1. mixed-effects models (Does diversity change with age of sample?)
    
microbiome numbers over time

microbiome species over time