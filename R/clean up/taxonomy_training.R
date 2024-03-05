# Using Ribosomal Database Project training set 19
dada2:::makeTaxonomyFasta_RDP(file.path("data", "raw",
                                        "RDPTraining",
                                        "trainset19_072023.fa"),
                              file.path("data", "raw",
                                        "RDPTraining",
                                        "trainset19_db_taxid.txt"),
                              "data/modified/rdp_train_set_19.fa.gz")