import numpy as np



#this method takes selects k/2 genes from the f-score and k/2 from the DPP
def fs_dpp_select_k(k,training_set,cov_mat):


    f_score_genes = gene_selection(training_set).select_k(np.floor(k/2))
    f_score_genes = set(f_score_genes)
    gene_indices = [i for i in range(training_set.shape[1])]

    dpp_k_features = dpp.sample_k(gene_indices,cov_mat,k, max_nb_iterations=20000)


    for gene in dpp_k_features:
            #if we still don't have k genes, then add genes into
            if gene not in f_score_genes and len(f_score_genes) != k:
                f_score_genes.add(gene)

            #we have reached k genes thus we can stop
            elif len(f_score_genes) == k:
                break


    return list(f_score_genes).sort()



