import numpy as np
import scipy.stats as stats
#Input format: Samples as Rows, Genes as Columns and last column is the class/group label

class GeneSelection2:

    def __init__(self,features_mat):
        self.features_mat = features_mat

    def select_k(self,k):
        labels = set()
        for lab in self.features_mat[:,-1]:
            labels.add(lab)

        #iterate over each gene/col
        gene_f_vals = []
        for gene_col in range(self.features_mat.shape[1]-1):
            tmp_classes = []
            for label in labels:
                tmp_classes.append(self.features_mat[self.features_mat[:,-1]==label][:,gene_col])

            test = stats.f_oneway(*tmp_classes)
            gene_f_val = test[0]
            gene_f_vals.append((gene_col,gene_f_val))

        #sort according to f-value
        gene_f_vals.sort(key=lambda x: x[1], reverse=True)
        #take top k
        f_vals = gene_f_vals[0:k]
        #resort according to their gene/column
        genes_selected_by_f_test = sorted(f_vals, key=lambda x: x[0])
        #returns only the genes
        genes = [x[0] for x in genes_selected_by_f_test]
        return genes


#a = np.random.rand(50,1001)
#a[:,-1]= 10*a[:,-1]
#a[:,-1] = np.floor(a[:,-1])
#print a
#b = GeneSelection2(a)
#c = b.select_k(10)
#print c