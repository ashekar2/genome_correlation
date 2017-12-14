import numpy as np
import scipy.stats as stats
#This method normalizes a matrix of gene expresson data to be between 0 and 1
#Input format: Samples as Rows, Genes as Columns and last column is the class/group label



class Group():
    def __init__(self,group_id):
        self.group_id = group_id
        self.rows = []

    def add_row(self,item):
        self.rows.append(item)
        sorted(self.rows)

    def get_group_submatrix(self,original_mat):
        self.group_mat = original_mat[sorted(self.rows)]

    def get_size(self):
        self.size = len(self.rows)


#This class performs Gene Selection using different methods from the literature
class GeneSelection:
    def __init__(self,Gene_exp_mat,k):
        self.GeneExpMat = Gene_exp_mat
        self.main()
        self.select_k(k)


    #Use Fisher Score/ F-Test and take the top ~100
    #Source: Gene selection for cancer classification with the help of bees

    def main(self):

        class_hashmap = dict()
        self.gene_f_vals = []
        self.gene_kruskal_vals = []
        #goes thru all columns since each one is a gene
        for gene_Col in range(self.GeneExpMat.shape[1]-1): # except last one since it is the labels for each row
            #each row is a sample
            for row in range (self.GeneExpMat.shape[0]):
                #add the class if it is not already in the set
                if self.GeneExpMat[row,-1] not in class_hashmap.keys():
                    class_hashmap[self.GeneExpMat[row,-1]] = Group(self.GeneExpMat[row,-1])
                    class_hashmap[self.GeneExpMat[row,-1]].add_row(row)

                else:
                    class_hashmap[self.GeneExpMat[row, -1]].add_row(row)

            #iterates over each class
            class_list = []
            for key in class_hashmap.keys():
                class_hashmap[key].get_group_submatrix(self.GeneExpMat[:,gene_Col])
                class_list.append(class_hashmap[key].group_mat)

            test = stats.f_oneway(*class_list)
            gene_f_val = test[0]
            self.gene_f_vals.append((gene_Col,gene_f_val))

            test2 = stats.kruskal(*class_list)
            gene_kruskal_stat = test2[0]
            self.gene_kruskal_vals.append((gene_Col,gene_kruskal_stat))





    def select_k(self,k):
        self.gene_f_vals.sort(key=lambda x: x[1],reverse=True)
        print self.gene_f_vals
        self.gene_kruskal_vals.sort(key=lambda x:x[1],reverse=True)
        f_vals = self.gene_f_vals[0:k]
        kruskal_vals = self.gene_kruskal_vals[0:k]
        self.genes_selected_by_f_test = sorted(f_vals,key=lambda x:x[0])
        self.genes_selected_by_kruskal = sorted(kruskal_vals,key=lambda x:x[0])

        self.genes = [x[0] for x in self.genes_selected_by_f_test]






#run it like:
# foo = GeneSelection('data_set_ALL_AML_independent.csv',20) # for k = 20
#foo.genes