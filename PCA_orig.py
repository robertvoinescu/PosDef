import numpy as np
from numpy.linalg import eigh as eig
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--table1",         required=True, help="name of input csv table")
    parser.add_argument("--table2",         required=True, help="name of input csv table")
    parser.add_argument("--work_directory", required=True, help="directory to write and read from")
    parser.add_argument("--name_col",       required=True, help="directory to write and read from")
    args = parser.parse_args()
    return args

def get_matrix(work_directory,table,name_col):
    df_input_matrix = pd.read_csv(work_directory+'\\'+table+'.csv')
    names = df_input_matrix.loc[pd.isna(df_input_matrix[name_col]) == False][name_col]
    input_matrix = df_input_matrix.loc[pd.isna(df_input_matrix[name_col]) == False, names].copy()
    input_matrix = input_matrix.fillna(value=0.0).to_numpy(dtype=np.float64)
    return input_matrix

def pca_info(cov_mat, confid_alpha, n_pc=-1):
     #function to calculate PCA
    eig_vals, eig_vecs = eig(cov_mat)
    eig_vals_trunc = eig_vals.copy()
    eig_vals_trunc=eig_vals_trunc.real
    eig_vecs = eig_vecs.real
    eig_vals_trunc[eig_vals_trunc < 0] = 1e-20
    ev_explained_ratio = eig_vals_trunc / sum(eig_vals_trunc)
    ev_cum_explained = np.cumsum(ev_explained_ratio)
    if n_pc==-1: # if no number of PC is given, select PC based on confid_alpha
        ev_values=eig_vals_trunc[ev_cum_explained<confid_alpha,]
    else: # select the first n_pc pricipal components
        ev_values=eig_vals_trunc[:n_pc]
    
    eig_vecs = eig_vecs[:,:ev_values.shape[0]]
    pca = np.dot(eig_vecs, np.diag(np.sqrt(ev_values)))

    return [ev_values, pca]

def pca_similarity(pca1,pca2,weights):
#main function to calculate similarity. For identical matrices it should produce 1. Anything less than 1 indicate that matrices are different 
    S_pca = 0
    for i in range(len(weights)):
        #S_pca = S_pca + weights[i]*np.abs(np.dot(pca1[:,i],pca2[:,i]))/(np.linalg.norm(pca1[:,i])*np.linalg.norm(pca2[:,i])) 
        S_pca = S_pca + weights[i]*np.abs(np.dot(pca1[:,i],pca2[:,i]))
    return S_pca
    #return S_pca/len(weights)

if __name__ == "__main__":
    #args = parse_args()
    #A = get_matrix(args.work_directory,args.table1,args.name_col)
    #B = get_matrix(args.work_directory,args.table2,args.name_col)

    print('BEGIN CALCULATION')
    print('#'*30)
    a = 0.01
    A = np.array([[1,0,0],[0,1,a],[0,a,1]])#+addM
    B = np.array([[1,a,0],[a,1,0],[0,0,1]])#+addM
    Alam, Aevec = eig(A)
    Blam, Bevec = eig(B)
    print('BEGIN CALCULATION')
    print('-'*60)
    print('\nA Matrix')
    print(A)
    print('\nEigenvalues for A')
    print(np.diag(Alam))
    print('\nEigenvectors for A')
    print(Aevec)

    print('-'*60)

    print('\nB Matrix')
    print(B)
    print('\nEigenvalues for B')
    print(np.diag(Blam))
    print('\nEigenvectors for B')
    print(Bevec)

    normA = np.sqrt(np.trace(np.matmul(A,A)))

    #Step 1: calculate PCAs for both input (original) and realized (simulations)
    A_ev, A_pca = pca_info(A, 1, n_pc=A.size)
    B_ev, B_pca = pca_info(B, 1, n_pc=B.size)
    print('-'*60)
    print('\nPCA MATRIX A')
    print(A_pca)
    print(np.diag(A_ev))
    print('\nPCA MATRIX B')
    print(B_pca)
    print(np.diag(B_ev))
    
    #Step 2 : Calculate the coefficient of similarity:
    weights = [1/len(A) for x in A]
    coefficient_of_similarity = pca_similarity(A_pca,B_pca,weights) # pca similarity defined by abs

    print('#'*30)
    print(np.matmul(A_pca.T,B_pca))
    print('#'*30)
    print('-'*60)
    print("Coefficient of Similarity: ",coefficient_of_similarity)

# python .\PCA_orig.py --name_col x --table1 outmat --table2 outmat --work_directory .\tabular_data\
