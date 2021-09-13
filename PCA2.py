import numpy as np
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

def pca_info(cov_mat, confid_alpha, n_pc=-1):
     #function to calculate PCA
    eig_vals, eig_vecs = np.linalg.eigh(cov_mat)
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
        S_pca = S_pca + weights[i]*np.abs(np.dot(pca1[i,:],pca2[i,:]))/(np.linalg.norm(pca1[i,:])*np.linalg.norm(pca2[i,:]))
        print('SPCA: ',S_pca)
        print(pca1[i,:])
        print(pca2[i,:])
    return S_pca/len(weights)



def get_matrix(work_directory,table,name_col):
    df_input_matrix = pd.read_csv(work_directory+'\\'+table+'.csv')
    names = df_input_matrix.loc[pd.isna(df_input_matrix[name_col]) == False][name_col]
    input_matrix = df_input_matrix.loc[pd.isna(df_input_matrix[name_col]) == False, names].copy()
    input_matrix = input_matrix.fillna(value=0.0).to_numpy(dtype=np.float64)
    return input_matrix

if __name__ == "__main__":
    args = parse_args()
    #Step 1: calculate PCAs for both input (original) and realized (simulations)
    #A = get_matrix(args.work_directory,args.table1,args.name_col)
    #B = get_matrix(args.work_directory,args.table2,args.name_col)
    a = .01
    A = np.array([[1,0,0],[0,1,a],[0,a,1]])#+addM
    B = np.array([[1,a,0],[a,1,0],[0,0,1]])#+addM
    Aevec, Alam = np.linalg.eigh(A)
    Bevec, Blam = np.linalg.eigh(B)
    print('A')
    print(Aevec)
    print(Alam)
    print('B')
    print(Bevec)
    print(Blam)

    orig_ev, orig_pca       = pca_info(A, 1, n_pc=A.size)
    np_sim_ev, np_sim_pca   = pca_info(B, 1, n_pc=B.size)
    print('PCA MATRIX')
    print(orig_pca)
    print(np_sim_pca)
    
    #Step 2 : Calculate the coefficient of similarity:
    np_weights = [1 for x in orig_ev]
    S_pca_np = pca_similarity(orig_pca,np_sim_pca,np_weights) # pca similarity defined by abs
    print("Coefficient of Similarity: ",S_pca_np)

    normA = np.sqrt(np.trace(np.matmul(A,A)))
    normB = np.sqrt(np.trace(np.matmul(B,B)))
    overlapfrob = np.trace(np.matmul(A,B.T))/normA/normB
    print("Frobenius Overlap: ",overlapfrob)
