import numpy as np
import pandas as pd
import argparse
from numpy.linalg import eigh as eig

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

if __name__=="__main__":
    #args = parse_args()
    # convert csv to numpy array
    #A = get_matrix(args.work_directory,args.table1,args.name_col)
    #B = get_matrix(args.work_directory,args.table2,args.name_col)

    a = 1
    n = 100
    A = np.eye(n,n,dtype=np.float64)
    B = np.eye(n,n,dtype=np.float64)
    C = np.random.rand(n,n)
    C = .5*abs(C-C.T)*1e-15
    print(C)
    #A[1,0], A[0,1] = a, a
    #B[2,1], B[1,2] = a, a
    B += C

    Alam, Aevec = eig(A)
    Blam, Bevec = eig(B)
    assert np.linalg.norm(np.matmul(np.matmul(Aevec,np.diag(Alam)),Aevec.T)-A) < 1e-10, 'Eigenvectors are rows'
    #assert np.linalg.norm(np.matmul(np.matmul(Aevec.T,np.diag(Alam)),Aevec)-A) < 1e-10, 'Eigenvectors are Rows'

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
    normB = np.sqrt(np.trace(np.matmul(B,B)))
    overlapfrob = np.trace(np.matmul(A,B.T))/normA/normB

    ###############
    num_pc = n-1
    k = Aevec.shape[1]-num_pc
    Ap = Aevec[:,k:]
    Bp = Bevec[:,k:]
    overlappaper = np.trace(np.matmul(np.matmul(Ap.T,Bp),np.matmul(Bp.T,Ap)))/Ap.shape[1]
    ###############
    overlapnrg = np.trace(abs(np.matmul(Ap.T,Bp)))/Ap.shape[1]

    print('-'*60)
    print(Ap)
    print(Bp)
    print('Coefficient of Similarity (NRG): ',overlapnrg)
    print('Coefficient of Simillarity (Paper): ',overlappaper)
    print('Frobenius Overlap: ',overlapfrob)
    print('-'*60)

# python .\PosDefRunIt_mat.py --max_iter 1000 --name_col x --input_table paper --output_table out --work_directory .\tabular_data\
