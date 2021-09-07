import numpy as np
import argparse
from numpy.linalg.linalg import eig
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max_iter",        required=True, help="number of projection steps",type=int)
    parser.add_argument("--input_table",     required=True, help="name of input csv table")
    parser.add_argument("--output_table",    required=True, help="name of output csv table")
    parser.add_argument("--work_directory",  required=True, help="directory to write and read from")
    parser.add_argument("--name_col",        required=True, help="directory to write and read from")
    args = parser.parse_args()
    return args

def isit_corr(C):
    # temporarily set this to false and to see what happens
    psd = False 
    eigval = np.linalg.eigvals(C)
    if  np.linalg.norm(np.diag(C)-1)>1e-6:
        #print('STATUS: DIAGONALS ARE NOT 1!')
        psd = False
    if np.linalg.norm(.5*(C+C.T)-C)>1e-6:
        #print('STATUS: NOT SYMMETRIC!')
        psd = False
    if np.product(C<=1+1e-6) != 1 or np.product(C>=-1-1e-6)!=1:
        #print('STATUS: VALUES NOT BETWEEN NEGATIVE 1 AND 1!')
        psd = False
    if np.product(eigval>0) == 0:
        #print('STATUS: EIGENVALUES ARE NOT POSITIVE!')
        print(min(eigval))
        psd = False
    if psd:
        print('It\'s a correlation matrix!')
    return psd

def matrix_func(lam, evec, f):
    lam = np.diag([f(x) for x in lam])
    A = np.matmul(np.matmul(evec,lam),evec.T)
    return A

def plus_spectral(A):
    lam, evec = np.linalg.eigh(A)
    return matrix_func(lam, evec, lambda x: max(x,0))

def Pu(A):
    proj = 1*A
    np.fill_diagonal(proj,1)
    return proj


def near_corr(A,W,max_iter):
    S = 0*A
    Y = 1*A
    for i in range(max_iter):
        if isit_corr(Y):
            print(f'Converged in {i} iterations.')
            break
        print(f'iteration: {i}')
        R = Y - S
        X = plus_spectral(R)
        S = X - R
        Y = Pu(X)
    return Y


if __name__=="__main__":
    args = parse_args()
    # convert csv to numpy array
    df_input_matrix = pd.read_csv(args.work_directory+'\\'+args.input_table+'.csv')
    names = df_input_matrix.loc[pd.isna(df_input_matrix[args.name_col]) == False][args.name_col]
    input_matrix = df_input_matrix.loc[pd.isna(df_input_matrix[args.name_col]) == False, names].copy()
    input_matrix = input_matrix.fillna(value=0.0).to_numpy(dtype=np.float64)

    # convert input_matrix to a correlation matrix if it isn't
    if  np.linalg.norm(np.diag(input_matrix)-1)>1e-6:
        print('You\'ve inputed a covariance matrix.')
        n = input_matrix.shape[0]
        var_list = np.array([np.sqrt(input_matrix[i,i]) for i in range(n)])
        corr = np.array([[input_matrix[i, j]/(var_list[i]*var_list[j]) for i in range(n)] for j in range(n)])
        is_input_corr = False
    else:
        print('You\'ve inputed a correlation matrix.')
        corr = input_matrix
        is_input_corr = True
 
    corr_out = near_corr(corr,np.eye(len(corr)),args.max_iter)

    # convert back to input_matrix
    if is_input_corr:
        input_matrix_out = corr_out
    else:
        input_matrix_out = np.array([[corr_out[i, j]*(var_list[i]*var_list[j]) for i in range(n)] for j in range(n)])
    # convert back to input format
    df_input_matrix.loc[pd.isna(df_input_matrix[args.name_col]) == False,names] = input_matrix_out
    df_input_matrix.to_csv(args.work_directory+'\\'+args.output_table+'.csv',index=False)
    
# Example Run: python simphigh.py --max_iter 1000 --name_col x --input_table pypos --output_table out --work_directory C:\Users\rvo67\Desktop\test_posdef\new_code
