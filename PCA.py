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

def get_matrix(work_directory,table,name_col):
    df_input_matrix = pd.read_csv(work_directory+'\\'+table+'.csv')
    names = df_input_matrix.loc[pd.isna(df_input_matrix[name_col]) == False][name_col]
    input_matrix = df_input_matrix.loc[pd.isna(df_input_matrix[name_col]) == False, names].copy()
    input_matrix = input_matrix.fillna(value=0.0).to_numpy(dtype=np.float64)
    return input_matrix

if __name__=="__main__":
    args = parse_args()
    # convert csv to numpy array
    A = get_matrix(args.work_directory,args.table1,args.name_col)
    B = get_matrix(args.work_directory,args.table2,args.name_col)

    #a = 0.001
    #b = 1 
    #addM = np.array([[0,b,b],[b,0,b],[b,b,0]])
    #A = np.array([[1,0,0],[0,1,a],[0,a,1]])#+addM
    #B = np.array([[1,a,0],[a,1,0],[0,0,1]])#+addM

    #A = np.array([[1,0,0],[0,1,0],[0,0,1]])
    #B = np.array([[1,b,b],[b,1,b],[b,b,1]])
    Alam, Aevec = np.linalg.eigh(A)
    Blam, Bevec = np.linalg.eigh(B)
    #print(A)
    #print(B)
    #print('#'*30)
    #print(f'Eigenvalues A: {Alam}')
    #print(Aevec)
    #print(f'Eigenvalues B: {Blam}')
    #print(Bevec)

    normA = np.sqrt(np.trace(np.matmul(A,A)))
    normB = np.sqrt(np.trace(np.matmul(B,B)))
    overlapfrob = np.trace(np.matmul(A,B.T))/normA/normB
    overlaporg = abs(np.trace(np.matmul(Aevec,Bevec.T))/Aevec.shape[0])
    print('OVERLAP ORG: ',overlaporg)
    print('OVERLAP FROB: ',overlapfrob)

# python .\PosDefRunIt_mat.py --max_iter 1000 --name_col x --input_table paper --output_table out --work_directory .\tabular_data\
