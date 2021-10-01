import numpy as np
import pandas as pd
import sys
import logging

def isit_corr(C,eigval):
    psd = True
    if  np.linalg.norm(np.diag(C)-1)>1e-6:
        #print('STATUS: DIAGONALS ARE NOT 1!')
        psd = False
    if np.linalg.norm(.5*(C+C.T)-C)>1e-6:
        #print('STATUS: NOT SYMMETRIC!')
        psd = False
    if np.product(C<=1+1e-6) != 1 or np.product(C>=-1-1e-6)!=1:
        #print('STATUS: VALUES NOT BETWEEN NEGATIVE 1 AND 1!')
        psd = False
    if min(eigval) < 0:
        #print('STATUS: EIGENVALUES ARE NOT POSITIVE!')
        psd = False
    return psd

def PSD_reb(C,eigen_lt,eigen_replace):
   C_sym = .5*(C+C.T)
   eigval, eigvec = np.linalg.eigh(C_sym)
   eigval_cutoff = np.diag([eigen_replace if val<eigen_lt else val for val in eigval])

   diagonal_scaling = np.diag(1/np.einsum('ij,jj->i',eigvec**2,eigval_cutoff))
   decomposed_matrix = np.matmul(np.sqrt(diagonal_scaling),np.matmul(eigvec,np.sqrt(eigval_cutoff)))
   C_prime = np.matmul(decomposed_matrix,decomposed_matrix.T)

   return(C_prime)

def PSD_loop(corr,eigen_lt,eigen_replace,max_iter,final_reb):
    corr = .5*(corr+corr.T)
    for i in range(0,max_iter):
        eigval, eigvec = np.linalg.eigh(corr)
        if isit_corr(corr,eigval):
            logging.debug(f'Looping converged in {i} iterations.')
            return corr
  
        eigval_cutoff = np.diag([eigen_replace if val<eigen_lt else val for val in eigval])
        corr = np.matmul(eigvec,np.matmul(eigval_cutoff,eigvec.T))
        np.fill_diagonal(corr,1)
        corr = np.clip(corr,-1,1)

    if final_reb == 1:
        logging.debug(f'Looping didn\'t converge in {max_iter} iterations, applying Rebanato\'s method.')
        corr = PSD_reb(corr,eigen_lt=eigen_lt,eigen_replace=eigen_replace)

    return corr

def PSD_approx(file_name,work_path,name_col,output_table,eigen_lt,eigen_replace,max_iter,method,final_reb):
    # convert to numpy array
    df_input_matrix = pd.read_csv(work_path+'\\'+file_name+'.csv')

    # ensure our name_col has the right case
    col_list =  df_input_matrix.columns.to_list()
    name_col = col_list[[str(col).lower() for col in col_list].index(str.lower(name_col))]

    names = df_input_matrix.loc[pd.isna(df_input_matrix[name_col]) == False][name_col]
    input_matrix = df_input_matrix.loc[pd.isna(df_input_matrix[name_col]) == False, names].copy()
    input_matrix = input_matrix.fillna(value=0.0).to_numpy()

    # convert input_matrix to a correlation matrix if it isn't
    if  np.linalg.norm(np.diag(input_matrix)-1)>1e-6:
        logging.debug('You\'ve inputed a covariance matrix.')
        n = input_matrix.shape[0]
        var_list = np.array([np.sqrt(input_matrix[i,i]) for i in range(n)])
        corr = np.array([[input_matrix[i, j]/(var_list[i]*var_list[j]) for i in range(n)] for j in range(n)])
        is_input_corr = False
    else:
        logging.debug('You\'ve inputed a correlation matrix.')
        corr = input_matrix
        is_input_corr = True
    
    if method == 1:
        corr_out = PSD_loop(corr,eigen_lt=eigen_lt,eigen_replace=eigen_replace,max_iter=max_iter,final_reb=final_reb)
    else:
        corr_out = PSD_reb(corr,eigen_lt=eigen_lt,eigen_replace=eigen_replace)

    # convert back to input_matrixariance
    if is_input_corr:
        input_matrix_out = corr_out
    else:
        input_matrix_out = np.array([[corr_out[i, j]*(var_list[i]*var_list[j]) for i in range(n)] for j in range(n)])
    # convert back to input format
    df_input_matrix.loc[pd.isna(df_input_matrix[name_col]) == False,names] = input_matrix_out
    df_input_matrix.to_csv(work_path+'\\'+output_table+'.csv',index=False)

    return df_input_matrix



if __name__ == "__main__":
    input_file_name  = sys.argv[1]
    param_name       = sys.argv[2]
    work_path        = sys.argv[3]
    log_file         = sys.argv[4]

    try:
        logging.basicConfig(level=logging.DEBUG,filename=log_file,  format='%(levelname)s - %(message)s')
    except:
        logging.basicConfig(level=logging.DEBUG,filename=log_file, filemode='w', format='%(levelname)s - %(message)s')
    logging.debug('\n\nRUNNING POSDEF\n\n')

    try:
        params = pd.read_csv(work_path+'\\'+param_name+'.csv')
        # make naming conventions uniform
        params.columns = map(str.lower, params.columns)
        params = params.rename(columns=dict(
            outputtable  = 'output_table',
            namecol      = 'name_col',
            eigenlt      = 'eigen_lt',
            eigenreplace = 'eigen_replace',
            maxiter      = 'max_iter',
            method       = 'method',
            finalreb    = 'final_reb'
            ))
        
        params.max_iter=int(params.max_iter)
        params.method=int(params.method)
        params.final_reb=int(params.final_reb)

        logging.debug('System Arguments: \n'+str(sys.argv))
        logging.debug('Params: \n'+str(params))

        PSD_approx(input_file_name,work_path,**params.to_dict('r')[0])
    except Exception:
        logging.exception('Failed to run main PosDef function.')
    
# Example Run:
# python .\PosDefRunIt.py NewCorrMatrix PosDefParams .\tabular_data\ .\tabular_data\posdef.log
