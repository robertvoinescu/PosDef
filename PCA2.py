def pca_info(cov_mat, confid_alpha, n_pc=-1):
     #function to calculate PCA
    eig_vals, eig_vecs = np.linalg.eig(cov_mat)
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
        S_pca = S_pca + weights[i]*np.abs(np.dot(pca1[:,i],pca2[:,i]))/(np.linalg.norm(pca1[:,i])*np.linalg.norm(pca2[:,i]))

    return S_pca/len(weights)

#Step 1: calculate PCAs for both input (original) and realized (simulations)
orig_ev, orig_pca = pca_info(cov_mat_df, 0.95)
np_sim_ev, np_sim_pca = pca_info(np_sim_cor, 0.95)
#Step 2 : Calculate the coefficient of similarity:

np_weights = [1/len(orig_ev) for x in orig_ev]
S_pca_np = pca_similarity(orig_pca,np_sim_pca,np_weights) # pca similarity defined by abs
