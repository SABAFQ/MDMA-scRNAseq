import pickle
import numpy as np 
import pandas as pd

def save_obj(obj, name ):
    with open(name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

        
def load_obj(name ):
    with open(name, 'rb') as f:
        return pickle.load(f)

    
<<<<<<< HEAD
=======
def normalization(adata,lognorm=False):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    if lognorm: 
        sc.pp.log1p(adata)


def clustering(adata,lognorm=False,rep=None, n_neighbor=30,n_comps=50,res=None,key="leiden"):
    sc.pp.pca(adata,n_comps=n_comps)
    sc.pp.neighbors(adata, use_rep=rep, n_neighbors=n_neighbor)   
    # why can't we just work with the default neighbors?
    sc.tl.umap(adata)
    sc.tl.leiden(adata,resolution=res,key_added=key)


def rank_genes_to_df(adata,n=50):
    result = adata.uns['rank_genes_groups']

    groups = result['names'].dtype.names

    df = pd.DataFrame(
        {group + '_' + key: result[key][group]
        for group in groups for key in ['names', 'scores']}).head(n)
    
    return df


>>>>>>> bd825b58299e6e3f64c31eec357c17b1f40fe81c
def find_top(df,value, value_thr, stat, stat_thr, drop_dup=False,silent=False):
    """
    Select rows (genes) which has value >= value_thr & stat < stat_thr 
    """
    up = df.iloc[
        [i for i,l in enumerate(
            np.array([
                np.array(df.loc[:,value] >= value_thr),
                np.array(df.loc[:,stat] < stat_thr)]).all(axis=0)) if l == 1]
            ,:] 
    dn = df.iloc[
        [i for i,l in enumerate(
            np.array([
                np.array(df.loc[:,value] <= -1*(value_thr)),
                np.array(df.loc[:,stat] < stat_thr)]).all(axis=0)) if l == 1]
            ,:] 
        
    if drop_dup==True:
        up = up.sort_values(stat).drop_duplicates(subset='gene_id', keep="last")
        dn = dn.sort_values(stat).drop_duplicates(subset='gene_id', keep="last")
    
    if not silent:
        print ('up: ', up.shape[0])
        print ('down:', dn.shape[0])

    return up, dn
