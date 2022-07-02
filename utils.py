import pickle
import scanpy as sc


def save_obj(obj, name ):
    with open(name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

        
def load_obj(name ):
    with open(name, 'rb') as f:
        return pickle.load(f)

    
def clustering(adata,lognorm=False,rep=None, n_neighbor=30,n_comps=50,res=None):
    if lognorm: 
        sc.pp.log1p(adata)
    sc.pp.pca(adata,n_comps=n_comps)
    sc.pp.neighbors(adata, use_rep=rep, n_neighbors=n_neighbor)   # why can't we just work with the default neighbors?
    sc.tl.umap(adata)
    sc.tl.leiden(adata,resolution=res)