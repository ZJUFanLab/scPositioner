import pandas as pd
import numpy as np
import random
import ot
import scipy
import scanpy as sc
from scanpy import AnnData
from scipy.sparse import issparse
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.neighbors import NearestNeighbors
# from tqdm import tqdm
from datetime import datetime
from typing import Optional, List


def set_seed(seed: int = 0):
    """
    Set random seed.

    Args:
        seed (int, optional): Random seed. Defaults to 0.

    Returns:
        None
    """
    random.seed(seed)
    np.random.seed(seed)


def filtering(adata: AnnData):
    """
    Filter cells and genes.

    Args:
        adata (AnnData): The AnnData object.

    Returns:
        AnnData: The filtered AnnData object.
    """
    sc.pp.filter_cells(adata, min_genes=1)
    sc.pp.filter_genes(adata, min_cells=1)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.obs['raw_id'] = adata.obs_names
    return adata


def normalize_log(adata: AnnData):
    """
    Normalize and log-transform the data.

    Args:
        adata (AnnData): The AnnData object.

    Returns:
        AnnData: The normalized AnnData object.
    """
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    return adata


def estimate_cell_number(adata: AnnData, mean_cell_numbers: int = 5, normalize: bool = True):
    """
    Estimate the number of cells per spot.
    Linear fitting based on the UMIs.
    From CytoSpace: https://github.com/digitalcytometry/cytospace/blob/main/cytospace/cytospace.py

    Args:
        adata (AnnData): The AnnData object.
        mean_cell_numbers (int, optional): The average number of cells across all spots. Defaults to 5.
        normalize (bool, optional): Whether to normalize data. Defaults to True.

    Returns:
        np.ndarray: The estimated number of cells per spot.
    """
    if mean_cell_numbers < 0:
        raise ValueError('mean_cell_numbers must be non-negative!')
    adata_use = adata.copy()
    if normalize:
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Estimating cell numbers per spot based on normalized data...")
        adata_use = normalize_log(adata_use)

    if issparse(adata_use.X):
        adata_use.X = adata_use.X.toarray()

    expr = adata_use.X.T.astype(float)

    # Set up fitting problem
    RNA_reads = np.sum(expr, axis=0, dtype=float)
    mean_RNA_reads = np.mean(RNA_reads)
    min_RNA_reads = np.min(RNA_reads)

    min_cell_numbers = 1 if min_RNA_reads > 0 else 0

    x = np.array([min_RNA_reads, mean_RNA_reads])
    y = np.array([min_cell_numbers, mean_cell_numbers])
    fit_parameters = np.polyfit(x, y, 1)
    polynomial = np.poly1d(fit_parameters)
    estimated_cell_number = np.round(polynomial(RNA_reads)).astype(int)

    return estimated_cell_number


def run_cell2location(
        adata_sc: AnnData,
        adata_st: AnnData,
        cell_num_per_spot: int = 5,
        celltype_key: Optional[str] = 'Cell_type',
        use_gpu: Optional[bool] = False,
        sc_epochs: Optional[int] = 1000,
        st_epochs: Optional[int] = 30000):
    
    #import cell2location
    import cell2location
    from cell2location.models import RegressionModel

    # sc
    # prepare anndata for the regression model
    cell2location.models.RegressionModel.setup_anndata(adata=adata_sc, labels_key=celltype_key)
    # create the regression model
    mod = RegressionModel(adata_sc)
    mod.train(max_epochs=sc_epochs, use_gpu=use_gpu)

    adata_sc = mod.export_posterior(
        adata_sc, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': use_gpu}
    )
    # export estimated expression in each cluster
    if 'means_per_cluster_mu_fg' in adata_sc.varm.keys():
        inf_aver = adata_sc.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                                             for i in adata_sc.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_sc.var[[f'means_per_cluster_mu_fg_{i}'
                                 for i in adata_sc.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_sc.uns['mod']['factor_names']

    # st
    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(adata=adata_st)
    # create and train the model
    # cell_num_per_spot = np.round(np.mean(adata_st.obs['Estimate_cell_num'])).astype(int)
    mod = cell2location.models.Cell2location(
        adata_st,
        cell_state_df=inf_aver,
        N_cells_per_location=cell_num_per_spot,
        detection_alpha=20
    )
    mod.train(max_epochs=st_epochs, batch_size=None, train_size=1, use_gpu=use_gpu)

    adata_st = mod.export_posterior(
        adata_st, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': use_gpu}
    )

    res = adata_st.obsm['q05_cell_abundance_w_sf']
    # normalize
    res = res.div(res.sum(axis=1), axis='rows')
    column_name = res.columns.tolist()
    column_name = [column_name[i].replace('q05cell_abundance_w_sf_', '') for i in range(len(column_name))]
    res.columns = column_name

    return res


def check_deconvolution_results(S: AnnData, R: AnnData, deconv_res: pd.DataFrame, celltype_key: str):
    """
    Check and adjust the format of the deconvolution results.

    Args:
        S (AnnData): Single cell AnnData object.
        R (AnnData): Spatial reference AnnData object.
        deconv_res (pd.DataFrame): Deconvolution results of the spatial reference, shape = (spot, celltype).
        celltype_key (str): The column name of `cell type` in `S.obs`. Defaults to 'celltype'.

    Returns:
        pd.DataFrame: The adjusted deconvolution results.
    """
    # deconv_res -> spot * celltype
    spot_d = list(deconv_res.index)
    celltype_d = list(deconv_res.columns)

    celltype_order = sorted(list(set(S.obs[celltype_key])))
    if sorted(celltype_d) != sorted(celltype_order):
        raise ValueError('Cell types in deconvolution results are not consistent with the input single cell AnnData!')
    if sorted(spot_d) != sorted(R.obs_names):
        raise ValueError('Spots in deconvolution results are not consistent with the input spatial AnnData!')
    deconv_res = deconv_res[celltype_order]
    deconv_res = deconv_res.loc[R.obs_names, :]

    return deconv_res


def integer_allocation(prop: np.ndarray, counts: np.ndarray):
    """
    Integer allocation for each spot.
    Ensure that the total number of cells obtained based on the deconvolution results matches the estimated number of cells.

    Args:
        prop (np.ndarray): Deconvolution results of the spatial reference.
        counts (np.ndarray): The estimated number of cells per spot.

    Returns:
        np.ndarray: Statistical results of the cell type composition within each spot, shape = (spot, celltype).
    """
    expected_cells = prop * counts[:, None]
    int_cells = np.floor(expected_cells).astype(int)
    frac_cells = expected_cells - int_cells
    remaining_cells = counts - int_cells.sum(axis=1)
    for i in range(len(counts)):
        frac_order = np.argsort(-frac_cells[i])
        for j in range(remaining_cells[i]):
            int_cells[i, frac_order[j]] += 1
    return int_cells


def adjust_abundance(S: AnnData, R: AnnData, celltype_key: str, seed: int = 0):
    """
    Adjust cell aboundance in single cell AnnData

    Args:
        S (AnnData): Single cell AnnData object.
        R (AnnData): Spatial reference AnnData object.
        celltype_key (str): The column name of `cell type` in `S.obs`.
        seed (int, optional): Random seed. Defaults to 0.

    Returns:
        AnnData: The adjusted single cell AnnData object.
    """
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Adjusting cell aboundance in single cell AnnData...")
    celltype_order = sorted(list(set(S.obs[celltype_key])))
    celltype_order2 = ['deconv_' + str(i) for i in celltype_order]
    deconv_res = np.array(R.obs[celltype_order2].copy())

    estimated_cell_number = np.array(R.obs['estimated_cell_number'])
    map_target = integer_allocation(deconv_res, estimated_cell_number)

    target_number_list = map_target.sum(axis=0)
    original_number_list = np.array(S.obs[celltype_key].value_counts()[celltype_order])
    diff_number_list = original_number_list - target_number_list

    adata_list = []
    for i in range(len(celltype_order)):
        adata_tmp = S[S.obs[celltype_key] == celltype_order[i]].copy()
        adata_list.append(adata_tmp)

    for i in range(len(celltype_order)):
        target_number = target_number_list[i]
        original_number = original_number_list[i]
        adjust_number = diff_number_list[i]
        adata_tmp = adata_list[i].copy()
        if adjust_number >= 0:
            # set seed
            set_seed(seed)
            selected_idx = np.random.choice(adata_tmp.shape[0], size=target_number, replace=False)
            adata_tmp = adata_tmp[selected_idx]
        elif adjust_number < 0:
            fold = np.abs(adjust_number) / original_number
            if fold > 1:
                fold_int = int(fold)
                selected_idx = list(range(adata_tmp.shape[0])) * fold_int

                set_seed(seed)
                selected_idx2 = list(
                    np.random.choice(
                        adata_tmp.shape[0], size=(np.abs(adjust_number) - fold_int * adata_tmp.shape[0]), replace=False
                        )
                    )
                selected_idx.extend(selected_idx2)
                selected_idx = np.array(selected_idx)
            else:
                set_seed(seed)
                selected_idx = np.random.choice(adata_tmp.shape[0], size=np.abs(adjust_number), replace=False)

            adata_tmp_replicate = adata_tmp[selected_idx].copy()
            adata_tmp = sc.concat([adata_tmp, adata_tmp_replicate]).copy()

        adata_list[i] = adata_tmp.copy()

    S_new = sc.concat(adata_list).copy()
    S_new.obs_names_make_unique()

    return S_new


def simulate_gene_dropout(sc_adata, dropout_rate, seed=0):
    set_seed(seed)
    
    n_genes = sc_adata.n_vars
    n_drop = int(round(n_genes * dropout_rate))

    gene_names = sc_adata.var_names.to_numpy()
    
    dropped_genes = np.random.choice(
        gene_names,
        size=n_drop,
        replace=False
    )

    kept_mask = ~sc_adata.var_names.isin(dropped_genes)
    new_adata = sc_adata[:, kept_mask].copy()

    return new_adata


def run_OT(S: AnnData, R: AnnData, numItermax: int = 1e6):
    """
    Optimal transport.

    Args:
        S (AnnData): Single cell AnnData object.
        R (AnnData): Spatial reference AnnData object.
        numItermax (int, optional): The maximum number of iterations before stopping the optimization algorithm if it has not converged. Defaults to 1e6.

    Returns:
        np.ndarray: The transport matrix.
    """
    Xs = S.X.copy()
    Xt = R.X.copy()

    if issparse(Xs):
        Xs = Xs.toarray()
    if issparse(Xt):
        Xt = Xt.toarray()

    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Calculating cost matrix (cosine similarity)...")
    cosine_sim_matrix = cosine_similarity(Xs, Xt)
    M = 1 - cosine_sim_matrix
    M /= M.max()

    # weights
    a = np.ones((S.shape[0],))
    b = np.array(R.obs['estimated_cell_number'])

    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Running optimal transport...")
    transport_matrix = ot.emd(a, b, M, numItermax=numItermax)

    return transport_matrix


def jitter_coord(coord: np.ndarray, seed: int = 0):
    """
    Jitter the coordinates

    Args:
        coord (np.ndarray): spatial coordinates of all spots.
        seed (int, optional): Random seed. Defaults to 0.

    Returns:
        np.ndarray: The jittered spatial coordinates.
    """
    set_seed(seed)
    num = coord.shape[0]
    
    if num == 1:
        tiny_noise = np.random.normal(0, 1e-8, size=(1, 2))
        return coord + tiny_noise
    
    # min distance
    coord_unique = np.unique(coord, axis=0)
    n_unique = coord_unique.shape[0]
    
    if n_unique < 2:
        tiny_noise = np.random.normal(0, 1e-8, size=(num, 2))
        return coord + tiny_noise
    
    n_neighbors = min(2, n_unique)
    nbrs = NearestNeighbors(n_neighbors=n_neighbors).fit(coord_unique)
    distances, _ = nbrs.kneighbors(coord_unique)

    if n_neighbors == 1:
        min_distance = 1.0
    else:
        positive_distances = distances[:, -1][distances[:, -1] > 0]
        if len(positive_distances) == 0:
            min_distance = 1.0
        else:
            min_distance = min(positive_distances)

    x_list = list(coord[:, 0])
    y_list = list(coord[:, 1])

    length = np.random.uniform(0, min_distance, num)
    radius = np.pi * np.random.uniform(0, 2, num)

    x_list_new = x_list + length * np.cos(radius)
    y_list_new = y_list + length * np.sin(radius)
    coord_new = np.array([[x_list_new[i], y_list_new[i]] for i in range(num)])

    return coord_new


def post_process(S: AnnData, R: AnnData, transport_matrix: np.ndarray, celltype_key: str):
    """
    Assign single cells based on the transport matrix and create new AnnData object.

    Args:
        S (AnnData): Single cell AnnData object.
        R (AnnData): Spatial reference AnnData object.
        transport_matrix (np.ndarray): The transport matrix.
        celltype_key (str): The column name of `cell type` in `S.obs`.

    Returns:
        AnnData: The mapped AnnData object.
    """
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Assigning cells...")
    cell_to_spot = []
    estimated_cell_number = R.obs['estimated_cell_number']
    for i in range(transport_matrix.shape[1]):
        k = int(estimated_cell_number[i])
        top_k_cells = np.argsort(-transport_matrix[:, i])[:k]
        cell_to_spot.append(top_k_cells)

    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Creating new AnnData object...")
    if issparse(S.X):
        S.X = S.X.toarray()

    # original
    original_spot = list(R.obs_names)
    original_cell = list(S.obs_names)
    original_celltype = list(S.obs[celltype_key])
    original_x = list(R.obsm['spatial'][:, 0])
    original_y = list(R.obsm['spatial'][:, 1])
    original_expr = S.X

    # new
    spot_list = []
    cell_list = []
    celltype_list = []
    x_list = []
    y_list = []
    expr_list = []

    for i, indices in enumerate(cell_to_spot):
        spot_list.extend([original_spot[i]] * len(indices))
        cell_list.extend(original_cell[idx] for idx in indices)
        celltype_list.extend(original_celltype[idx] for idx in indices)

        x_list.extend([original_x[i]] * len(indices))
        y_list.extend([original_y[i]] * len(indices))

        expr_list.extend(original_expr[indices])

    new_idx_list = ['CID' + str(i + 1) for i in range(len(cell_list))]
    new_obs = pd.DataFrame(
        {
            'CellID_new': new_idx_list,
            'CellID_original': cell_list,
            'CellType': celltype_list,
            'SpotID': spot_list,
            'X': x_list,
            'Y': y_list,
        }
    )
    new_obs.index = new_idx_list
    new_expr = np.array(expr_list)
    new_expr = scipy.sparse.csr_matrix(new_expr)

    coord = np.array(new_obs[['X', 'Y']])
    coord_jitter = jitter_coord(coord)
    new_obs['X_jitter'] = coord_jitter[:, 0]
    new_obs['Y_jitter'] = coord_jitter[:, 1]

    # new AnnData object
    adata = sc.AnnData(new_expr)
    adata.obs = new_obs
    adata.obs_names = new_idx_list
    adata.var_names = S.var_names
    adata.obsm['spatial_raw'] = coord
    adata.obsm['spatial'] = coord_jitter

    return adata


def get_celltype_idx(S: AnnData, R: AnnData, S_celltype_key: str, R_celltype_key: str):
    """
    Get the indices of cells, divided by cell tyoe.

    Args:
        S (AnnData): Single cell AnnData object.
        R (AnnData): Spatial reference AnnData object.
        S_celltype_key (str): The column name of `cell type` in `S.obs`.
        R_celltype_key (str): The column name of `cell type` in `R.obs`.

    Returns:
        2D-list: The cell indices for single cell AnnData object (S_idx_list), and for spatial reference AnnData object (R_idx_list).
        e.g., if cells in single cell AnnData objec `S` is [celltype A, celltype B, celltype C, celltype A], the `S_idx_list` should be [[0, 3], [1], [2]].
    """
    celltype_order = sorted(list(set(S.obs[S_celltype_key])))
    S_obs = S.obs.copy()
    R_obs = R.obs.copy()
    S_obs = S_obs.reset_index(drop=True)
    R_obs = R_obs.reset_index(drop=True)

    S_idx_list = []
    R_idx_list = []
    for celltype in celltype_order:
        S_idx = list(S_obs[S_obs[S_celltype_key] == celltype].index)
        R_idx = list(R_obs[R_obs[R_celltype_key] == celltype].index)
        S_idx_list.append(S_idx)
        R_idx_list.append(R_idx)

    return S_idx_list, R_idx_list


def find_celltype(idx_list: list, MAX_NUM: int = 50000):
    """
    Get cell types with cell counts exceeding the threshold (MAX_NUM) and return their indices.
    To prevent out of memory in the OT step caused by excessive cell counts in any cell type during high-resolution spatial mapping.

    Args:
        idx_list (list): The list of cell counts of all cell types.
        MAX_NUM (int, optional): The cell count threshold. Defaults to 50000.

    Returns:
        List: The indices for cell type with excessive cell counts.
    """
    return [i for i, num in enumerate(idx_list) if num > MAX_NUM]


def split_list(idx_list: list, n: int, shuffle: bool = True, seed: int = 0):
    """
    Split the list of cell indices for cell type with excessive cell counts.
    To prevent out of memory in the OT step caused by excessive cell counts in any cell type during high-resolution spatial mapping.

    Args:
        idx_list (list): The list of cell indices for cell type with excessive cell counts.
        n (int): Number of parts to split.
        shuffle (bool, optional): Whether to shuffle the indices before splitting. Defaults to True.
        seed (int, optional): Random seed. Defaults to 0.

    Returns:
        2D-list: The splitted cell indices. e.g., if `idx_list` = [1, 2, 3, 4, 5], `n` = 2, `shuffle` = False, the `result` should be [[1, 2, 3], [4, 5]].
    """
    if shuffle:
        set_seed(seed)
        random.shuffle(idx_list)

    size = len(idx_list) // n
    remainder = len(idx_list) % n
    result = [idx_list[i * size + min(i, remainder):(i + 1) * size + min(i + 1, remainder)] for i in range(n)]
    return result


def create_new_list(S_idx_list: list, R_idx_list: list, celltype_idx: list, ratio_list: list):
    """
    Create new cell indices for single cell AnnData object and spatial reference AnnData object after splitting.
    For `R_idx_list`, split the cell type (in `celltype_idx`) into N parts (in `ratio_list`).
    For `S_idx_list`, duplicate the corresponding cell type.
    Cell types  that are not in `celltype_idx` remain the same.

    Args:
        S_idx_list (list): The original cell indices for single cell AnnData object. Output of the `get_celltype_idx` function.
        R_idx_list (list): The original cell indices for spatial reference AnnData object. Output of the `get_celltype_idx` function.
        celltype_idx (list): The indices for cell type with excessive cell counts. Output of the `find_celltype` function.
        ratio_list (list): The number of parts to split for each cell type with excessive cell counts.
    Returns:
        2D-list: The new cell indices for single cell AnnData object (S_idx_list_new), and for spatial reference AnnData object (R_idx_list_new).
        e.g., if `R_idx_list = [[1, 2, 3, 4], [5, 6], [7, 8]], `S_idx_list` = [[1, 5, 7], [2, 3, 6, 9], [4]], `celltype_idx` = [0, 2], `ratio_list` = [3, 2], `shuffle` = False,
        the `R_idx_list_new` should be [[1, 2], [3], [4], [5, 6], [7], [8]], `S_idx_list_new` should be [[1, 5, 7], [1, 5, 7], [1, 5, 7], [2, 3, 6, 9], [4], [4]].
    """
    S_idx_list_new = []
    R_idx_list_new = []
    ratio_dict = dict(zip(celltype_idx, ratio_list))

    for idx in range(len(R_idx_list)):
        if idx in ratio_dict.keys():
            split_part = split_list(idx_list=R_idx_list[idx].copy(), n=ratio_dict[idx])
            R_idx_list_new.extend(split_part)
            S_idx_list_new.extend([S_idx_list[idx].copy()] * ratio_dict[idx])
        else:
            S_idx_list_new.append(S_idx_list[idx].copy())
            R_idx_list_new.append(R_idx_list[idx].copy())

    return S_idx_list_new, R_idx_list_new
