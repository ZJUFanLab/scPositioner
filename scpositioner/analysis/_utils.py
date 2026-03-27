import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
import scipy.stats as stats
from scanpy import AnnData
from datetime import datetime
from sklearn.neighbors import KDTree


def knn(data: np.ndarray, query: np.ndarray, k: int = 5):
    """
    Find K-Nearest Neighbor.

    Args:
        data (np.ndarray): Data.
        query (np.ndarray): Query data.
        k (int, optional): Number of nearest neighbors. Defaults to 5.

    Returns:
        np.ndarray: dist: distance, ind: indices
    """
    tree = KDTree(data)
    dist, ind = tree.query(query, k)

    return dist, ind


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


def assign_group(
        adata: AnnData,
        select_celltype: list,
        background_celltype: list,
        celltype_key: str,
        spatial_key: str = 'spatial_raw',
        cutoff_type: str = 'median',
        n_neighbor: int = 5):
    """
    Assign cell type A (select_celltype) to "Close" and "Far" groups based on the distance from another cell type B (background_celltype).
    The distance for each cell of cell type A is defined as the average distance to the n nearest cells of cell type B.

    Args:
        adata (AnnData): AnnData object.
        select_celltype (list): The selected cell type to assign. e.g., ['cell type A', ...]
        background_celltype (list): The background cell type. e.g., ['cell type B', ...]
        celltype_key (str): The column name of `cell type` in `adata.obs`.
        spatial_key (str, optional): The key of spatial coordinates in `adata.obsm`. For the output of scPositioner,
                                     'spatial_raw' is the original coordinates of spatial reference, and 'spatial' is the jittered coordinates.
                                     Defaults to 'spatial_raw'.
        cutoff_type (str, optional): Method for setting cutoff values. if `cutoff_type` = 'median', the cutoff value is set to the median of the average distance of all cells,
                                     else if `cutoff_type` = 'mean', the cutoff value is set to the mean of the average distance of all cells. Defaults to 'median'.
        n_neighbor (int, optional): Number of nearest neighbors. Defaults to 5.

    Returns:
        adata_select_new: AnnData object of `select_celltype`. The assigned group is stored in `adata_select_new.obs['Assigned_group']`.
    """

    print(
        f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Assigning {select_celltype} "
        f"based on the average distance to the {n_neighbor} nearest {background_celltype}..."
    )

    coord1 = np.array(adata[adata.obs[celltype_key].isin(select_celltype)].obsm[spatial_key]).copy()
    coord2 = np.array(adata[adata.obs[celltype_key].isin(background_celltype)].obsm[spatial_key]).copy()

    dist, _ = knn(data=coord2, query=coord1, k=n_neighbor)
    dist_mean = dist.mean(axis=1)

    if cutoff_type == 'median':
        cutoff = np.median(dist_mean)
    elif cutoff_type == 'mean':
        cutoff = np.mean(dist_mean)
    else:
        raise ValueError("cutoff_type must be 'mean' or 'median'!")

    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] The cutoff value is the {cutoff_type} of the average distance of all cells: {cutoff}")

    idx1 = np.argwhere(dist_mean > cutoff).flatten().tolist()
    idx2 = np.argwhere(dist_mean <= cutoff).flatten().tolist()

    adata_select = adata[adata.obs[celltype_key].isin(select_celltype)].copy()
    adata_select.obs['Assigned_group'] = 'Others'
    adata_select.obs.loc[adata_select.obs.index[idx1], 'Assigned_group'] = 'Far'
    adata_select.obs.loc[adata_select.obs.index[idx2], 'Assigned_group'] = 'Close'

    adata1 = adata_select[adata_select.obs['Assigned_group'].isin(['Close'])].copy()
    adata2 = adata_select[adata_select.obs['Assigned_group'].isin(['Far'])].copy()
    adata_select_new = sc.concat([adata1, adata2])
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Done!")

    return adata_select_new


def run_GSEA(
        adata: AnnData,
        geneset: dict,
        group_key: str = 'Assigned_group',
        normalize: bool = True,
        permutation_num: int = 1000,
        permutation_type: str = 'phenotype',
        method: str = 's2n',  # signal_to_noise,
        threads: int = 32,
        min_size: int = 1,
        max_size: int = 500,
        seed: int = 123):
    """
    Run Gene Set Enrichment Analysis.
    See https://gseapy.readthedocs.io/en/latest/_modules/gseapy.html#gsea for details.

    Args:
        adata (AnnData): AnnData object.
        geneset (dict): Geneset. e.g., {'geneset1': ['gene1', 'gene2'], 'geneset2': ['gene3', 'gene4'], ...}
        group_key (str): The column name of `assigned group` in `adata.obs`. Defaults to 'Assigned_group'.
        normalize (bool, optional): Whether to normalize data. Defaults to True.
        permutation_num (int, optional): Number of permutations. Defaults to 1000.
        permutation_type (str, optional): Type of permutation reshuffling. Defaults to 'phenotype'.
        method (str, optional): The method used to calculate a correlation or ranking. Defaults to 's2n'.
        threads (int, optional): Number of threads. Defaults to 32.
        min_size (int, optional): Minimum allowed number of genes from gene set also the data set. Defaults to 1.
        max_size (int, optional): Maximum allowed number of genes from gene set also the data set. Defaults to 500.
        seed (int, optional): Random. Defaults to 123.

    Returns:
        adata: AnnData object. The GSEA result is stored in `adata.uns['GSEA_results']`.
    """
    adata_use = adata.copy()
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Running GSEA ...")
    if normalize:
        adata_use = normalize_log(adata_use)

    # run GSEA
    expr = adata_use.to_df().T
    anno = adata_use.obs[group_key]
    select_celltype = anno[0]
    background_celltype = anno[anno != select_celltype][0]

    GSEA_res = gp.gsea(
        data=expr,
        gene_sets=geneset,
        cls=anno,
        permutation_num=permutation_num,
        permutation_type=permutation_type,
        method=method,
        threads=threads,
        min_size=min_size,
        max_size=max_size,
        seed=seed)

    res = {}
    res['GSEA_obj'] = GSEA_res
    res['compare'] = f"{str(select_celltype)} v.s. {str(background_celltype)}"

    adata.uns['GSEA_results'] = res

    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Done!")

    return adata


def calculate_score(
        adata: AnnData,
        signature: list,
        select_group: str = 'Close',
        background_group: str = 'Far',
        group_key: str = 'Assigned_group',
        normalize: bool = True,
        calculate_p: bool = True):
    """
    Calculate gene signature score.

    Args:
        adata (AnnData): AnnData object.
        signature (list): gene signature list.
        select_group (str, optional):  The selected group. Defaults to 'Close'.
        background_group (str, optional): The background group. Defaults to 'Far'.
        group_key (str, optional): The column name of `assigned group` in `adata.obs`. Defaults to 'Assigned_group'.
        normalize (bool, optional): Whether to normalize data. Defaults to True.
        calculate_p (bool, optional): Whether to calculate p-values. Defaults to True.

    Returns:
        adata: AnnData object. The gene signature score result is stored in `adata.uns['score_results']`.
    """
    adata_use = adata.copy()
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]Calculating gene signature score ...")
    if normalize:
        adata_use = normalize_log(adata_use)

    signature = list(set(signature) & set(adata_use.var_names))
    sc.tl.score_genes(adata_use, signature, )
    score1 = list(adata_use[adata_use.obs[group_key] == select_group].obs['score'])
    score2 = list(adata_use[adata_use.obs[group_key] == background_group].obs['score'])

    p_value = None
    method = None
    if calculate_p:
        _, p_value = stats.mannwhitneyu(score1, score2, alternative='two-sided')
        method = 'two sided Mann-Whitney U test'

    group1 = [select_group] * len(score1)
    group2 = [background_group] * len(score2)

    score_all = score1 + score2
    group_all = group1 + group2

    df = pd.DataFrame({'Group': group_all, 'Score': score_all})

    res = {}
    res['score_df'] = df
    res['calculate_p'] = calculate_p
    res['method'] = method
    res['p_value'] = p_value
    res['signature'] = signature

    adata.uns['score_results'] = res

    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Done!")

    return adata


def calculate_dist(
        adata: AnnData,
        select_celltype: list,
        background_celltype: list,
        celltype_key: str,
        spatial_key: str = 'spatial',
        dist_type: str = 'mean',
        reverse: bool = False,
        calculate_linear_corr: bool = False):
    """
    Calculate Euclidean distance from cell type A (select_celltype) to cell type B (background_celltype) based on the spatial coordinates.

    Args:
        adata (AnnData): AnnData object.
        select_celltype (list): The selected cell type to assign. e.g., ['cell type A', ...]
        background_celltype (list): The background cell type. e.g., ['cell type B', ...]
        celltype_key (str): The column name of `cell type` in `adata.obs`.
        spatial_key (str, optional): The key of spatial coordinates in `adata.obsm`. For the output of scPositioner,
                                     'spatial_raw' is the original coordinates of spatial reference, and 'spatial' is the jittered coordinates.
                                     Defaults to 'spatial'.
        dist_type (str, optional): Type of distance. if `dist_type` = 'mean', the distance for each cell in cell type A (select_celltype) is defined as the mean
        distance to all cells in cell type B (background_celltype), else if `dist_type` = 'median', the distance is defined as the median distance to all cells in cell type B,
        and else if `dist_type` = 'minimum', the distance is defined as the minimum distance to all cells in cell type B. Defaults to 'mean'.
        reverse (bool, optional): Whether to reverse the distribution order. if `reverse` = False, the distribution order is [1, 2, 3, ..., len(select_celltype)],
                                  else if `reverse` = True, the distribution order is [len(select_celltype), len(select_celltype) - 1, ..., 1]. Defaults to False.
        calculate_linear_corr (bool, optional): Whether to calculate the linear correlation between mean Euclidean distance and distribution order. Defaults to False.

    Returns:
       adata: AnnData object. The distance result is stored in `adata.uns['distance_results']`.
    """

    if dist_type not in ['mean', 'median', 'minimum']:
        raise ValueError("dist_type must be 'mean', 'median', or 'minimum'!")

    print(
        f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Calculating the {dist_type} distance from each cell in {select_celltype} to {background_celltype}..."
    )

    dist_mean = []
    dist_cell = []
    coord2 = np.array(adata[adata.obs[celltype_key].isin(background_celltype)].obsm[spatial_key]).copy()
    # for i in tqdm(select_celltype):
    for i in select_celltype:
        coord1 = np.array(adata[adata.obs[celltype_key].isin([i])].obsm[spatial_key]).copy()
        if coord1.shape[0] > 0:
            dist, _ = knn(data=coord2, query=coord1, k=coord2.shape[0])
            if dist_type == 'mean':
                dist1 = np.mean(dist, axis=1)
            elif dist_type == 'median':
                dist1 = np.median(dist, axis=1)
            elif dist_type == 'minimum':
                dist1 = np.min(dist, axis=1)
            dist_mean1 = dist1.mean()
        else:
            dist1 = np.nan
            dist_mean1 = np.nan
        dist_mean.append(dist_mean1)
        dist_cell.append(dist1)

    order_list = [i + 1 for i in range(len(select_celltype))]
    if reverse:
        order_list = order_list[::-1]

    df = pd.DataFrame({'Selected cell type': select_celltype, 'Mean distance to background': dist_mean, 'Distribution order': order_list, })

    p_value = None
    r = None
    if calculate_linear_corr:
        dist_list1 = list(df['Mean distance to background'])
        order_list1 = list(df['Distribution order'])
        idx = np.where(~np.isnan(np.array(dist_list1)))[0]
        dist_list1 = [dist_list1[i] for i in idx]
        order_list1 = [order_list1[i] for i in idx]

        r, p_value = stats.pearsonr(order_list1, dist_list1)

    res = {}
    res['mean_distance_df'] = df
    res['distance_cell_list'] = dist_cell
    res['select_celltype'] = select_celltype
    res['order_list'] = order_list
    res['background_celltype'] = background_celltype
    res['dist_type'] = dist_type
    res['r'] = r
    res['p_value'] = p_value

    adata.uns['distance_results'] = res

    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Done!")

    return adata
