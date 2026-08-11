import pandas as pd
import numpy as np
import scanpy as sc
from ._utils import *
from concurrent.futures import ThreadPoolExecutor
from scanpy import AnnData
from datetime import datetime
from typing import Optional, List
from sklearn.neighbors import NearestNeighbors


class LowResMapper:
    def __init__(
            self,
            S: AnnData,
            R: AnnData,
            deconv_res: pd.DataFrame,
            celltype_key: str = 'celltype',
            estimate_cell_number_list: Optional[list] = None,
            mean_cell_numbers: int = 5,
            normalize: bool = True,
            numItermax: int = 1e6,
            seed: int = 0,
    ):
        """
        Instantiate the scPositioner mapper (low resolution spatial reference).

        Args:
            S (AnnData): Single cell AnnData object.
            R (AnnData): Spatial reference AnnData object.
            deconv_res (pd.DataFrame): Deconvolution results of the spatial reference, shape = (spot, celltype).
            celltype_key (str, optional): The column name of `cell type` in `S.obs`. Defaults to 'celltype'.
            estimate_cell_number_list (Optional[list], optional): The number of cells per spot. Defaults to None.
            mean_cell_numbers (int, optional): The average number of cells across all spots. Defaults to 5.
            normalize (bool, optional): Whether to normalize data. Defaults to True.
            numItermax (int, optional): The maximum number of iterations before stopping the optimization algorithm if it has not converged. Defaults to 1e6.
            seed (int, optional): Random seed. Defaults to 0.

        Returns:
            AnnData: The mapped AnnData object.
        """
        self.S = S
        self.R = R
        self.deconv_res = deconv_res
        self.celltype_key = celltype_key
        self.estimate_cell_number_list = estimate_cell_number_list
        self.mean_cell_numbers = mean_cell_numbers
        self.normalize = normalize
        self.numItermax = numItermax
        self.seed = seed

    def _check(self):
        print(
            f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Checking single cell AnnData: "
            f"{self.S.shape[0]} cells * {self.S.shape[1]} features, "
            f"{len(set(self.S.obs[self.celltype_key]))} cell types."
        )

        # check estimate_cell_number_list
        if self.estimate_cell_number_list is not None:
            if min(self.estimate_cell_number_list) < 0:
                raise ValueError('Estimated cell numbers must be non-negative!')
            if len(self.estimate_cell_number_list) != self.R.shape[0]:
                raise ValueError('The length of estimate_cell_number_list must be equal to the number of spots in the input spatial AnnData!')
            self.R.obs['estimated_cell_number'] = self.estimate_cell_number_list

        print(
            f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Checking spatial reference AnnData: "
            f"{self.R.shape[0]} spots * {self.R.shape[1]} features; "
            f"'estimate_cell_number_list' is provided: {self.estimate_cell_number_list is not None}."
        )

        # check deconv_res, if it doesn't exsist, automatically run cell2location first
        if self.deconv_res is not None:
            self.deconv_res = check_deconvolution_results(S=self.S, R=self.R, deconv_res=self.deconv_res, celltype_key=self.celltype_key)
        else:
            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Deconvolution result not provided, automatically run cell2location...")
            self.deconv_res = run_cell2location(adata_sc=self.S, adata_st=self.R, cell_num_per_spot=self.mean_cell_numbers, celltype_key=self.celltype_key)
            self.deconv_res = check_deconvolution_results(S=self.S, R=self.R, deconv_res=self.deconv_res, celltype_key=self.celltype_key)
            
        self.deconv_res.columns = ['deconv_' + str(i) for i in self.deconv_res.columns]
        R_obs = self.R.obs.copy()
        R_obs_new = pd.concat([R_obs, self.deconv_res], axis=1)
        self.R.obs = R_obs_new

    def run(self):
        # check
        self._check()

        # filtering: features and cells
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Filtering features and cells...")
        self.S = filtering(self.S)
        self.R = filtering(self.R)

        # copy
        R_raw = self.R.copy()
        S_raw = self.S.copy()

        # estimate cell number
        if 'estimated_cell_number' not in R_raw.obs.keys():
            self.estimate_cell_number_list = estimate_cell_number(
                adata=R_raw, mean_cell_numbers=self.mean_cell_numbers, normalize=self.normalize
            )
            R_raw.obs['estimated_cell_number'] = self.estimate_cell_number_list

        # adjust cell abundance
        S_raw = adjust_abundance(S=S_raw, R=R_raw, celltype_key=self.celltype_key, seed=self.seed)

        # common features
        common_feat = S_raw.var_names.intersection(R_raw.var_names)
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {len(common_feat)} common features are used for mapping...")
        S_use = S_raw[:, common_feat].copy()
        R_use = R_raw[:, common_feat].copy()

        # normalize
        if self.normalize:
            S_use = normalize_log(S_use)
            R_use = normalize_log(R_use)

        # optimal transport
        transport_matrix = run_OT(S=S_use, R=R_use, numItermax=self.numItermax)

        # post-processing
        adata = post_process(S=S_raw, R=R_raw, transport_matrix=transport_matrix, celltype_key=self.celltype_key)
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Done!")

        return adata

    
class HighResMapper:
    def __init__(
            self,
            S: AnnData,
            R: AnnData,
            S_celltype_key: str = 'celltype',
            R_celltype_key: str = 'celltype',
            normalize: bool = True,
            max_number: int = 20000,
            numItermax: int = 1e6,
            max_workers: int = 10,
            seed: int = 0,
    ):
        """
        Instantiate the scPositioner mapper (high resolution spatial reference).

        Args:
            S (AnnData): Single cell AnnData object.
            R (AnnData): Spatial reference AnnData object.
            S_celltype_key (str, optional): The column name of `cell type` in `S.obs`. Defaults to 'celltype'.
            R_celltype_key (str, optional): The column name of `cell type` in `R.obs`. Defaults to 'celltype'.
            normalize (bool, optional): Whether to normalize data. Defaults to True.
            max_number (int, optional): The maximum number of cell number in each chunk. Defaults to 2e5.
            numItermax (int, optional): The maximum number of iterations before stopping the optimization algorithm if it has not converged. Defaults to 1e6.
            max_workers (int, optional): Number of threads. Defaults to 10.
            seed (int, optional): Random seed. Defaults to 0.

        Returns:
            AnnData: The mapped AnnData object.
        """
        self.S = S
        self.R = R
        self.S_celltype_key = S_celltype_key
        self.R_celltype_key = R_celltype_key
        self.normalize = normalize
        self.numItermax = numItermax
        self.max_workers = max_workers
        self.seed = seed
        self.max_number = max_number

    def _check(self):
        print(
            f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Checking single cell AnnData: "
            f"{self.S.shape[0]} cells * {self.S.shape[1]} features, "
            f"{len(set(self.S.obs[self.S_celltype_key]))} cell types."
        )

        print(
            f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Checking spatial reference AnnData: "
            f"{self.R.shape[0]} cells * {self.R.shape[1]} features; "
            f"{len(set(self.R.obs[self.R_celltype_key]))} cell types."
        )

        common_celltypes = set(self.S.obs[self.S_celltype_key]).intersection(set(self.R.obs[self.R_celltype_key]))
        if len(common_celltypes) == 0:
            raise ValueError('No common cell types found between single cell and spatial reference AnnData!')
        print(
            f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {len(common_celltypes)} common cell types found between single cell and spatial reference AnnData."
        )

        self.max_workers = min(self.max_workers, len(common_celltypes))
        print(
            f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {self.max_workers} cores are used."
        )

    def _process_part(self, i, S, R, S_idx_list, R_idx_list):
        # print(f"------- Part {i} -------")
        S_part_raw = S[S_idx_list[i]].copy()
        R_part_raw = R[R_idx_list[i]].copy()

        # adjust cell abundance
        S_part_raw = adjust_abundance(S=S_part_raw, R=R_part_raw, celltype_key=self.S_celltype_key, seed=self.seed)

        # common features
        common_feat = S_part_raw.var_names.intersection(R_part_raw.var_names)
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {len(common_feat)} common features are used for mapping...")
        S_part_use = S_part_raw[:, common_feat].copy()
        R_part_use = R_part_raw[:, common_feat].copy()

        # normalize
        if self.normalize:
            S_part_use = normalize_log(S_part_use)
            R_part_use = normalize_log(R_part_use)

        # optimal transport
        transport_matrix = run_OT(S=S_part_use, R=R_part_use, numItermax=self.numItermax)

        # post-processing
        adata = post_process(S=S_part_raw, R=R_part_raw, transport_matrix=transport_matrix, celltype_key=self.S_celltype_key)

        return adata

    def run(self):
        # check
        self._check()

        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Filtering features and cells...")
        self.S = filtering(self.S)
        self.R = filtering(self.R)

        R_raw = self.R.copy()
        S_raw = self.S.copy()

        # estimate cell number
        R_raw.obs['estimated_cell_number'] = [1] * R_raw.shape[0]

        # add deconv_res
        celltype_order = sorted(list(set(S_raw.obs[self.S_celltype_key])))
        deconv_res = pd.get_dummies(R_raw.obs[self.R_celltype_key])
        deconv_res = deconv_res[celltype_order]
        deconv_res = deconv_res.loc[R_raw.obs_names, :]
        deconv_res.columns = ['deconv_' + str(i) for i in deconv_res.columns]
        R_obs = R_raw.obs.copy()
        R_obs_new = pd.concat([R_obs, deconv_res], axis=1)
        R_raw.obs = R_obs_new

        # cell type idx list
        S_idx_list, R_idx_list = get_celltype_idx(S=S_raw, R=R_raw,
                                                 S_celltype_key=self.S_celltype_key,
                                                 R_celltype_key=self.R_celltype_key)

        if self.max_number is not None:
            celltype_num_list = [len(R_idx_list[i]) for i in range(len(R_idx_list))]
            celltype_over = find_celltype(idx_list=celltype_num_list, MAX_NUM=self.max_number)
            if len(celltype_over) > 0:
                ratio = [int(np.ceil(celltype_num_list[i] / self.max_number)) for i in celltype_over]
                S_idx_list, R_idx_list = create_new_list(S_idx_list=S_idx_list,
                                                         R_idx_list=R_idx_list,
                                                         celltype_idx=celltype_over,
                                                         ratio_list=ratio)

        # run
        adata_list = []
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            results = [
                executor.submit(self._process_part, i, S_raw, R_raw, S_idx_list, R_idx_list)
                for i in range(len(S_idx_list))
            ]
            for res in results:
                adata_list.append(res.result())

        adata = sc.concat(adata_list).copy()
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Done!")

        return adata
    

class StabilityAnalysis:
    def __init__(
            self,
            S: sc.AnnData,
            R: sc.AnnData,
            deconv_res: pd.DataFrame = None,
            celltype_key: str = 'celltype',
            estimate_cell_number_list: Optional[List] = None,
            mean_cell_numbers: int = 5,
            normalize: bool = True,
            numItermax: int = 1e6,
            dropout_rates: List[float] = [0, 0.1, 0.2, 0.3, 0.5, 0.8, 0.9, 0.92, 0.95, 0.97, 0.98, 0.99],
            seeds: List[int] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    ):
        self.S = S
        self.R = R
        self.deconv_res = deconv_res
        self.celltype_key = celltype_key
        self.estimate_cell_number_list = estimate_cell_number_list
        self.mean_cell_numbers = mean_cell_numbers
        self.normalize = normalize
        self.numItermax = numItermax
        self.dropout_rates = dropout_rates
        self.seeds = seeds
        self.mapping_df = None

    def _check(self):
        print(
            f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Checking single cell AnnData: "
            f"{self.S.shape[0]} cells * {self.S.shape[1]} features, "
            f"{len(set(self.S.obs[self.celltype_key]))} cell types."
        )
        print(
            f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Checking spatial reference AnnData: "
            f"{self.R.shape[0]} spots * {self.R.shape[1]} features."
        )

        if self.estimate_cell_number_list is not None:
            if len(self.estimate_cell_number_list) != self.R.shape[0]:
                raise ValueError("estimate_cell_number_list length mismatch")

        # check deconv_res, if it doesn't exsist, automatically run cell2location first
        if self.deconv_res is not None:
            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Validating provided deconvolution results...")
            self.deconv_res = check_deconvolution_results(S=self.S, R=self.R, deconv_res=self.deconv_res, celltype_key=self.celltype_key)
        else:
            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Deconvolution result not provided, automatically run cell2location...")
            self.deconv_res = run_cell2location(adata_sc=self.S, adata_st=self.R, cell_num_per_spot=self.mean_cell_numbers, celltype_key=self.celltype_key)
            self.deconv_res = check_deconvolution_results(S=self.S, R=self.R, deconv_res=self.deconv_res, celltype_key=self.celltype_key)

        self.deconv_res.columns = ['deconv_' + str(c) for c in self.deconv_res.columns]
        self.R.obs = pd.concat([self.R.obs, self.deconv_res], axis=1)
        print(
            f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Deconvolution results merged into R.obs "
            f"({self.deconv_res.shape[1]} columns added)."
        )

    def _calculate_stability_scores(self):
        total_runs = len(self.dropout_rates) * len(self.seeds)
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Calculating per-dropout-rate stability scores...")

        for rate in self.dropout_rates:
            cols = [f'dropout_rate_{rate}_seed_{s}' for s in self.seeds]

            def mode_func(row):
                vals = row[cols].dropna().values
                if len(vals) == 0:
                    return np.nan, 0.0
                ser = pd.Series(vals)
                mode_val = ser.mode().iloc[0]
                freq = (ser == mode_val).sum() / len(self.seeds)
                return mode_val, freq

            res = self.mapping_df.apply(mode_func, axis=1)
            self.mapping_df[f'dropout_rate_{rate}_assigned_spot'] = res.apply(lambda x: x[0])
            self.mapping_df[f'dropout_rate_{rate}_stability'] = res.apply(lambda x: x[1])

        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Calculating overall stability scores...")

        all_cols = [
            f'dropout_rate_{r}_seed_{s}'
            for r in self.dropout_rates for s in self.seeds
        ]

        def overall_mode(row):
            vals = row[all_cols].dropna().values
            if len(vals) == 0:
                return np.nan, 0.0
            ser = pd.Series(vals)
            mode_val = ser.mode().iloc[0]
            freq = (ser == mode_val).sum() / total_runs
            return mode_val, freq

        res = self.mapping_df.apply(overall_mode, axis=1)
        self.mapping_df['final_assigned_spot'] = res.apply(lambda x: x[0])
        self.mapping_df['overall_stability'] = res.apply(lambda x: x[1])
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Stability scores computed for {len(self.mapping_df)} cells.")

    def _assign_spatial_coordinates(self):
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Assigning spatial coordinates with jittering...")
        spot_coords = self.R.obsm['spatial']
        spot_dict = dict(zip(self.R.obs_names, spot_coords))

        n_cells = len(self.mapping_df)
        jittered = np.full((n_cells, 2), np.nan, dtype=np.float32)

        for spot in self.mapping_df['final_assigned_spot'].dropna().unique():
            idx = np.where(self.mapping_df['final_assigned_spot'] == spot)[0]
            if len(idx) == 0:
                continue

            base = spot_dict.get(spot, None)
            if base is None or np.isnan(base).any():
                continue

            seed = int(hash(spot) % (2**31))
            inp = np.tile(base, (len(idx), 1))
            jit = jitter_coord(inp, seed=seed)

            if jit.ndim == 1:
                jit = jit.reshape(-1, 2)

            jittered[idx] = jit

        assert jittered.shape == (n_cells, 2), "Coordinate shape mismatch!"
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Spatial coordinates assigned for {n_cells} cells.")
        return jittered

    def run(self):
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] ===== Starting stability analysis =====")
        self._check()

        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Filtering features and cells...")
        self.S = filtering(self.S)
        self.R = filtering(self.R)

        R_raw = self.R.copy()
        S_raw = self.S.copy()

        if 'estimated_cell_number' not in R_raw.obs.columns:
            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Estimating cell numbers per spot...")
            R_raw.obs['estimated_cell_number'] = estimate_cell_number(
                adata=R_raw,
                mean_cell_numbers=self.mean_cell_numbers,
                normalize=self.normalize
            )

        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Adjusting cell abundance...")
        S_raw = adjust_abundance(
            S=S_raw,
            R=R_raw,
            celltype_key=self.celltype_key,
            seed=self.seeds[0]
        )

        self.mapping_df = pd.DataFrame(index=S_raw.obs_names)

        total_runs = len(self.dropout_rates) * len(self.seeds)
        print(
            f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Running {total_runs} mappings "
            f"({len(self.dropout_rates)} dropout rates × {len(self.seeds)} seeds)..."
        )

        run_count = 0
        for ri, rate in enumerate(self.dropout_rates):
            for si, seed in enumerate(self.seeds):
                run_count += 1
                print(
                    f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] "
                    f"[{run_count}/{total_runs}] dropout_rate={rate}, seed={seed}",
                    flush=True
                )
                S_do = simulate_gene_dropout(S_raw, rate, seed)
                common = S_do.var_names.intersection(R_raw.var_names)
                if len(common) == 0:
                    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]   -> no common features, skip.")
                    continue

                S_use = S_do[:, common].copy()
                R_use = R_raw[:, common].copy()

                if self.normalize:
                    S_use = normalize_log(S_use)
                    R_use = normalize_log(R_use)

                tm = run_OT(S=S_use, R=R_use, numItermax=self.numItermax)
                adata_tmp = post_process(
                    S=S_raw, R=R_raw,
                    transport_matrix=tm,
                    celltype_key=self.celltype_key
                )

                if 'CellID_original' in adata_tmp.obs and 'SpotID' in adata_tmp.obs:
                    mp = adata_tmp.obs.set_index('CellID_original')['SpotID']
                    self.mapping_df[f'dropout_rate_{rate}_seed_{seed}'] = \
                        mp.reindex(self.mapping_df.index)

        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] All mappings done, calculating stability scores...")
        self._calculate_stability_scores()

        final_adata = S_raw.copy()

        if not self.mapping_df.index.equals(final_adata.obs.index):
            self.mapping_df = self.mapping_df.reindex(final_adata.obs.index)

        final_adata.obs = pd.concat([final_adata.obs, self.mapping_df], axis=1)

        jittered_coords = self._assign_spatial_coordinates()
        jittered_coords = jitter_coord(jittered_coords)
        final_adata.obsm['spatial'] = jittered_coords

        final_adata.uns['stability_analysis'] = {
            'dropout_rates': self.dropout_rates,
            'seeds': self.seeds,
            'total_runs': len(self.dropout_rates) * len(self.seeds),
            'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        }

        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] ===== Stability analysis finished =====")
        return final_adata
