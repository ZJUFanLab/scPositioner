import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
from anndata import AnnData
from typing import Optional, Union, Sequence
from matplotlib.gridspec import GridSpec
import seaborn as sns
import scanpy as sc
import warnings
import operator
import matplotlib.transforms as transforms
from matplotlib.colors import ListedColormap
from . import palettes


# from cell2location: https://github.com/BayraktarLab/cell2location/blob/master/cell2location/plt/plot_spatial.py
def _create_colormap(R, G, B, white_spacing):
    """
    Creat colormap.
    """
    spacing = int(white_spacing * 2.55)

    N = 255
    M = 3

    alphas = np.concatenate([[0] * spacing * M, np.linspace(0, 1.0, (N - spacing) * M)])

    vals = np.ones((N * M, 4))
    for i, color in enumerate([R, G, B]):
        vals[:, i] = color / 255
    vals[:, 3] = alphas

    return ListedColormap(vals)


def _get_rgb_function(cmap, min_value, max_value):
    """
    Generate a function to map continous values to RGB values using colormap between min_value & max_value.
    """
    if min_value >= max_value:
        raise ValueError("Max_value should be greater or than min_value.")

    if min_value == max_value:
        warnings.warn(
            "Max_color is equal to min_color. It might be because of the data or bad parameter choice. "
            "If you are using plot_contours function try increasing max_color_quantile parameter and"
            "removing cell types with all zero values."
        )

        def func_equal(x):
            factor = 0 if max_value == 0 else 0.5
            return cmap(np.ones_like(x) * factor)

        return func_equal

    def func(x):
        return cmap((np.clip(x, min_value, max_value) - min_value) / (max_value - min_value))

    return func


def _hex_to_rgb(hex_color):
    """
    Convert a hexadecimal color code to RGB.

    Args:
        hex_color (str): A hexadecimal color code, e.g., "#FF0000" for red.

    Returns:
        tuple: A tuple of (R, G, B) values, each in the range 0-255.
    """
    # Remove the leading '#' if present
    hex_color = hex_color.lstrip('#')

    # Convert the hexadecimal code to integer
    rgb_int = int(hex_color, 16)

    # Extract the red, green, and blue components
    r = (rgb_int >> 16) & 0xFF
    g = (rgb_int >> 8) & 0xFF
    b = rgb_int & 0xFF

    # Normalize
    r = r / 255
    g = g / 255
    b = b / 255

    return r, g, b


def _rgb_to_ryb(rgb):
    """
    Converts colours from RGB colorspace to RYB

    Args:
        rgb (np.ndarray): Nx3

    Returns:
        np.ndarray: Nx3
    """
    rgb = np.array(rgb)
    if len(rgb.shape) == 1:
        rgb = rgb[np.newaxis, :]

    white = rgb.min(axis=1)
    black = (1 - rgb).min(axis=1)
    rgb = rgb - white[:, np.newaxis]

    yellow = rgb[:, :2].min(axis=1)
    ryb = np.zeros_like(rgb)
    ryb[:, 0] = rgb[:, 0] - yellow
    ryb[:, 1] = (yellow + rgb[:, 1]) / 2
    ryb[:, 2] = (rgb[:, 2] + rgb[:, 1] - yellow) / 2

    mask = ~(ryb == 0).all(axis=1)
    if mask.any():
        norm = ryb[mask].max(axis=1) / rgb[mask].max(axis=1)
        ryb[mask] = ryb[mask] / norm[:, np.newaxis]

    return ryb + black[:, np.newaxis]


def _ryb_to_rgb(ryb):
    """
    Converts colours from RYB colorspace to RGB.

    Args:
        ryb (np.ndarray): Nx3

    Returns:
        np.ndarray: Nx3
    """
    ryb = np.array(ryb)
    if len(ryb.shape) == 1:
        ryb = ryb[np.newaxis, :]

    black = ryb.min(axis=1)
    white = (1 - ryb).min(axis=1)
    ryb = ryb - black[:, np.newaxis]

    green = ryb[:, 1:].min(axis=1)
    rgb = np.zeros_like(ryb)
    rgb[:, 0] = ryb[:, 0] + ryb[:, 1] - green
    rgb[:, 1] = green + ryb[:, 1]
    rgb[:, 2] = (ryb[:, 2] - green) * 2

    mask = ~(ryb == 0).all(axis=1)
    if mask.any():
        norm = rgb[mask].max(axis=1) / ryb[mask].max(axis=1)
        rgb[mask] = rgb[mask] / norm[:, np.newaxis]

    return rgb + white[:, np.newaxis]


def _set_palette(length):
    """
    Set palettes.
    """

    if length <= 10:
        palette = palettes.default_10
    elif length <= 20:
        palette = palettes.default_20
    elif length <= 28:
        palette = palettes.default_28
    elif length <= 57:
        palette = palettes.default_57
    elif length <= len(palettes.default_102):  # 103 colors
        palette = palettes.default_102
    else:
        palette = ['grey' for _ in range(length)]
        print(
            'the obs value has more than 103 categories. Uniform '
            "'grey' color will be used for all categories."
        )

    return palette


def _get_prop(adata: AnnData, ):
    """
    Get cell type aboundance of each spot.

    Args:
        adata (AnnData): AnnData object of scPositioner output.

    Returns:
        adata_sub(AnnData): AnnData of spots; prop(pd.DataFrame): cell type aboundance of each spot.
    """
    obs = adata.obs.copy()
    # cell type proportion of each spot
    # 'SpotID', 'CellType': Mapper output
    prop = obs.groupby('SpotID')['CellType'].value_counts().unstack(fill_value=0)

    # unique
    obs_unique = obs.drop_duplicates(subset=['SpotID'], keep='first')
    idx = list(obs_unique.index)
    adata_sub = adata[idx]
    prop.index = idx
    obs_sub = adata_sub.obs.copy()
    obs_sub_new = pd.concat([obs_sub, prop], axis=1)
    adata_sub.obs = obs_sub_new

    return adata_sub, prop


def _pie_marker(
        ratios: Sequence[float],
        res: int = 50,
        direction: str = "+",
        start: float = 0.0,):
    """
    Create each slice of pie as a separate marker.

    Args:
        ratios (Sequence[float]): List of ratios that add up to 1.
        res (int, optional): Number of points around the circle.. Defaults to 50.
        direction (str, optional): '+' for counter-clockwise, or '-' for clockwise. Defaults to "+".
        start (float, optional): Starting position in radians.. Defaults to 0.0.

    Returns:
        xys(list): list of xy points of each slice; ss(list): list of size of each slice.
    """
    if np.abs(np.sum(ratios) - 1) > 0.01:
        print("Warning: Ratios do not add up to 1.")

    if direction == '+':
        op = operator.add
    elif direction == '-':
        op = operator.sub

    xys = []  # list of xy points of each slice
    ss = []  # list of size of each slice
    start = float(start)
    for ratio in ratios:
        # points on the circle including the origin (0,0) and the slice
        end = op(start, 2 * np.pi * ratio)
        n = round(ratio * res)  # number of points forming the arc
        x = [0] + np.cos(np.linspace(start, end, n)).tolist()
        y = [0] + np.sin(np.linspace(start, end, n)).tolist()
        xy = np.column_stack([x, y])
        xys.append(xy)
        ss.append(np.abs(xy).max())
        start = end

    return xys, ss


def cell_scatterplot(
        adata: AnnData,
        show: Union[str, list],
        spatial_key: str = 'spatial',
        palette: Optional[Union[list, dict]] = None,
        save: bool = False,
        save_dir: str = '',
        kwargs: dict = {}):
    """
    Scatter plot for single cell on the cell coordinates.

    Args:
        adata (AnnData): AnnData object of scPositioner output.
        show (Union[str, list]): The column name of the object to be shown in `adata.obs`.
        spatial_key (str, optional): The key of spatial coordinates in `adata.obsm`.. Defaults to 'spatial'.
        palette (Optional[Union[list, dict]], optional): Palette. Defaults to None.
        save (bool, optional): Whether to save. Defaults to False.
        save_dir (str, optional): Directory to save. Defaults to ''.
        kwargs (dict, optional): Other parameters. Defaults to {}.
    """
    if palette is None:
        if isinstance(show, str):
            length = len(set(adata.obs[show]))
        elif isinstance(show, list):
            length_list = [len(set(adata.obs[i])) for i in show]
            length = max(length_list)
        palette = _set_palette(length=length)

    sc.pl.embedding(adata, basis=spatial_key, color=show, palette=palette, show=False, **kwargs)
    if save:
        plt.savefig(save_dir, format='svg')


def spot_scatterpie(
        adata: AnnData,
        spatial_key: str = 'spatial_raw',
        palette: Optional[list] = None,
        point_size: float = 45,
        save: bool = False,
        save_dir: str = '',
        kwargs: dict = {}):
    """
    Scatterpie plot for each cell type abundance on the spot coordinates.

    Args:
        adata (AnnData): AnnData object of scPositioner output.
        spatial_key (str, optional): The key of spatial coordinates in `adata.obsm`. Defaults to 'spatial_raw'.
        palette (Optional[list], optional): Palette. Defaults to None.
        point_size (float, optional): Size of point. Defaults to 45.
        save (bool, optional): Whether to save. Defaults to False.
        save_dir (str, optional): Directory to save. Defaults to ''.
        kwargs (dict, optional): Other parameters. Defaults to {}.
    """

    adata_sub, prop = _get_prop(adata)
    coord = adata_sub.obsm[spatial_key]
    prop = prop.div(prop.sum(axis=1), axis=0)
    celltype = list(prop.columns)

    length = len(celltype)
    if palette is None:
        palette = _set_palette(length=length)
        palette = palette[0:length]
    else:
        if len(palette) != length:
            raise ValueError(f"{len(celltype)} cell types exist, please provide a sufficient number of colors")

    palette_use = {}
    for i in range(len(palette)):
        R, G, B = _hex_to_rgb(palette[i])
        palette_use[celltype[i]] = (R, G, B, 1.0)

    ratios = prop[celltype].to_records(index=False).tolist()
    colors = [palette_use[cat] for cat in celltype]

    # plot
    fig, ax = plt.subplots()
    # make pie marker for each unique set of ratios
    df = pd.DataFrame({'x': coord[:, 0], 'y': coord[:, 1], 'ratios': ratios})
    df.ratios = df.ratios.apply(tuple)
    gb = df.groupby("ratios")
    for ratio in gb.groups:
        group = gb.get_group(ratio)
        xys, ss = _pie_marker(ratio)
        for xy, s, color in zip(xys, ss, colors):
            # plot non-zero slices
            if s != 0:
                ax.scatter(group.x, group.y, marker=xy, s=[s * s * point_size], facecolor=color, **kwargs)

    handles = [plt.scatter([], [], color=palette_use[i], label=i) for i in celltype]
    ax.legend(handles=handles, bbox_to_anchor=(1, 0.5), loc='center left',
              ncol=(1 if len(celltype) <= 14 else 2 if len(celltype) <= 30 else 3), borderaxespad=0, frameon=False)
    # ax.axis('off')

    if save:
        plt.savefig(save_dir, format='svg')


def spot_scatterplot(
        adata: AnnData,
        show: Union[str, list],
        spatial_key: str = 'spatial_raw',
        cmap: str = 'magma',
        save: bool = False,
        save_dir: str = '',
        kwargs: dict = {}):
    """
    Scatter plot for single cell type abundance on the spot coordinates.

    Args:
        adata (AnnData): AnnData object of scPositioner output.
        show (Union[str, list]): Cell type to plot. e.g., ['Cell type A', 'Cell type B', 'Cell type C', ...]. Different cell types are displayed separately.
        spatial_key (str, optional): The key of spatial coordinates in `adata.obsm`. Defaults to 'spatial_raw'.
        cmap (str, optional): Cmap. Defaults to 'magma'.
        save (bool, optional): Whether to save. Defaults to False.
        save_dir (str, optional): Directory to save. Defaults to ''.
        kwargs (dict, optional): Other parameters. Defaults to {}.
    """
    adata_sub, _ = _get_prop(adata)
    sc.pl.embedding(adata_sub, basis=spatial_key, color=show, cmap=cmap, show=False, **kwargs)
    if save:
        plt.savefig(save_dir, format='svg')


def spot_scatterplot_multi(
        adata: AnnData,
        show: list,
        spatial_key: str = 'spatial_raw',
        palette: Optional[list] = None,
        white_spacing: float = 20.0,
        colorbar_tick_size: float = 10.0,
        save: bool = False,
        save_dir: str = '',
        kwargs: dict = {}):
    """
    Scatter plot for multi cell type abundance on the spot coordinates.

    Args:
        adata (AnnData): AnnData object of scPositioner output.
        show (list): Cell type to plot. e.g., ['Cell type A', 'Cell type B', 'Cell type C', ...]
        spatial_key (str, optional): The key of spatial coordinates in `adata.obsm`. Defaults to 'spatial_raw'.
        palette (Optional[list], optional): Palette. Defaults to None.
        white_spacing (float, optional): Threshold of cell counts, spots below this value are visualized as white. Defaults to 20.0.
        colorbar_tick_size (float, optional): Font size on the colorbar. Defaults to 10.0.
        save (bool, optional): Whether to save. Defaults to False.
        save_dir (str, optional): Directory to save. Defaults to ''.
        kwargs (dict, optional): Other parameters. Defaults to {}.
    """
    adata_sub, prop = _get_prop(adata)
    prop = prop[show]
    coord = adata_sub.obsm[spatial_key]

    length = len(show)
    if palette is None:
        palette = _set_palette(length=length)
    else:
        if len(palette) != length:
            raise ValueError("The number of provided palettes should be same as the number of cell types to be showed")

    rgb_select = []
    for i in range(length):
        R, G, B = _hex_to_rgb(palette[i])
        rgb_select.append([R, G, B])

    cmaps = []
    for i in range(len(rgb_select)):
        R = int(rgb_select[i][0] * 255)
        G = int(rgb_select[i][1] * 255)
        B = int(rgb_select[i][2] * 255)
        cmap_tmp = _create_colormap(R, G, B, white_spacing=white_spacing)
        cmaps.append(cmap_tmp)

    # plot
    fig = plt.figure()
    colorbar_grid = (length, 1)
    shape = {"vertical_gaps": 2.0, "horizontal_gaps": 0.1, "width": 0.3, "height": 0.2}

    gs = GridSpec(
        nrows=colorbar_grid[0] + 2,
        ncols=colorbar_grid[1] + 1,
        width_ratios=[1, *[shape["width"]] * colorbar_grid[1]],
        height_ratios=[1, *[shape["height"]] * colorbar_grid[0], 1],
        hspace=shape["vertical_gaps"],
        wspace=shape["horizontal_gaps"],
    )
    ax = fig.add_subplot(gs[:, 0], aspect="equal", rasterized=True)

    cbar_axes = []
    for row in range(1, colorbar_grid[0] + 1):
        for column in range(1, colorbar_grid[1] + 1):
            cbar_axes.append(fig.add_subplot(gs[row, column]))

    n_excess = colorbar_grid[0] * colorbar_grid[1] - length
    if n_excess > 0:
        for i in range(1, n_excess + 1):
            cbar_axes[-i].set_visible(False)

    counts = prop.values.copy()

    # plot spots as circles
    c_ord = list(np.arange(0, counts.shape[1]))

    colors = np.zeros((*counts.shape, 4))
    weights = np.zeros(counts.shape)

    for c in c_ord:
        min_color_intensity = counts[:, c].min()
        max_color_intensity = counts[:, c].max()

        rgb_function = _get_rgb_function(cmap=cmaps[c], min_value=min_color_intensity, max_value=max_color_intensity)

        color = rgb_function(counts[:, c])

        norm = matplotlib.colors.Normalize(vmin=min_color_intensity, vmax=max_color_intensity)

        cbar_ticks = [
            min_color_intensity,
            np.mean([min_color_intensity, max_color_intensity]),
            max_color_intensity,
        ]
        cbar_ticks = np.array(cbar_ticks)

        if max_color_intensity > 13:
            cbar_ticks = cbar_ticks.astype(np.int32)
        else:
            cbar_ticks = cbar_ticks.round(2)

        cbar = fig.colorbar(
            matplotlib.cm.ScalarMappable(norm=norm, cmap=cmaps[c]),
            cax=cbar_axes[c],
            orientation="horizontal",
            extend="both",
            ticks=cbar_ticks,
        )

        cbar.ax.tick_params(labelsize=colorbar_tick_size)
        max_color = rgb_function(max_color_intensity / 1.5)
        cbar.ax.set_title(show[c], **{**{"size": colorbar_tick_size, "color": max_color, "alpha": 1}})

        colors[:, c] = color
        weights[:, c] = np.clip(counts[:, c] / (max_color_intensity + 1e-10), 0, 1)
        weights[:, c][counts[:, c] < min_color_intensity] = 0

    colors_ryb = np.zeros((*weights.shape, 3))

    for i in range(colors.shape[0]):
        colors_ryb[i] = _rgb_to_ryb(colors[i, :, :3])

    def kernel(w):
        return w ** 2

    kernel_weights = kernel(weights[:, :, np.newaxis])
    weighted_colors_ryb = (colors_ryb * kernel_weights).sum(axis=1) / kernel_weights.sum(axis=1)

    weighted_colors = np.zeros((weights.shape[0], 4))
    weighted_colors[:, :3] = _ryb_to_rgb(weighted_colors_ryb)
    weighted_colors[:, 3] = colors[:, :, 3].max(axis=1)

    ax.scatter(x=coord[:, 0], y=coord[:, 1], c=weighted_colors, **kwargs)
    # ax.axis('off')

    if save:
        plt.savefig(save_dir, format='svg')


def dist_regplot(
        adata: AnnData,
        save: bool = False,
        save_dir: str = '',
        kwargs: dict = {}):
    """
    Regplot for distance results. This function visualize the mean Euclidean distance of each cell type in select_celltype to background_celltype.

    Args:
        adata (AnnData): AnnData object of scPositioner output. The distance result is stored in `adata.uns['distance_results']`.
        save (bool, optional): Whether to save. Defaults to False.
        save_dir (str, optional): Directory to save. Defaults to ''.
        kwargs (dict, optional): Other parameters. Defaults to {}.
    """
    res = adata.uns['distance_results'].copy()
    mean_distance_df = res['mean_distance_df']
    select_celltype = res['select_celltype']
    background_celltype = res['background_celltype']
    order_list = res['order_list']

    r = res['r']
    p_value = res['p_value']

    title = 'Mean distance to Cell type: ' + str(background_celltype)

    # plot
    ax = sns.regplot(data=mean_distance_df, x="Distribution order", y="Mean distance to background", scatter_kws={'color': 'black'}, **kwargs)
    ax.set_xticks(order_list)
    ax.set_xticklabels(select_celltype)
    plt.xticks(rotation=90)

    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Cell type")
    ax.set_ylabel("Mean distance")

    if r is not None and p_value is not None:
        r = 'Pearson correlation r: ' + str(np.round(r, 2))
        p_value = 'p-value: ' + "{0:.1e}".format(p_value)

        ax.text(x=0.05, y=0.95, s=r, transform=ax.transAxes, ha='left', va='top')
        ax.text(x=0.05, y=0.85, s=p_value, transform=ax.transAxes, ha='left', va='top')

    if save:
        plt.savefig(save_dir, format='svg')


def dist_boxplot(
        adata: AnnData,
        palette: Optional[list] = None,
        save: bool = False,
        save_dir: str = '',
        kwargs: dict = {}):
    """
    Boxplot for distance results. This function visualize the mean Euclidean distance of each cell in select_celltype to background_celltype.

    Args:
        adata (AnnData): AnnData object of scPositioner output. The distance result is stored in `adata.uns['distance_results']`.
        palette (Optional[list], optional): Palette. Defaults to None.
        save (bool, optional): Whether to save. Defaults to False.
        save_dir (str, optional): Directory to save. Defaults to ''.
        kwargs (dict, optional): Other parameters. Defaults to {}.
    """
    res = adata.uns['distance_results'].copy()
    distance_cell_list = res['distance_cell_list']
    select_celltype = res['select_celltype']
    background_celltype = res['background_celltype']

    title = 'Mean distance to Cell type: ' + str(background_celltype)

    length = len(distance_cell_list)
    if palette is None:
        palette = _set_palette(length=length)

    # plot
    ax = sns.boxplot(data=distance_cell_list, showfliers=False, gap=0.2, linecolor='black', showmeans=True, palette=palette, **kwargs)
    ax.set_xticklabels(select_celltype)
    plt.xticks(rotation=90)

    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Cell type")
    ax.set_ylabel("Mean distance")

    if save:
        plt.savefig(save_dir, format='svg')


def gsea_plot(
        adata: AnnData,
        term: str,
        color: Optional[str] = None,
        save: bool = False,
        save_dir: str = '',
        kwargs: dict = {}):
    """
    Enrichment score plot for Gene Set Enrichment Analysis (GSEA) results

    Args:
        adata (AnnData): AnnData object of scPositioner output. The GSEA result is stored in `adata.uns['GSEA_results']`.
        term (str): Name of geneset to plot.
        color (Optional[str], optional): Color. Defaults to None.
        save (bool, optional): Whether to save. Defaults to False.
        save_dir (str, optional): Directory to save. Defaults to ''.
        kwargs (dict, optional): Other parameters. Defaults to {}.
    """
    res = adata.uns['GSEA_results'].copy()
    GSEA_obj = res['GSEA_obj']
    GSEA_obj_df = GSEA_obj.res2d
    compare = res['compare']

    hit = GSEA_obj._results['gsea'][term]['hits']
    RES = GSEA_obj._results['gsea'][term]['RES']
    RES = np.array(RES)
    x_length = np.arange(len(RES))

    record = GSEA_obj_df.iloc[0]
    NES = record['NES']
    p_value = record["NOM p-val"]

    NES = 'NES: ' + "{:.3f}".format(float(NES))
    p_value = 'p-value: ' + "{:.3e}".format(float(p_value))

    title = term + ' (' + compare + ')'
    color = "#88C544" if color is None else color

    left = 0.1
    width = 0.8
    bottom = 0.1
    height = 0

    stat_height_ratio = 0.4
    hits_height_ratio = 0.05

    base = 0.8 / (stat_height_ratio + hits_height_ratio)

    height = hits_height_ratio * base

    fig = plt.figure()
    # Gene hits plot
    ax1 = fig.add_axes([left, bottom, width, height])
    # the x coords of this transformation are data, and the y coord are axes
    trans1 = transforms.blended_transform_factory(ax1.transData, ax1.transAxes)
    # to make axes shared with same x cooridincates, make the vlines same ranges to x
    ax1.vlines(
        [x_length[0], x_length[-1]],
        0,
        1,
        linewidth=0.5,
        transform=trans1,
        color="white",
        alpha=0,  # alpha 0 to transparency
        )
    # add hits line
    ax1.vlines(
        hit, 0, 1, linewidth=0.5, transform=trans1, color="black"
    )

    ax1.tick_params(
        axis="both",
        which="both",
        bottom=True,
        top=False,
        right=False,
        left=False,
        labelbottom=True,
        labelleft=False,
    )
    ax1.set_xlabel("Gene Rank")
    ax1.spines["bottom"].set_visible(True)

    # Enrichment score plot
    bottom += height
    height = stat_height_ratio * base
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.plot(x_length, RES, linewidth=2.5, color=color)
    ax2.text(0.1, 0.2, p_value, transform=ax2.transAxes)
    ax2.text(0.1, 0.3, NES, transform=ax2.transAxes)

    # the y coords of this transformation are data, and the x coord are axes
    trans2 = transforms.blended_transform_factory(ax2.transAxes, ax2.transData)
    ax2.hlines(0, 0, 1, linewidth=1, transform=trans2, color="grey")
    ax2.set_ylabel("Enrichment Score")
    # ax2.set_xlim(min(self._x), max(self._x))
    ax2.tick_params(
        axis="both",
        which="both",
        bottom=False,
        top=False,
        right=False,
        labelbottom=False,
    )
    ax2.locator_params(axis="y", nbins=5)
    # FuncFormatter need two argument, I don't know why. this lambda function used to format yaxis tick labels.
    ax2.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda tick_loc, tick_num: "{:.1f}".format(tick_loc))
    )

    plt.title(title)

    if save:
        plt.savefig(save_dir, format='svg')


def score_boxplot(
        adata: AnnData,
        palette: Optional[list] = None,
        save: bool = False,
        save_dir: str = '',
        kwargs: dict = {}):
    """
    Boxplot for gene signature score results.

    Args:
        adata (AnnData): AnnData object of scPositioner output. The gene signature score result is stored in `adata.uns['score_results']`.
        palette (Optional[list], optional): Palette. Defaults to None.
        save (bool, optional): Whether to save. Defaults to False.
        save_dir (str, optional): Directory to save. Defaults to ''.
        kwargs (dict, optional): Other parameters. Defaults to {}.
    """
    res = adata.uns['score_results'].copy()
    score_df = res['score_df']

    method = res['method']
    p_value = res['p_value']

    title = 'Gene signature score'

    length = len(set(score_df['Group']))
    if palette is None:
        palette = _set_palette(length=length)

    # plot
    ax = sns.boxplot(data=score_df, x="Group", y="Score", hue="Group", showfliers=False, gap=0.2, linecolor='black', showmeans=True, palette=palette, **kwargs)
    plt.xticks(rotation=90)

    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Group")
    ax.set_ylabel("Score")

    if method is not None and p_value is not None:
        p_value = 'p-value: ' + "{0:.2e}".format(p_value)

        ax.text(x=0.05, y=0.95, s=method, transform=ax.transAxes, ha='left', va='top')
        ax.text(x=0.05, y=0.85, s=p_value, transform=ax.transAxes, ha='left', va='top')

    if save:
        plt.savefig(save_dir, format='svg')
