import multiprocessing
from functools import partial

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

from src.CellAbundance import CellAbundance
from src.utils.equations import alfa2rgb, alfa2color, color_mapper

colors = ['#629563', '#044E75', '#CA8F04', '#645D0D',
          '#43BC52', '#B25E89', '#2E3790', '#F118BE',
          '#50974E', '#3273D6', '#0AF24B', '#A3F159',
          '#933835', '#CEB134', '#226BCF', '#856218']


def is_in_square(x, y, x_min, x_max, y_min, y_max):
    return x_min < x < x_max and y_min < y < y_max


def get_segmentation_matrix(data, cell_types, pca_obj, archetype_obj, color_fun, h=800, w=800, radius=100, granularity=200):
    m = np.zeros((h, w, 3), dtype='uint8')
    x_span = np.arange(0, h+granularity, granularity)
    y_span = np.arange(0, w+granularity, granularity)

    y = data['x'].to_numpy()
    x = data['y'].to_numpy()
    t = data['cell_type'].to_numpy()

    for x_i in range(len(x_span)-1):
        for y_i in range(len(y_span)-1):
            site = []
            for c_idx in range(len(x)):
                x_center = (x_span[x_i+1] + x_span[x_i]) / 2
                y_center = (y_span[y_i+1] + y_span[y_i]) / 2
                if CellAbundance.is_in_cirle(x[c_idx], y[c_idx], x_center, y_center, radius):
                    site.append((x[c_idx], y[c_idx], t[c_idx]))
            site = np.array(site)
            counts = CellAbundance.calculate_cells_count(site, cell_types)

            new_pc = pca_obj.transform(counts.reshape(1, -1))
            _, alfa = archetype_obj.transform(new_pc[:, :3])
            color_submatrix = np.array([color_fun(alfa).T for _ in range(granularity)])
            m[x_span[x_i]:x_span[x_i+1], y_span[y_i]:y_span[y_i+1], :] = color_submatrix

    return m


def site_generation_and_count(x, y, t, x_span, y_span, radius, grid_point):
    site = []
    for c_idx in range(len(x)):
        x_center = (x_span[grid_point[0]+1] + x_span[grid_point[0]]) / 2
        y_center = (y_span[grid_point[1]+1] + y_span[grid_point[1]]) / 2
        if CellAbundance.is_in_cirle(x[c_idx], y[c_idx], x_center, y_center, radius):
            site.append((x[c_idx], y[c_idx], t[c_idx]))
    site = np.array(site)

    return site, grid_point


def get_segmentation_matrix_parallel(data, cell_types, pca_obj, archetype_obj, color_fun, h=800, w=800, radius=100, granularity=200):
    m = np.zeros((h, w, 3), dtype='uint8')
    x_span = np.arange(0, h+granularity, granularity)
    y_span = np.arange(0, w+granularity, granularity)

    y = data['x'].to_numpy()
    x = data['y'].to_numpy()
    t = data['cell_type'].to_numpy()

    mesh = np.array(np.meshgrid(range(len(x_span)-1), range(len(y_span)-1)))
    combinations = list(mesh.T.reshape(-1, 2))
    sites_pool = multiprocessing.Pool(20)
    f = partial(site_generation_and_count, x, y, t, x_span, y_span, radius)
    results = sites_pool.map(f, combinations)

    for site, grid_point in results:
        if len(site) > 0:
            counts = CellAbundance.calculate_cells_count(site, cell_types)
        else:
            counts = np.zeros(len(cell_types))
        new_pc = pca_obj.transform(counts.reshape(1, -1))
        _, alfa = archetype_obj.transform(new_pc[:, :3])
        color_submatrix = np.tile(color_fun(alfa[:, 0]), granularity * granularity).reshape((granularity, granularity, -1))
        m[x_span[grid_point[0]]:x_span[grid_point[0]+1], y_span[grid_point[1]]:y_span[grid_point[1]+1], :] = color_submatrix

    return m


def plot_cells_positions(data, cell_types, segment_image=False, segmentation_type='hard', granularity=25, pca_obj=None, AA_obj=None, to_plot='all'):
    groups = data.groupby('cell_type')

    plt.figure(figsize=(10, 10))
    if segment_image is True:
        if pca_obj is None or AA_obj is None:
            raise ValueError("To segment the image pca and archetypes objects are needed")
        if segmentation_type == 'hard':
            color_fun = alfa2color
        else:
            color_vector = np.array([[255, 0, 0], [0, 255, 0], [0, 0, 255], [255, 255, 0]]).T
            color_fun = partial(color_mapper, color_vector)
        z = get_segmentation_matrix_parallel(data, cell_types, pca_obj, AA_obj, color_fun, granularity=granularity)
        plt.imshow(z, origin='lower')

    for (name, group), col in zip(groups, colors):
        if to_plot == 'all' or name in to_plot:
            plt.scatter(group['x'], group['y'], marker="o", s=5, label=name, c=col)

    plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
    plt.show()


def plot_all_cells_positions(patient_ids, cell_types, segment_image=False, segmentation_type='hard', granularity=25, pca_obj=None, AA_obj=None, to_plot='all'):
    plt.figure(figsize=(30, 50))

    for i, patientID in enumerate(patient_ids):
        plt.subplot(8, 5, i+1)
        data = pd.read_csv("../../output/cell_positions_data/patient{}_cell_positions.csv".format(patientID))
        groups = data.groupby('cell_type')

        if segment_image is True:
            if pca_obj is None or AA_obj is None:
                raise ValueError("To segment the image pca and archetypes objects are needed")
            if segmentation_type == 'hard':
                color_fun = alfa2color
            else:
                color_vector = np.array([[255, 0, 0], [0, 255, 0], [0, 0, 255], [255, 255, 0]]).T
                color_fun = partial(color_mapper, color_vector)
            z = get_segmentation_matrix_parallel(data, cell_types, pca_obj, AA_obj, color_fun, granularity=granularity)
            plt.imshow(z, origin='lower')

        for (name, group), col in zip(groups, colors):
            if to_plot == 'all' or name in to_plot:
                plt.scatter(group['x'], group['y'], marker="o", s=5, label=name, c=col)


    #plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
    plt.show()


def plot_scatter_pca(principal_components, evr, labels=None, original_axis=None, archetypes=None, cell_type=None):
    fig, ax = plt.subplots(figsize=(10,10))

    if labels is not None:
        ax.scatter(principal_components[:, 0], principal_components[:, 1],
                   color=labels, alpha=0.9, cmap=plt.cm.nipy_spectral)
    else:
        ax.scatter(principal_components[:, 0], principal_components[:, 1], alpha=0.5)

    if original_axis is not None:
        origin = original_axis[len(original_axis)-1, :]

        red_patch = []
        colors = 'bgrcmyk'
        for i in range(len(original_axis)-1):
            plt.plot([origin[0], original_axis[i, 0]], [origin[1], original_axis[i, 1]], 'k-', color='b', alpha=0.5)
            #plt.arrow(origin[0], origin[1], original_axis[i, 0], original_axis[i, 1], color=colors[i%len(colors)])
            plt.text(original_axis[i, 0], original_axis[i, 1], i, color='r', ha='center', va='center')
            red_patch.append(mpatches.Patch(color=colors[i%len(colors)], label="{} = {}".format(i, cell_type[i])))
            plt.legend(handles=red_patch, title='Cell Type', bbox_to_anchor=(1.0, 1), loc='upper left')

    if archetypes is not None:
        plt.scatter(archetypes[0,:], archetypes[1,:], marker="^", s=100, color='orange')


    tot_evr = np.sum(evr[0:2])
    plt.title("{:.2f}% Total Exp. Var.".format(tot_evr*100))
    ax.set_xlabel("PC1 ({:.2f}% Exp. Var.)".format(evr[0] * 100))
    ax.set_ylabel("PC2 ({:.2f}% Exp. Var.)".format(evr[1] * 100))
    if original_axis is None:
        minimum, maximum = _calculate_pca_max_min(principal_components[:, 0:2])
        plt.xlim(minimum, maximum)
        plt.ylim(minimum, maximum)
    plt.show()


def plot_scatter_NNMF(nnmf_features, labels=None):
    plt.figure(figsize=(10, 10))
    xs = nnmf_features[:, 0]
    ys = nnmf_features[:, 1]

    minimum, maximum = _calculate_pca_max_min(nnmf_features[:, 0:2])
    plt.xlim(minimum, maximum)
    plt.ylim(minimum, maximum)
    # Scatter plot plt.scatter(xs, ys, alpha=0.5)
    # Annotate the points
    if labels:
        plt.scatter(xs, ys, c=labels, cmap=plt.cm.nipy_spectral, alpha=0.5)
    else:
        plt.scatter(xs, ys, alpha=0.5)

    plt.show()


def plot_3Dscatter_pca(principal_components, evr, labels, archetypes=None):
    tot_evr = np.sum(evr[0:3])
    print("{:.2f}% Total Exp. Var.".format(tot_evr))
    fig = plt.figure(1, figsize=(8, 8))
    plt.clf()
    ax = Axes3D(fig, rect=(0, 0, 1, 1), elev=30, azim=134)

    minimum, maximum = _calculate_pca_max_min(principal_components[:, 0:3])
    ax.set_xlim(minimum, maximum)
    ax.set_ylim(minimum, maximum)
    ax.set_zlim(minimum, maximum)
    ax.set_xlabel("PC1 ({:.2f}% Exp. Var.)".format(evr[0] * 100))
    ax.set_ylabel("PC2 ({:.2f}% Exp. Var.)".format(evr[1] * 100))
    ax.set_zlabel("PC3 ({:.2f}% Exp. Var.)".format(evr[2] * 100))

    ax.scatter(principal_components[:, 0], principal_components[:, 1], principal_components[:, 2],
               c=labels, cmap=plt.cm.nipy_spectral, edgecolor='k')

    if archetypes is not None:
        ax.scatter(archetypes[0, :], archetypes[1, :], archetypes[2, :],
                   marker = "^", s = 100, color = 'orange')

    plt.show()

def plot_3Dscatter_NNMF(nnmf_features, labels=None):
    fig = plt.figure(1, figsize=(8, 8))
    plt.clf()
    ax = Axes3D(fig, rect=(0, 0, 1, 1), elev=30, azim=134)

    minimum, maximum = _calculate_pca_max_min(nnmf_features[:, 0:3])
    ax.set_xlim(minimum, maximum)
    ax.set_ylim(minimum, maximum)
    ax.set_zlim(minimum, maximum)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_zlabel("PC3")

    ax.scatter(nnmf_features[:, 0], nnmf_features[:, 1], nnmf_features[:, 2],
               c=labels, cmap=plt.cm.nipy_spectral, edgecolor='k')

    plt.show()


def biplot_PCA(score, coeff, labels=None, cell_type=None):
    plt.figure(figsize=(10, 10))

    minimum, maximum = _calculate_pca_max_min(score[:, 0:2])
    plt.xlim(minimum, maximum)
    plt.ylim(minimum, maximum)
    n = coeff.shape[0]
    red_patch = []
    for i in range(n):
        plt.arrow(0, 0, coeff[i, 0], coeff[i, 1], color='b', alpha=0.5)
        plt.text(coeff[i, 0] * 1.15, coeff[i, 1] * 1.15, i, color='r', ha='center', va='center')
        red_patch.append(mpatches.Patch(color='red', label="{} = {}".format(i, cell_type[i])))

    if labels is not None:
        plt.scatter(score[:, 0], score[:, 1],
                   c=labels, cmap=plt.cm.nipy_spectral, alpha=0.5)
    else:
        plt.scatter(score[:, 0], score[:, 1], alpha=0.5)

    plt.legend(handles=red_patch, title='Cell Type', bbox_to_anchor=(1.0, 1), loc='upper left')
    plt.xlabel("PC{}".format(1))
    plt.ylabel("PC{}".format(2))
    plt.grid()


def plot_cumulative_explained_variance(evr):
    plt.plot(range(len(np.cumsum(evr))), np.cumsum(evr))
    v = np.argmax(np.cumsum(evr) > 0.95)
    plt.axvline(x=v, color='red', linestyle='--')
    plt.ylim(0, 1)
    plt.xlabel('number of components')
    plt.ylabel('cumulative explained variance')
    plt.show()


def plot_PCAs(cell_ab_list, scale_pca=True, fix_axis=False):
    plt.figure(figsize=(30, 50))

    pca_list = [ca.perform_PCA(scale=scale_pca)[:, 0:2] for ca in cell_ab_list]
    minimum, maximum = 0, 0
    if fix_axis:
        minimum, maximum = _calculate_pca_max_min(pca_list)

    for i, (ca, pc) in enumerate(zip(cell_ab_list, pca_list)):
        plt.subplot(8, 5, i+1)
        evr = ca.pca.explained_variance_ratio_
        tot_evr = (evr[0] + evr[1])*100
        plt.title("Patient ID: {}, {:.2f}% Total Exp. Var.".format(ca.patient_id, tot_evr))

        if not fix_axis:
            minimum, maximum = _calculate_pca_max_min(pc[:, 0:2])

        plt.ylim(minimum, maximum)
        plt.xlim(minimum, maximum)
        plt.xlabel("PC1 ({:.2f}% Exp. Var.)".format(evr[0] * 100))
        plt.ylabel("PC2 ({:.2f}% Exp. Var.)".format(evr[1] * 100))
        plt.scatter(pc[:, 0], pc[:, 1], alpha=0.5)
    plt.show()


def plot_NNMFs(cell_ab_list, model, fix_axis=False):
    plt.figure(figsize=(30, 50))

    nnmf_list = [model.fit_transform(ca.abundance_matrix)[:, 0:2] for ca in cell_ab_list]
    minimum, maximum = 0, 0
    if fix_axis:
        minimum, maximum = _calculate_pca_max_min(nnmf_list)

    for i, (ca, pc) in enumerate(zip(cell_ab_list, nnmf_list)):
        plt.subplot(8, 5, i+1)
        plt.title("Patient ID: {}".format(ca.patient_id))

        if not fix_axis:
            minimum, maximum = _calculate_pca_max_min(pc[:, 0:2])

        plt.ylim(minimum, maximum)
        plt.xlim(minimum, maximum)
        #plt.xlabel("PC1 ({:.2f}% Exp. Var.)".format(evr[0] * 100))
        #plt.ylabel("PC2 ({:.2f}% Exp. Var.)".format(evr[1] * 100))
        plt.scatter(pc[:, 0], pc[:, 1], alpha=0.5)
    plt.show()


def _calculate_pca_max_min(pca):
    if isinstance(pca, list):
        PC_joined = [j for pc in pca for j in pc.flatten()]
    else:
        PC_joined = pca.flatten()
    return min(PC_joined), max(PC_joined)


def plot_CVEs(cell_ab_list, scale_pca=True):
    plt.figure(figsize=(10, 10))
    v = []
    for i, ca in enumerate(cell_ab_list):
        _ = ca.perform_PCA(scale=scale_pca)
        evr = ca.pca.explained_variance_ratio_
        plt.plot(np.cumsum(evr))
        plt.ylim(0, 1)
        v.append(np.argmax(np.cumsum(evr) > 0.95))
    plt.axvline(x=np.mean(v), color='red', linestyle='--')
    plt.show()


def plot_exp_variance(cell_ab_list, scale_pca=False):
    plt.figure(figsize=(10, 10))
    for i, ca in enumerate(cell_ab_list):
        _ = ca.perform_PCA(scale=scale_pca)
        evr = np.sqrt(ca.pca.explained_variance_)
        #alpha = 1 if evr[0] < 0.5 else 0.1
        lab = "Patient {}".format(ca.patient_id) if evr[0] < 0.5 else None
        plt.plot(evr, label="Patient {}".format(ca.patient_id))
    #plt.ylim(0, 0.8)
    plt.xlim(0, 6)
    plt.xlabel("# of components")
    plt.ylabel("Exp. Variance")
    plt.legend(title='Patients PC1 Expr. Var < 0.5', bbox_to_anchor=(1.0, 1), loc='upper left')
    plt.show()


def plot_stacked_variance_PCA(cell_ab_list, scale_pca=False):

    patients_idx = [ca.patient_id for ca in cell_ab_list]

    evr_matrix = np.zeros((17, len(cell_ab_list)))

    for i, ca in enumerate(cell_ab_list):
        _ = ca.perform_PCA(scale=scale_pca)
        evr = ca.pca.explained_variance_
        evr_matrix[:, i] = evr

    data = {"PC{}".format(i+1): evr_matrix[i, :] for i in range(17)}
    df = pd.DataFrame(data)

    index = pd.Index(patients_idx, name='Patient ID')
    df = pd.DataFrame(data, index=index)
    s = df.sort_values(by=['PC1', 'PC2'], ascending=False)

    ax = s.plot(kind='bar', stacked=True, figsize=(20, 10))
    plt.axhline(y=0.4, color='r', linestyle='-')
    plt.axhline(y=0.5, color='r', linestyle='-')
    plt.axhline(y=0.6, color='r', linestyle='-')
    ax.set_ylabel('Stacked PC')
    plt.legend(title='labels', bbox_to_anchor=(1.0, 1), loc='upper left')


def plot_cev_radius(cumulative_expr_dict, patient):
    plt.figure(figsize=(10, 10))
    for k, cev in cumulative_expr_dict.items():
        plt.plot(cev, label=k)
    plt.xlabel('number of components')
    plt.ylabel('cumulative explained variance')
    plt.title("Patient {}".format(patient))
    plt.legend(title='labels', bbox_to_anchor=(1.0, 1), loc='upper left')
    plt.show()


def plot_stacked_var_radius(cumulative_expr_dict):
    radius = [r for r in cumulative_expr_dict]

    evr_matrix = np.zeros((17, len(radius)))

    for i, r in enumerate(radius):
        evr = cumulative_expr_dict[r]
        evr_matrix[:, i] = evr

    data = {"PC{}".format(i + 1): evr_matrix[i, :] for i in range(17)}
    df = pd.DataFrame(data)

    index = pd.Index(radius, name='Radius')
    df = pd.DataFrame(data, index=index)
    #s = df.sort_values(by=['PC1', 'PC2'], ascending=False)
    s = df

    ax = s.plot(kind='bar', stacked=True, figsize=(20, 10))
    #plt.axhline(y=0.4, color='r', linestyle='-')
    #plt.axhline(y=0.5, color='r', linestyle='-')
    #plt.axhline(y=0.6, color='r', linestyle='-')
    ax.set_ylabel('Stacked PC')
    plt.legend(title='labels', bbox_to_anchor=(1.0, 1), loc='upper left')


def tumor_misfit_barplot(tumor_id, errors):
    plt.figure(figsize=(10, 5))
    y_pos = np.arange(len(tumor_id))
    plt.bar(y_pos, errors)
    plt.xticks(y_pos, tumor_id)
    plt.xlabel('Patient ID')
    plt.ylabel("Relative Error")
    plt.title("Relative misfit error for every patient")
    #plt.ylim(0, 1)

    plt.show()


def archetypes_bar_plot(cell_number_archetypes, cell_types, colors, y_axis='count', radius=100):
    if y_axis not in ['log', 'count', 'density']:
        raise ValueError("Wrong parameter: y_axis paramenter must be log, count or density")

    fig = plt.figure(figsize=(20, 10))
    y_pos = 3*np.arange(len(cell_number_archetypes[0]))
    ax = fig.add_axes([0, 0, 1, 1])
    width = 0.50

    data = cell_number_archetypes
    if y_axis == 'log':
        data = [np.log(d) for d in cell_number_archetypes]
    elif y_axis == 'density':
        data = [100 * (d / (radius*radius*np.pi)) for d in cell_number_archetypes]

    ax.bar(y_pos - 2*width, data[0], color=colors[0], width=width, label="Arch 0")
    ax.bar(y_pos - width, data[1], color=colors[1], width=width, label="Arch 1")
    ax.bar(y_pos, data[2], color=colors[2], width=width, label="Arch 2")
    ax.bar(y_pos + width, data[3], color=colors[3], width=width, label="Arch 3")

    plt.legend()

    plt.xticks(y_pos, cell_types, rotation=45)
    plt.show()


def archetype_simple_plot(cell_number_archetypes, archetype_id, colors, cell_types):
    plt.figure(figsize=(20, 5))
    y_pos = np.arange(len(cell_number_archetypes))
    plt.bar(y_pos, cell_number_archetypes, color=colors[archetype_id])
    plt.xticks(y_pos, cell_types, rotation=45)

    plt.xlabel('Cells Type')
    plt.ylabel("Count")
    plt.title("Archetype {} cells count".format(archetype_id))
