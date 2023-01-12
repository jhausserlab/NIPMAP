import multiprocessing
from functools import partial

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
from scipy.spatial import ConvexHull
import mpl_toolkits.mplot3d as a3
import matplotlib as mpl
import matplotlib.colors
import numpy as np
import seaborn as sns
#import image_seaborn as isns
import pandas as pd
import matplotlib.patheffects as path_effects

from src.CellAbundance import CellAbundance
from src.utils.equations import alfa2rgb, alfa2color, color_mapper

colors = ['#629563', '#044E75', '#CA8F04', '#645D0D',
          '#43BC52', '#B25E89', '#2E3790', '#F118BE',
          '#50974E', '#3273D6', '#0AF24B', '#A3F159',
          '#933835', '#CEB134', '#226BCF', '#856218','#831CCB']


def is_in_square(x, y, x_min, x_max, y_min, y_max):
    return x_min < x < x_max and y_min < y < y_max


def dist_from_borders(x, y, h, w):
    return x, h-x, y, w-y


def calculate_outside_area(x_center, y_center, radius, h, w, circle_out_area_dict):
    x_left, x_right, y_bottom, y_upper = dist_from_borders(x_center, y_center, h, w)
    outsite_area = 0
    if min(x_left, x_right) < radius and min(y_bottom, y_upper) < radius:
        # is an angle
        try:
            outsite_area = circle_out_area_dict[(min(x_left, x_right), min(y_bottom, y_upper))]
        except KeyError:
            outsite_area = calculate_out_area(x_center, y_center, radius, 'angle')
            circle_out_area_dict[(min(x_left, x_right), min(y_bottom, y_upper))] = outsite_area
    elif min(x_left, x_right, y_bottom, y_upper) < radius:
        try:
            outsite_area = circle_out_area_dict[min(x_left, x_right, y_bottom, y_upper)]
        except KeyError:
            # print(x_center, y_center, radius)
            outsite_area = calculate_out_area(400, min(x_left, x_right, y_bottom, y_upper), radius, 'segment')
            circle_out_area_dict[min(x_left, x_right, y_bottom, y_upper)] = outsite_area

    return outsite_area


def site_generation(x, y, t, x_span, y_span, radius, grid_point, h, w, circle_out_area_dict):
    """
    Generate the sites using a circle with the center corresponding in the center of the square.
    """
    x_center = (x_span[grid_point[0]+1] + x_span[grid_point[0]]) / 2
    y_center = (y_span[grid_point[1]+1] + y_span[grid_point[1]]) / 2

    area_outside = calculate_outside_area(x_center, y_center, radius, h, w, circle_out_area_dict)

    idx = np.where((x - x_center) * (x - x_center) + (y - y_center) * (y - y_center) <= radius * radius)
    site = np.array([(x[i], y[i], t[i]) for i in idx])
    return site, area_outside


def calculate_out_area(x_center, y_center, radius, type='segment'):
    if type not in ['segment', 'angle']:
        raise ValueError("wrong type, must be segment or angle")
    area = 0
    if type == 'segment':
        x_1 = x_center + np.sqrt(radius**2 - y_center**2)
        x_2 = x_center - np.sqrt(radius**2 - y_center**2)
        alfa = np.arccos((2*(radius**2) - np.abs(x_2-x_1)**2)/(2*(radius**2)))
        area = 0.5 * (alfa - np.sin(alfa)) * radius**2
    elif type == 'angle':
        n_sample = 10000
        r = np.random.uniform(0, radius, n_sample)
        theta = np.random.uniform(0, 360, n_sample)
        x = r*np.cos(theta) + x_center
        y = r*np.sin(theta) + y_center
        in_points = np.sum((x >= 0) & (y >= 0))
        area = ((n_sample-in_points) * (np.pi * radius ** 2)) / n_sample

    return area


def xfrange(start, stop, step):
    i = 0
    while start + i * step < stop:
        yield start + i * step
        i += 1


def get_segmentation_matrix(data, cell_types, pca_obj, archetype_obj, color_fun, h=800, w=800, counting_type='abs', radius=100, granularity=25):
    """
    Calculate the segmentation matrix for a given dataset.

    @param data: csv containing the data of the sample
    @param cell_types: cell types vector
    @param pca_obj: the sklearn pca trained model
    @param archetype_obj: the archetype analysis trained model
    @param color_fun: color function to be used
    @param h: {int} height of the image in pixels (Default=800)
    @param w: {int} width of the image in pixels (Default=800)
    @param counting_type: counting method of cells within sites (Default='abs')
    @param radius: radius of the sites generated in the images in micrometer (Default=100)
    @param granularity: segmentation granularity (dimension of the square sites)
    
    @return: m {matrix} matrix of colors at each position x,y
    """
    m = np.empty((h, w, 3), dtype='uint8')
    x_span = np.arange(0, h+granularity, granularity)
    y_span = np.arange(0, w+granularity, granularity)
    #print(data)
    y = data['x'].to_numpy()
    x = data['y'].to_numpy()
    t = data['cell_type'].to_numpy()

    mesh = np.array(np.meshgrid(range(len(x_span)-1), range(len(y_span)-1)))
    combinations = list(mesh.T.reshape(-1, 2)) # creating all grid points
    circle_out_area = {}
    print(archetype_obj.alfa.shape[0]-1)
    total_site_area = np.pi * radius ** 2
    for grid_point in combinations:
        site, area_outside = site_generation(x, y, t, x_span, y_span, radius, grid_point, h, w, circle_out_area)

        # based on the counting type we generate the counts matrix for every locations
        if len(site) > 0:
            if counting_type == 'abs':
                counts = CellAbundance.calculate_cells_count(site, cell_types, radius)
            elif counting_type == 'gaussian':
                x_center = (x_span[grid_point[0] + 1] + x_span[grid_point[0]]) / 2
                y_center = (y_span[grid_point[1] + 1] + y_span[grid_point[1]]) / 2
                center = np.vstack((x_center, y_center))
                std_coeff = 2
                counts, _ = CellAbundance.calculate_gaussian_count(site, cell_types, center, std_coeff*radius)
            else:
                raise ValueError("abs or gaussian are the only allowed counting types")
            counts = counts / (1 - (area_outside / total_site_area))
        else:
            counts = np.zeros(len(cell_types))
        #print(pca_obj.n_features_,pca_obj.n_samples_,pca_obj.n_components_)
        new_pc = pca_obj.transform(counts.reshape(1, -1))
        _, alfa = archetype_obj.transform(new_pc[:, :archetype_obj.alfa.shape[0]-1])#transform(new_pc[:, :3])
        # using the alfa given by the archetype model for the given square
        x_min = x_span[grid_point[0]]
        x_max = np.minimum(x_span[grid_point[0]+1], h)
        y_min = y_span[grid_point[1]]
        y_max = np.minimum(y_span[grid_point[1]+1], w)
        x_dim = x_max - x_min
        y_dim = y_max - y_min
        #print(color_fun(alfa[:, 0]))
        color_submatrix = np.tile(
            color_fun(alfa[:, 0]),
            x_dim * y_dim).reshape((x_dim, y_dim, -1))
        m[x_min:x_max, y_min:y_max, :] = color_submatrix
    #print(m.shape)
    return m


def plot_cells_positions(data, cell_types, segment_image=False, segmentation_type='hard', color_vector=None,counting_type='abs', h=800, w=800, granularity=25, radius=25, pca_obj=None, AA_obj=None, to_plot='all',path_fig=None):
    '''
    plots cells positions in MIBI image of TNBC + saves it in a .svg image
    
    @param cell_types: {list} of str, cell type labels in the dataset
    @param segment_image: {boolean}, if True, color image from results of AA, if False, plot only cells (Default=False)
    @param segmentation_type:{str} segment pixels of images and color them as TMENs proportions (default='hard')
    @param color_vector:{list} list of colors in RGB format (0 to 255). Default if matplotlib colormap 'Dark2'.
    @param counting_type:{str} counting type of cells within the sites:'abs' for absolute, 'log' fro log normalization and 'gaussian' for gaussian density (Default='abs')
    @param h:{int} height of image in pixels (Default=800)
    @param w:{int} width of image in pixels (Default=800)
    @param radius: {int} radius of sites generated in images in micrometer (Default=100)
    @param granularity: {int} granularity of colors (Default=25)
    @param pca_obj: {PCA obj} object of PCA on sites generated in the images (Default=None)
    @param AA_obj: {AA object} object of Archetype Analysis
    @param to_plot:{str} plot all images or not (Default='all')
    
    @return: None
    
    '''
    groups = data.groupby('cell_type')
    cells_cols = {cell_types[i]:colors[i] for i in range(len(cell_types))}
    plt.figure(figsize=(8, 8))
    if segment_image is True:
        if pca_obj is None or AA_obj is None:
            raise ValueError("To segment the image pca and archetypes objects are needed")
        
        if color_vector is None:
            colormap = mpl.cm.Dark2.colors
            
            color_vector=np.array(colormap[0:AA_obj.alfa.shape[0]])*255
            print(color_vector)
        else:
            color_vector =  np.array([[255, 0, 223],[255,0,0],[70,203,236],[0,0,0]])#np.array([[255, 0, 0], [0, 153, 51], [0, 0, 255], [255, 255, 0]])
        if segmentation_type == 'hard': #color pixel by 1 of the colors defining TMENs (discrete)
            color_fun = partial(alfa2color, color_vector)

        else: #color image by continuous spectrum of colors, depending on granularity 
            color_fun = partial(color_mapper, color_vector.T)
            
        z = get_segmentation_matrix(data, cell_types, pca_obj, AA_obj, color_fun, counting_type=counting_type, radius=radius, granularity=granularity, h=h, w=w)#get_segmentation_matrix(data, cell_types, pca_obj, AA_obj, color_fun, counting_type=counting_type, radius=radius, granularity=granularity, h=h, w=w)
        plt.imshow(z, origin='lower')

    for (name, group), col in zip(groups, colors):
        if to_plot == 'all' or name in to_plot:
            plt.scatter(group['x'], group['y'], marker="o", s=5, label=name, c=cells_cols[name])

    plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
    plt.xlim(0, w)
    plt.ylim(0, h)
    if path_fig!=None:
        plt.savefig(path_fig,format="svg")
    return None #plt


def plot_cells_markers_tmens(patient_id,cell_types,path_data, data_markers_path,cell_type, marker,symbols,col=None,segment_image=False,
                             segmentation_type="hard", counting_type="gaussian",h=800,w=800,granularity=20,radius=25,
                             pca_obj=None,AA_obj=None, to_plot=None,path_fig=None, intOutQuant=0):
    '''
    plots cells positions in MIBI image of TNBC overlayed with TMENs colors + saves it in a .svg image
    @param patient_id:{integer} id of the image patient (number here)
    @param cell_types: {list} of str, cell type labels in the dataset
    @param path_data: {string} path of data of patient (cells positions and cell types in an image)
    @param data_markers_path:{string} path to csv file containing markers expressions for each cells in all images
    @param cell_type:{string} cell type to display
    @param marker:{string} marker expressed in a cell type
    @param symbols:{list} of symbols for each marker to plot
    @param col:{string} color of a symbol to distinguish with the other symbol
    @param segment_image: {boolean}, if True, color image from results of AA, if False, plot only cells (Default=False)
    @param segmentation_type:{str} ssegment pixels of images and color them as TMENs proportions (default='hard')
    @param counting_type:{str} counting type of cells within the sites:'abs' for absolute, 'log' fro log normalization and 'gaussian' for gaussian density (Default='abs')
    @param h:{int} height of image in pixels (Default=800)
    @param w:{int} width of image in pixels (Default=800)
    @param radius: {int} radius of sites generated in images in micrometer (Default=100)
    @param granularity: {int} granularity of colors (Default=25)
    @param pca_obj: {PCA obj} object of PCA on sites generated in the images (Default=None)
    @param AA_obj: {AA object} object of Archetype Analysis
    @param to_plot:{str} plot all images or not (Default='all')
    
    @return: plot
    '''    
    data = pd.read_csv(path_data+"/patient{}_cell_positions.csv".format(patient_id))
   # data = data.loc[data["cell_type"] == cell_type]
    groups = data.groupby('cell_type')
    print(type(marker)==str)
    df_markers = pd.read_csv(data_markers_path)
    plt.figure(figsize=(8, 8))
    if segment_image is True:
        if pca_obj is None or AA_obj is None:
                raise ValueError("To segment the image pca and archetypes objects are needed")
        color_vector =  np.array([[255, 0, 223],[255,0,0],[70,203,236],[0,0,0]])#np.array([[255, 0, 0], [0, 153, 51], [0, 0, 255], [255, 255, 0]])
        if segmentation_type == 'hard': #color pixel by 1 of the colors defining TMENs (discrete)
                color_fun = partial(alfa2color, color_vector)

        else: #color image by continuous spectrum of colors, depending on granularity 
            color_fun = partial(color_mapper, color_vector.T)

        z = get_segmentation_matrix(data, cell_types, pca_obj, AA_obj, color_fun, counting_type=counting_type, radius=radius, granularity=granularity, h=h, w=w)
        data = data.loc[data["cell_type"] == cell_type]
        df_markers = df_markers.loc[df_markers["SampleID"] == patient_id]
        df_markers.rename(columns = {"cellLabelInImage":"label"},inplace=True)
        data_CM = pd.merge(data,df_markers,on = "label",how = "left")
        print(data.shape)
        #print(data_CM)
        plt.imshow(z, alpha=.5,origin='lower')#isns.imgplot(z)
        #cm = plt.cm.get_cmap("YlGn")#('RdYlBu')
        #cm = plt.cm.get_cmap("summer")#('RdYlBu') #fails on DC Keratin6
        #cm = plt.cm.get_cmap("viridis")#('RdYlBu') #fails on DC Keratin6
        #cm = plt.cm.get_cmap("bone")#('RdYlBu') #fails on DC Keratin6 becuase high Ker6 DCs diseapear in cancer niche
        cm = plt.cm.get_cmap("pink")#('RdYlBu')
        #cm = plt.cm.get_cmap("hot")#('RdYlBu') #works
        #cm = plt.cm.get_cmap("copper")#('RdYlBu') #ok, but DC Ker6 not so visible in cancer
        if type(marker)==str:
            print("ok")
            maxIntensity = data_CM[marker].quantile(1-intOutQuant)
            data_CM.loc[data_CM[marker] > maxIntensity,marker ] = maxIntensity
            plt.scatter(data_CM['x'], data_CM['y'], marker="o", s=50, c=data_CM[marker],cmap=cm)
            plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
            plt.colorbar()
        else:
            #print(data_CM["Keratin6"].tolist())
            plt.scatter(data_CM['x'], data_CM['y'], marker="o",facecolors="none", s=90, edgecolors="white",alpha=.1)
            for m,s in zip(marker,symbols):##
                #print(m)
                namePhen = m+'_phen'
                data_CM[namePhen]="o" # for each marker, if cells are positive to m, set a symbol and plot it
                data_CM.loc[data_CM[m]>0,namePhen] = s
                #if m=="CD45RO":
                #    #print(data_CM[namePhen])
                for i in range(len(data_CM[namePhen].tolist())):
                    if m=="Keratin6"and data_CM[m+'_phen'][i]==s:
                        #print(data_CM[namePhen].tolist())
                        #print(m,s)
                        plt.scatter(data_CM['x'][i],data_CM['y'][i],marker="o",s=100,facecolors='none',edgecolors=col,linewidths=3,label=m,alpha=.8)
                    else:
                        plt.scatter(data_CM['x'][i],data_CM['y'][i],marker=data_CM[namePhen].tolist()[i],s=90,facecolors='none',edgecolors="white",linewidths=1,label=m,alpha=.9)
                #sns.scatterplot(data_CM['x'],data_CM['y'],markers=data_CM[namePhen].tolist(),s=14,color=".2")#
        # if cells are positive to multiple markers, stack the same points with different symbols in the plot       
            
        recs = []
        plt.legend(recs,marker, loc='upper left')#plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
    plt.xlim(0, w)
    plt.ylim(0, h)
    if path_fig!=None:
        plt.savefig(path_fig,format="svg")
    
    
    return plt

def plot_cells_TI_border(data, cell_types, patientID=1,h=800,w=800,
                              data_niches_path ='../../output/archetypes_sites_centered_cells_gaussian2.csv',
                              intf=["arch2","arch3"],normalize=False,quant=0.95,path_fig=None):
    plt.figure(figsize=(8, 8))
    df_niches = pd.read_csv(data_niches_path)
    #print(df_niches["arch2"])
    if normalize==True:
        ### Set fibrotic nice weight to 0
        df_niches.loc[df_niches["arch4"]!=0,"arch4"]=0
        ## Renormalize to 1
        df_niches["arch1"]=df_niches.apply(lambda row: row["arch1"]/(row["arch1"] + row["arch2"] + row["arch3"]),axis=1)
        df_niches["arch2"]=df_niches.apply(lambda row: row["arch2"]/(row["arch1"] + row["arch2"] + row["arch3"]),axis=1)
        df_niches["arch3"]=df_niches.apply(lambda row: row["arch3"]/(row["arch1"] + row["arch2"] + row["arch3"]),axis=1)
        #print(df_niches["arch2"]) 
    
    groups = data.groupby('cell_type')
    cells_cols = {cell_types[i]:colors[i] for i in range(len(cell_types))}
    data.rename(columns = {"label":"site_id"},inplace=True)
    #print(data)
    data_niches = pd.merge(data,df_niches,on="site_id",how="left")
    data_niches=data_niches.loc[data_niches["patient_id"]==patientID]
    
    if len(intf)==2:
        data_niches["interface"]= data_niches.apply(lambda row: row.arch2 * row.arch3,axis=1)
        #print(data_niches["interface"].describe())
    else:
        data_niches["interface"]=data_niches[intf]
        #print(data_niches["interface"].describe())
        
    if normalize==False and quant!=0:  
        #print("Ok")
        maxWeight = data_niches["interface"].quantile(quant)
        print(maxWeight)
        data_niches.loc[data_niches["interface"]> maxWeight,"interface"] = maxWeight
    
    for (name, group), col in zip(groups, colors):
        plt.scatter(group['x'], group['y'], marker="o", s=10, label=name,facecolors='none',edgecolors=cells_cols[name])

    plt.scatter(data_niches['x'], data_niches['y'], marker="o", s=5,
                c=data_niches["interface"],cmap=plt.cm.get_cmap("Greys"))

    plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
    plt.xlim(0, w)
    plt.ylim(0, h)
    plt.colorbar()
    if path_fig!=None:
        plt.savefig(path_fig,format="svg")
    return plt

def plot_cells_alfa_tmens(patient_id,cell_types,path_data, data_markers_path,cell_type, path_pred_alfas,segment_image=False,
                             segmentation_type="hard", counting_type="gaussian",h=800,w=800,granularity=20,radius=25,
                             pca_obj=None,AA_obj=None, to_plot=None,path_fig=None):
    '''
    plots cells positions in MIBI image of TNBC overlayed with TMENs colors + saves it in a .svg image
    @param patient_id:{integer} id of the image patient (number here)
    @param cell_types: {list} of str, cell type labels in the dataset
    @param path_data: {string} path of data of patient (cells positions and cell types in an image)
    @param data_markers_path:{string} path to csv file containing markers expressions for each cells in all images
    @param cell_type:{string} cell type to display
    @param path_pred_alfas:{string} paht to csv file with predicted values of niche alfa of cell type of interest
    @param segment_image: {boolean}, if True, color image from results of AA, if False, plot only cells (Default=False)
    @param segmentation_type:{str} ssegment pixels of images and color them as TMENs proportions (default='hard')
    @param counting_type:{str} counting type of cells within the sites:'abs' for absolute, 'log' fro log normalization and 'gaussian' for gaussian density (Default='abs')
    @param h:{int} height of image in pixels (Default=800)
    @param w:{int} width of image in pixels (Default=800)
    @param radius: {int} radius of sites generated in images in micrometer (Default=100)
    @param granularity: {int} granularity of colors (Default=25)
    @param pca_obj: {PCA obj} object of PCA on sites generated in the images (Default=None)
    @param AA_obj: {AA object} object of Archetype Analysis
    @param to_plot:{str} plot all images or not (Default='all')
    
    @return: plot
    '''    
    data = pd.read_csv(path_data+"/patient{}_cell_positions.csv".format(patient_id))
   # data = data.loc[data["cell_type"] == cell_type]
    groups = data.groupby('cell_type')
    pred_alfas = pd.read_csv(path_pred_alfas)
    df_markers = pd.read_csv(data_markers_path)
    plt.figure(figsize=(8, 8))
    if segment_image is True:
        if pca_obj is None or AA_obj is None:
                raise ValueError("To segment the image pca and archetypes objects are needed")
        color_vector =  np.array([[255, 0, 223],[255,0,0],[70,203,236],[0,0,0]])#np.array([[255, 0, 0], [0, 153, 51], [0, 0, 255], [255, 255, 0]])
        if segmentation_type == 'hard': #color pixel by 1 of the colors defining TMENs (discrete)
                color_fun = partial(alfa2color, color_vector)

        else: #color image by continuous spectrum of colors, depending on granularity 
            color_fun = partial(color_mapper, color_vector.T)

        z = get_segmentation_matrix(data, cell_types, pca_obj, AA_obj, color_fun, counting_type=counting_type, radius=radius, granularity=granularity, h=h, w=w)
        data = data.loc[data["cell_type"] == cell_type]
        df_markers = df_markers.loc[df_markers["SampleID"] == patient_id]
        df_markers.rename(columns = {"cellLabelInImage":"label"},inplace=True)
        data_CM = pd.merge(data,df_markers,on = "label",how = "left")
        data_CM_pred = pd.merge(data_CM,pred_alfas,on=["SampleID","label"],how="left")
        
        
        #print(data_CM)
        plt.imshow(z, origin='lower')#isns.imgplot(z)
        #cm = plt.cm.get_cmap("YlGn")#('RdYlBu')
        #cm = plt.cm.get_cmap("summer")#('RdYlBu') #fails on DC Keratin6
        #cm = plt.cm.get_cmap("viridis")#('RdYlBu') #fails on DC Keratin6
        #cm = plt.cm.get_cmap("bone")#('RdYlBu') #fails on DC Keratin6 becuase high Ker6 DCs diseapear in cancer niche
        #cm = plt.cm.get_cmap("pink")
        cm = matplotlib.colors.ListedColormap(["white","black"])
        cm2 = plt.cm.get_cmap("pink")#plt.cm.get_cmap("Greys")    
        gps_niches = data_CM_pred.groupby("inflammatory_niche_pred")
        alfas= list(gps_niches.groups)
        print(alfas)
        colors=["white","black"]
        cols = {alfas[i]:colors[i] for i in range(len(alfas))}
        plt.scatter(data_CM_pred['x'], data_CM_pred['y'], marker="o", s=50, c=data_CM_pred["value"],cmap=cm2,label="infl niche weight")
        #for (name, gp) in gps_niches:
        #    plt.scatter(gp['x'], gp['y'], marker="o", s=14, c = cols[name],label=name)
        plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
        #maxIntensity = data_CM[marker].quantile(1-intOutQuant)
        #data_CM.loc[data_CM[marker] > maxIntensity,marker ] = maxIntensity
        #plt.scatter(data_CM_pred['x'], data_CM_pred['y'], marker="o", s=14, c=data_CM_pred["inflammatory_niche_pred"],cmap=cm,label='inflammatory niche weight')

    plt.xlim(0, w)
    plt.ylim(0, h)
    plt.colorbar()
    if path_fig!=None:
        plt.savefig(path_fig,format="svg")
    
    
    return plt


def plot_cells_hclust_tmens(patient_id,cell_types,path_data, data_markers_path,cell_type,hclust_ids,clustIDs ,cols_clusts,segment_image=False,
                             segmentation_type="hard", counting_type="gaussian",h=800,w=800,granularity=20,radius=25,
                             pca_obj=None,AA_obj=None, to_plot=None,path_fig=None, intOutQuant=0):
    '''
    plots cells positions in MIBI image of TNBC overlayed with TMENs colors + saves it in a .svg image
    @param patient_id:{integer} id of the image patient (number here)
    @param cell_types: {list} of str, cell type labels in the dataset
    @param path_data: {string} path of data of patient (cells positions and cell types in an image)
    @param data_markers_path:{string} path to csv file containing markers expressions for each cells in all images
    @param cell_type:{string} cell type to display
    @param hclust_ids:{string} path to the .csv files with dataframe of cells and the cluster number assigned (from hclust) 
    @param clustIDs:{list} list of int of clusters IDs to visualize in colors
    @param cols_clusts:{list} list of colors for each cluster. Size = size(clustIDs)+1 ,the last color is the one for the
                                other cluster IDs
    @param segment_image: {boolean}, if True, color image from results of AA, if False, plot only cells (Default=False)
    @param segmentation_type:{str} ssegment pixels of images and color them as TMENs proportions (default='hard')
    @param counting_type:{str} counting type of cells within the sites:'abs' for absolute, 'log' fro log normalization and 'gaussian' for gaussian density (Default='abs')
    @param h:{int} height of image in pixels (Default=800)
    @param w:{int} width of image in pixels (Default=800)
    @param radius: {int} radius of sites generated in images in micrometer (Default=100)
    @param granularity: {int} granularity of colors (Default=25)
    @param pca_obj: {PCA obj} object of PCA on sites generated in the images (Default=None)
    @param AA_obj: {AA object} object of Archetype Analysis
    @param to_plot:{str} plot all images or not (Default='all')
    
    @return: plot
    '''    
    data = pd.read_csv(path_data+"/patient{}_cell_positions.csv".format(patient_id))
   # data = data.loc[data["cell_type"] == cell_type]
    groups = data.groupby('cell_type')
    df_markers = pd.read_csv(data_markers_path)
    plt.figure(figsize=(8, 8))
    if segment_image is True:
        if pca_obj is None or AA_obj is None:
                raise ValueError("To segment the image pca and archetypes objects are needed")
        color_vector =  np.array([[255, 0, 223],[255,0,0],[70,203,236],[0,0,0]])#np.array([[255, 0, 0], [0, 153, 51], [0, 0, 255], [255, 255, 0]])
        if segmentation_type == 'hard': #color pixel by 1 of the colors defining TMENs (discrete)
                color_fun = partial(alfa2color, color_vector)

        else: #color image by continuous spectrum of colors, depending on granularity 
            color_fun = partial(color_mapper, color_vector.T)

        z = get_segmentation_matrix(data, cell_types, pca_obj, AA_obj, color_fun, counting_type=counting_type, radius=radius, granularity=granularity, h=h, w=w)
        data = data.loc[data["cell_type"] == cell_type]
        df_markers = df_markers.loc[df_markers["SampleID"] == patient_id]
        df_markers.rename(columns = {"cellLabelInImage":"label"},inplace=True)
        data_CM = pd.merge(data,df_markers,on = "label",how = "left")
        # Get the hclust ids for each cell
        clusters = pd.read_csv(hclust_ids)
        data_CM_clusts = pd.merge(data_CM, clusters, on = ["SampleID","label"],how="left")
        #print(data_CM)
        plt.imshow(z, alpha=.4,origin='lower')
        #cm = plt.cm.get_cmap("pink")
        ### Change lcusts labels
        data_CM_clusts.loc[~data_CM_clusts["cluster"].isin(clustIDs),"cluster"]="others"
        
        hclustID = list(range(10))#list(set(data_CM_clusts["cluster"].tolist()))
        #colors = np.where(hclustID==6,"g",np.where(hclustID==3,"pink",np.where(hclustID==1,"brown","grey")))
        #labs = np.where(hclustID==6,6,np.where(hclustID==3,3,np.where(hclustID==1,1,"others")))
        
        #colors=["pink","white","green","grey"]#plt.get_cmap('Accent')
        print(len(colors))
        print(hclustID)
        #matplotlib.colors.ListedColormap(["pink","white","green","grey"])
        #sns.scatter('x', 'y',data=data_CM_clusts,hue="cluster",alpha=.7)
        groups = data_CM_clusts.groupby("cluster")
        #listNames = 
        
        hclustID = list(groups.groups)
        print(hclustID)
        cols = {hclustID[i]:cols_clusts[i] for i in range(len(hclustID))}
        
        for (name, gp) in groups:
            plt.scatter(gp['x'],gp['y'], marker="o", s=80,label=name,c=cols[name],alpha=.9)
        plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
    
    plt.xlim(0, w)
    plt.ylim(0, h)
    if path_fig!=None:
        plt.savefig(path_fig,format="svg")
    
    
    return plt


def plot_all_tumors_cell_positions(patient_ids, cell_types, segment_image=False, segmentation_type='hard', radius=100,
                                   granularity=25, pca_obj=None, AA_obj=None, to_plot='all',
                                   root_path='./data/cell_positions_data'):
    '''
    plots cells positions in MIBI image of tumors for all patients + saves in a .svg image
    :@param patient_ids: {list} of patient IDs (number for data from Keren et al,Cell 2018)
    :@param cell_types: {list} of str, cell type labels in the dataset
    :@param segment_image: {boolean}, if True, color image from results of AA, if False, plot only cells (Default=False)
    :@param segmentation_type:{str} (default='hard') boolean 
    :@param radius: {int} radius of sites generated in images (Default=100)
    :@param granularity: {int} granularity of colors (default=25)
    :@param pca_obj: {PCA object} object of PCA on sites generated in the images (Default=None)
    :@param AA_obj: {ArchetypeAnalysis object} object of Archetype Analysis
    :@param to_plot:{str} plot all images or not (Default='all')
    :@param root_path:{str} path of a folder containing csv data of cells (Default='./data/cell_positions_data')
    
    :@return None: 
    '''
    plt.figure(figsize=(30, 50))
    cells_cols = {cell_types[i]:colors[i] for i in range(len(cell_types))}#dict(zip(cell_types, colors)) 
    #key: cell type
    #Value: color
    #print(len(cells_cols.keys()))
    for i, patientID in enumerate(patient_ids):
        #print("Processing patient ID: {}".format(patientID))
        plt.subplot(8, 5, i+1)
        data = pd.read_csv("{}/patient{}_cell_positions.csv".format(root_path, patientID))
        groups = data.groupby('cell_type')
        plt.title("Patient ID: {}".format(patientID))
        if segment_image is True:
            if pca_obj is None or AA_obj is None:
                raise ValueError("To segment the image pca and archetypes objects are needed")
            color_vector = np.array([[255, 0, 223],[255,0,0],[70,203,236],[0,0,0]])#np.array([[255, 0, 0], [0, 153, 51], [0, 0, 255], [255, 255, 0]])
            if color_vector is None:
                colormap = mpl.cm.Dark2.colors
                
            if segmentation_type == 'hard': #color pixel by 1 of the colors defining TMENs (discrete)
                color_fun = partial(alfa2color, color_vector)
                
            else:#color image by continuous spectrum of colors, depending on granularity
                color_fun = partial(color_mapper, color_vector.T)
            try:
                h = data['size_y'][0]
                w = data['size_x'][0]
            except Exception: #default size of the image in pixels 
                h = 800
                w = 800
            z = get_segmentation_matrix(data, cell_types, pca_obj, AA_obj, color_fun, granularity=granularity, h=h, w=w,
                                        radius=radius)
            plt.imshow(z, origin='lower')
            
        for (name,group) in groups:
                    if to_plot == 'all' or name in to_plot:
                        #print(name)
                        plt.scatter(group['x'], group['y'], marker="o", s=5, label=name, c=cells_cols[name])
#        for (name, group), col in zip(groups, colors):
#            if to_plot == 'all' or name in to_plot:
#                plt.scatter(group['x'], group['y'], marker="o", s=5, label=name, c=col)
    plt.savefig("../../output/figs/all_patients_cells.svg",format="svg")
    plt.show()


def plot_scatter_pca(principal_components, evr, labels=None, original_axis=None, archetypes=None, cell_type=None):
    fig, ax = plt.subplots(figsize=(10,10))

    if labels is not None:
        ax.scatter(principal_components[:, 0], principal_components[:, 1],
                   c=labels, alpha=0.9, cmap=plt.cm.nipy_spectral)
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
            red_patch.append(mpatches.Patch(color=colors[i % len(colors)], label="{} = {}".format(i, cell_type[i])))
            plt.legend(handles=red_patch, title='Cell Type', bbox_to_anchor=(1.0, 1), loc='upper left')

    if archetypes is not None:
        plt.scatter(archetypes[0, :], archetypes[1, :], marker="^", s=100, color='orange')


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


def plot_3Dscatter_pca(principal_components, evr, labels, cmap=plt.cm.nipy_spectral, archetypes=None, halfspaces=None, original_axis = None, cell_type=None):
    tot_evr = np.sum(evr[0:3])
    print("{:.2f}% Total Exp. Var.".format(tot_evr))
    fig = plt.figure(1, figsize=(8, 8))
    plt.clf()
    ax = Axes3D(fig, rect=(0, 0, 1, 1), elev=21, azim=-58)

    minimum, maximum = _calculate_pca_max_min(principal_components[:, 0:3])
    ax.set_xlim(minimum, maximum)
    ax.set_ylim(minimum, maximum)
    ax.set_zlim(minimum, maximum)
    ax.set_xlabel("PC1 ({:.2f}% Exp. Var.)".format(evr[0] * 100))
    ax.set_ylabel("PC2 ({:.2f}% Exp. Var.)".format(evr[1] * 100))
    ax.set_zlabel("PC3 ({:.2f}% Exp. Var.)".format(evr[2] * 100))
    
    ax.xaxis.set_ticklabels([])#ax.set_xticks([])
    ax.yaxis.set_ticklabels([])#ax.set_yticks([])
    ax.zaxis.set_ticklabels([])#ax.set_zticks([])
    

    alpha = 0.1 if halfspaces is not None else 1
    alpha = 0 if original_axis is not None else 1
    ax.scatter(principal_components[:, 0], principal_components[:, 1], principal_components[:, 2],
               c=labels, cmap=cmap, edgecolor='k', alpha = alpha)

    if archetypes is not None:
        ax.scatter(archetypes[0, :], archetypes[1, :], archetypes[2, :],
                   marker = "^", s = 100, color='orange')

    if original_axis is not None:
        origin = original_axis[len(original_axis)-1, :]
        for i in range(len(original_axis)-1):
            plt.plot([origin[0], original_axis[i, 0]], [origin[1], original_axis[i, 1]], [origin[2], original_axis[i, 2]], 'k-', color='b', alpha=0.5)
            ax.text(original_axis[i, 0], original_axis[i, 1], original_axis[i, 2], cell_type[i], color='red')

    if halfspaces is not None:
        x, A, b = halfspaces
        verts = x.intersections
        hull = ConvexHull(verts)
        faces = hull.simplices

        for i, s in enumerate(faces):
            sq = [
                [verts[s[0], 0], verts[s[0], 1], verts[s[0], 2]],
                [verts[s[1], 0], verts[s[1], 1], verts[s[1], 2]],
                [verts[s[2], 0], verts[s[2], 1], verts[s[2], 2]]
            ]

            f = a3.art3d.Poly3DCollection([sq])

            eps=1e-6
            if np.linalg.norm(np.dot(A, np.array(sq).T) + b) < eps:
                f.set_color("#ff0000")
                f.set_alpha(1)
            else:
                f.set_color(matplotlib.colors.rgb2hex(np.random.rand(3)))
                f.set_alpha(0.1)
            f.set_edgecolor('k')
            ax.add_collection3d(f)

    return plt

def plot_3Dscatter_NNMF(nnmf_features, labels=None):
    fig = plt.figure(1, figsize=(8, 8))
    plt.clf()
    ax = Axes3D(fig, rect=(0, 0, 1, 1), elev=34, azim=-22)

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
    plt.figure(figsize=(3, 3))
    plt.plot(range(1, len(np.cumsum(evr))+1), np.cumsum(evr)*100)
    #v = np.argmax(np.cumsum(evr) > 0.95)
    plt.axvline(x=3, color='red', linestyle='--')
    plt.axhline(y=86, color='red', linestyle='--')
    plt.ylim(0, 100)
    plt.xlim(1, 10)
    plt.xlabel('# of principal components')
    plt.ylabel('% of cumulative variance explained')
    #plt.show()
    return plt


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


def archetypes_bar_plot(cell_number_archetypes, cell_types, colors, y_axis='count', radius=50,path_fig=None):
    if y_axis not in ['log', 'count', 'density']:
        raise ValueError("Wrong parameter: y_axis paramenter must be log, count or density")
    font = {'size'   : 14}
    matplotlib.rc('font', **font)
    fig = plt.figure(figsize=(10, 5))
    y_pos = 3*np.arange(len(cell_number_archetypes[0]))
    ax = fig.add_axes([0, 0, 1, 1])
    width = 0.50
    ## Get the colors
   
    
    data = cell_number_archetypes
    nbArch= data.shape[0]
    
    if y_axis == 'log':
        data = [np.log(d) for d in cell_number_archetypes]
    elif y_axis == 'density':
        data = [100 * (d / (radius*radius*np.pi)) for d in cell_number_archetypes]

    #print(y_pos)
    x = (len(cell_number_archetypes)//2) + len(cell_number_archetypes)%2
    #print(x)
    nbs = [4,2,1,3]
    #3-->1, 1-->4, 2-->2 , 4-->3
    #print(list(range(-x, x)))
    #print(data)
    if colors==None:
        colormap = mpl.cm.Dark2.colors 
        for d, idx, col in zip(data,list(range(-x, x)),colormap): #zip(data, colors, list(range(-x, x)),nbs):
            ax.bar(y_pos + idx * width, d, color=col, width=width, label="Arch "+str(idx+3))
    else:
        for d,c, idx, nb in zip(data, colors, list(range(-x, x)),nbs):
            ax.bar(y_pos + idx * width, data[nb-1], color=colors[nb-1], width=width, label="Arch "+str(nb))#ax.bar(y_pos + idx * width, d, color=c, width=width, label="Arch "+str(nb))
    

    #ax.bar(y_pos - 2*width, data[0], color=colors[0], width=width, label="Arch 1")
    #ax.bar(y_pos - width, data[1], color=colors[1], width=width, label="Arch 2")
    #ax.bar(y_pos, data[2], color=colors[2], width=width, label="Arch 3")
    #ax.bar(y_pos + width, data[3], color=colors[3], width=width, label="Arch 4")

    plt.xlabel('cell type')
    if y_axis == 'log':
        plt.ylabel('log(#cells)')
    elif y_axis == 'density':
        plt.ylabel('density [#cells per density area]')#'density [#cells / 100 um^2]'
    else:
        plt.ylabel('#cells')

    plt.legend()

    plt.xticks(y_pos, cell_types, rotation=90)
    #plt.savefig(path_fig,format="svg")
    plt.show()


def archetype_simple_plot(cell_number_archetypes, archetype_id, colors, cell_types, y_axis='count', radius=100):
    if y_axis not in ['log', 'count', 'density']:
        raise ValueError("Wrong parameter: y_axis paramenter must be log, count or density")

    if y_axis == 'log':
        cell_number_archetypes = np.log(cell_number_archetypes)
    elif y_axis == 'density':
        cell_number_archetypes = 100 * (cell_number_archetypes / (radius * radius * np.pi))

    plt.figure(figsize=(20, 5))
    y_pos = np.arange(len(cell_number_archetypes))
    plt.bar(y_pos, cell_number_archetypes, color=colors[archetype_id])
    plt.xticks(y_pos, cell_types, rotation=90)

    plt.xlabel('Cells Type')
    if y_axis == 'log':
        plt.ylabel('log(#cells)')
    elif y_axis == 'density':
        plt.ylabel('density [#cells / 100 um^2]')
    else:
        plt.ylabel('#cells')
    plt.title("Archetype {} {}".format(archetype_id+1, y_axis))


def get_explained_variance_matrix(X, Y, expl_var_ratio, cells_number=18):
    z = np.empty((X.shape[0], Y.shape[1]))
    for i, x in enumerate(X.T):
        print(x)
        if x[0] != 0:
            z[:, i] = np.insert(expl_var_ratio[x[0]], 0, 0.0)
        else:
            z[:, i] = np.zeros(cells_number)

    return z


def radius_pc_variance_contourf(patient_ids, expl_var_cum_ratio, cells_number=18):
    fig, axes = plt.subplots(nrows=8, ncols=5, figsize=(20, 40))
    axes = axes.flat
    for i, patientID in enumerate(patient_ids):
        x = [0] + list(expl_var_cum_ratio[patientID].keys())
        y = np.arange(0, cells_number)
        X, Y = np.meshgrid(x, y)
        Z = get_explained_variance_matrix(X, Y, expl_var_cum_ratio[patientID])
        im = axes[i].contourf(X, Y, Z)
        axes[i].axvline(x=100, color='r', linestyle='--')
        axes[i].axhline(y=3, color='r', linestyle='--')
        axes[i].set_ylim(0, 8)

        axes[i].set_title("Patient ID: {}".format(patientID))
        axes[i].set_xlabel("Radius")
        axes[i].set_ylabel("#PC")

        plt.colorbar(im, ax=axes[i])

    #plt.colorbar(im, ax=axes[0], orientation='horizontal')


    plt.show()


def radius_pc_heatmap(expl_var_ratio):
    radius = expl_var_ratio.keys()
    cumulative_var_exp_matrix = np.flip(np.vstack(list(expl_var_ratio.values())).T, 0)
    pcs_number = list(reversed(range(1, cumulative_var_exp_matrix.shape[0]+1)))

    fig, ax = plt.subplots(figsize=(12, 12))
    im = ax.imshow(cumulative_var_exp_matrix)

    ax.set_xticks(np.arange(len(radius)))
    ax.set_yticks(np.arange(len(pcs_number)))
    ax.set_xticklabels(radius)
    ax.set_yticklabels(pcs_number)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    plt.colorbar(im, fraction=0.03, pad=0.04)

    for i in range(len(pcs_number)):
        for j in range(len(radius)):
            text = ax.text(j, i, "{:.3f}".format(cumulative_var_exp_matrix[i, j]),
                           ha="center", va="center", color="w", fontsize=14)
            text.set_path_effects([path_effects.Stroke(linewidth=1.5, foreground='black'),
                                   path_effects.Normal()])

    #ax.set_title("Harvest of local farmers (in tons/year)")
    fig.tight_layout()
    plt.show()


def radius_pc_all_variance(expl_var_cum_ratio, radius_lim = 100,nPC_lim = 3,save_fig=False,path_fig=None,cells_number=18, pca_limit=8):
    plt.figure(figsize=(15, 15))
    x = [0] + list(expl_var_cum_ratio.keys())
    y = np.arange(0, cells_number)
    X, Y = np.meshgrid(x, y)
    expl_var_cum_ratio={int(i):k for i,k in expl_var_cum_ratio.items()}
    Z = get_explained_variance_matrix(X, Y, expl_var_cum_ratio, cells_number)
    im = plt.contourf(X, Y, Z)
    if radius_lim is not None and nPC_lim is not None:
        plt.axvline(x = radius_lim, color = 'r', linestyle = '--')  #x = 100
        plt.axhline(y = nPC_lim, color = 'r', linestyle = '--')  #y = 3
    plt.ylim(0, pca_limit)

    plt.xlabel("radius in micrometers")
    plt.ylabel("number of principal components")

    plt.colorbar()
    if save_fig==True:
        plt.savefig(path_fig)
    plt.show()


def plot_eigenvalues(poisson_ev, ev,path_fig=None):
    plt.figure(figsize=(10, 10))
    v = []
    if poisson_ev !=None:
        for pev in poisson_ev:
            #print(np.sum(pev),pev)
            plt.plot(list(range(1,18)),np.cumsum(pev)/np.sum(pev)*100, alpha=0.3)
    plt.plot(list(range(1,18)),np.cumsum(ev)/np.sum(ev)*100, color='black')
    plt.axvline(x=3,color="black",ls='--', lw=2)
    #print(ev)
    #plt.xticks(np.arange(1, len(ev), 1))
    plt.grid()
    plt.xlabel("Number of Principal Components (PCs)")
    plt.ylabel("% of variance explained")
    #plt.xlim(0, 18)
    if path_fig !=None:
        plt.savefig(path_fig,format="svg")
    plt.show()

