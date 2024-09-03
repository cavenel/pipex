import sys
import os
from skimage.io import imread
import numpy as np
import datetime
from tifffile import TiffFile, TiffWriter
import pandas as pd
import bigfish.detection as detection
from xml.etree import ElementTree
import anndata
import matplotlib.pyplot as plt

data_folder = os.environ.get('PIPEX_DATA')
spot_marker = ""
voxel_size = 103
spot_radius = 150
cluster_radius = 350
cluster_nb_min_spots = 4
dense_alpha = 0.7
dense_beta = 1
dense_gamma = 5

def options(argv):
    for arg in argv:
        if arg.startswith('-help'):
            print('Usage: \n\t-data=<optional /path/to/images/folder, defaults to /home/pipex/data> : example -> -data=/lab/projectX/images \n\t-spot_marker=<optional, name of the spot marker> : example -> -spot_marker=AMY2A \n\t-voxel_size=<optional, size of the voxel in nm, defaults to 103> : example -> -voxel_size=103 \n\t-spot_radius=<optional, radius of the spot in nm, defaults to 150> : example -> -spot_radius=150 \n\t-cluster_radius=<optional, radius of the cluster in nm, defaults to 350> : example -> -cluster_radius=350 \n\t-cluster_nb_min_spots=<optional, minimum number of spots in a cluster, defaults to 4> : example -> -cluster_nb_min_spots=4 \n\t-dense_alpha=<optional, alpha parameter for dense region decomposition, defaults to 0.7> : example -> -dense_alpha=0.7 \n\t-dense_beta=<optional, beta parameter for dense region decomposition, defaults to 1> : example -> -dense_beta=1 \n\t-dense_gamma=<optional, gamma parameter for dense region decomposition, defaults to 5> : example -> -dense_gamma=5', flush=True)
            sys.exit()
        elif arg.startswith('-data='):
            global data_folder
            data_folder = arg[6:]
        elif arg.startswith('-spot_marker='):
            global spot_marker
            spot_marker = arg[13:]
        elif arg.startswith('-voxel_size='):
            global voxel_size
            voxel_size = int(arg[12:])
        elif arg.startswith('-spot_radius='):
            global spot_radius
            spot_radius = int(arg[13:])
        elif arg.startswith('-cluster_radius='):
            global cluster_radius
            cluster_radius = int(arg[16:])
        elif arg.startswith('-cluster_nb_min_spots='):
            global cluster_nb_min_spots
            cluster_nb_min_spots = int(arg[33:])
        elif arg.startswith('-dense_alpha='):
            global dense_alpha
            dense_alpha = float(arg[13:])
        elif arg.startswith('-dense_beta='):
            global dense_beta
            dense_beta = float(arg[12:])
        elif arg.startswith('-dense_gamma='):
            global dense_gamma
            dense_gamma = int(arg[13:])

if __name__ == '__main__':
    options(sys.argv[1:])

    pidfile_filename = './RUNNING'
    if "PIPEX_WORK" in os.environ:
        pidfile_filename = './work/RUNNING'
    with open(pidfile_filename, 'w', encoding='utf-8') as f:
        f.write(str(os.getpid()))
        f.close()

    print(">>> Start time dot_extraction =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
    qptiff_files = [f for f in os.listdir(data_folder) if f.endswith(".qptiff")]
    if len(qptiff_files) > 0 or os.path.isfile(os.path.join(data_folder, spot_marker + ".tif")):
        if (not os.path.isfile(os.path.join(data_folder, spot_marker + ".tif"))):
            with TiffFile(os.path.join(data_folder, qptiff_files[0])) as tif:
                for page in tif.series[0].pages:
                    biomarker = ElementTree.fromstring(page.description).find('Biomarker').text
                    if biomarker == spot_marker:
                        if not os.path.isfile(os.path.join(data_folder, f'{biomarker}.tif')):
                            with TiffWriter(os.path.join(data_folder, f'{biomarker}.tif'), bigtiff=False) as tif:
                                tif.write(page.asarray())

        # read spot_marker image:
        img = imread(os.path.join(data_folder, spot_marker + ".tif"))
        img_name = spot_marker
        # Convert to int32
        img = img.astype(np.uint16)
        spots, threshold = detection.detect_spots(images=img, return_threshold=True, voxel_size=(voxel_size, voxel_size), spot_radius=(spot_radius, spot_radius))

        spots_post_decomposition, dense_regions, reference_spot = detection.decompose_dense(image=img, spots=spots, voxel_size=(voxel_size, voxel_size), spot_radius=(spot_radius, spot_radius), alpha=dense_alpha, beta=dense_beta, gamma=dense_gamma)

        spots_post_clustering, clusters = detection.detect_clusters(spots=spots_post_decomposition, voxel_size=(voxel_size, voxel_size), radius=cluster_radius, nb_min_spots=cluster_nb_min_spots)

        path = rf"{data_folder}/analysis/downstream/{img_name}_spots.csv"
        spots_df = pd.DataFrame(spots_post_clustering, columns=["y", "x", "cluster"])
        spots_df.to_csv(path, index=False)
        path = rf"{data_folder}/analysis/downstream/{img_name}_clusters.csv"
        clusters_df = pd.DataFrame(clusters, columns=["y", "x", "size", "cluster"])
        clusters_df.to_csv(path, index=False)

        # Add to adata object "F:\projects\SBPDA23\small\analysis\downstream\anndata_TissUUmaps.h5ad"
        adata = anndata.read_h5ad(rf"{data_folder}/analysis/downstream/anndata.h5ad")
        adata.uns["dots_spots"] = spots_df
        adata.uns["dots_clusters"] = clusters_df

        # Load segmentation data into numpy array format
        labels = np.load(os.path.join(data_folder, 'analysis', 'segmentation_data.npy'))
        df = pd.read_csv(os.path.join(data_folder, 'analysis', 'downstream', 'cell_data_norm.csv'))
        for col in ['dots_spots_count', 'dots_spots_density', 'dots_clusters_count', 'dots_clusters_sum']:
            if col in df.columns:
                df.drop(columns=[col], inplace=True)

        # Create a dataframe that will map each label to its count and size
        mapping_df = pd.DataFrame(df['cell_id'])
        mapping_df['dots_spots_count'] = 0
        mapping_df['dots_clusters_count'] = 0
        mapping_df['dots_clusters_sum'] = 0
        
        mapping_df.set_index('cell_id', inplace=True)

        # Update the count
        for i in spots_df.index:
            x = int(spots_df.loc[i, 'x'])
            y = int(spots_df.loc[i, 'y'])
            label = labels[y, x]
            if label in mapping_df.index:
                mapping_df.loc[label, 'dots_spots_count'] += 1
        for i in clusters_df.index:
            x = int(clusters_df.loc[i, 'x'])
            y = int(clusters_df.loc[i, 'y'])
            size = int(clusters_df.loc[i, 'size'])
            label = labels[y, x]
            if label in mapping_df.index:
                mapping_df.loc[label, 'dots_clusters_count'] += 1
                mapping_df.loc[label, 'dots_clusters_sum'] += size

        # Merge the count and size data back to the main dataframe
        df = df.merge(mapping_df, how='inner', left_on='cell_id', right_index=True)
        
        # Calculate dots density
        df['dots_spots_density'] = df['dots_spots_count'] / df['size']

        # Save df to file
        df.to_csv(os.path.join(data_folder, 'analysis', 'downstream', 'cell_data_norm.csv'), index=False)

        adata.obs['dots_spots_count'] = df['dots_spots_count'].values.astype(int)
        adata.obs['dots_spots_density'] = df['dots_spots_density'].values.astype(float)
        adata.obs['dots_clusters_count'] = df['dots_clusters_count'].values.astype(int)
        adata.obs['dots_clusters_sum'] = df['dots_clusters_sum'].values.astype(int)

        adata.write_h5ad(rf"{data_folder}/analysis/downstream/anndata.h5ad")
    
    print(">>> End time dot_extraction =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
