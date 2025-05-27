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

data_folder = os.environ.get('PIPEX_DATA')
spot_markers = []
voxel_size = 103
spot_radius = 150
cluster_radius = 350
cluster_nb_min_spots = 4
dense_alpha = 0.7
dense_beta = 1
dense_gamma = 5

def options(argv):
    global data_folder, spot_markers, voxel_size, spot_radius
    global cluster_radius, cluster_nb_min_spots, dense_alpha, dense_beta, dense_gamma

    for arg in argv:
        if arg.startswith('-help'):
            print('Usage: \n\t-data=<optional /path/to/images/folder, defaults to /home/pipex/data> : example -> -data=/lab/projectX/images \n\t-spot_markers=<optional, list of spot markers> : example -> -spot_markers=AMY2A,SST,GORASP2 \n\t-voxel_size=<optional, size of the voxel in nm, defaults to 103> : example -> -voxel_size=103 \n\t-spot_radius=<optional, radius of the spot in nm, defaults to 150> : example -> -spot_radius=150 \n\t-cluster_radius=<optional, radius of the cluster in nm, defaults to 350> : example -> -cluster_radius=350 \n\t-cluster_nb_min_spots=<optional, minimum number of spots in a cluster, defaults to 4> : example -> -cluster_nb_min_spots=4 \n\t-dense_alpha=<optional, alpha parameter for dense region decomposition, defaults to 0.7> : example -> -dense_alpha=0.7 \n\t-dense_beta=<optional, beta parameter for dense region decomposition, defaults to 1> : example -> -dense_beta=1 \n\t-dense_gamma=<optional, gamma parameter for dense region decomposition, defaults to 5> : example -> -dense_gamma=5', flush=True)
            sys.exit()
        elif arg.startswith('-data='):
            data_folder = arg[6:]
        elif arg.startswith('-spot_markers='):
            spot_markers = arg[14:].split(',')
        elif arg.startswith('-voxel_size='):
            voxel_size = int(arg[12:])
        elif arg.startswith('-spot_radius='):
            spot_radius = int(arg[13:])
        elif arg.startswith('-cluster_radius='):
            cluster_radius = int(arg[16:])
        elif arg.startswith('-cluster_nb_min_spots='):
            cluster_nb_min_spots = int(arg[22:])
        elif arg.startswith('-dense_alpha='):
            dense_alpha = float(arg[13:])
        elif arg.startswith('-dense_beta='):
            dense_beta = float(arg[12:])
        elif arg.startswith('-dense_gamma='):
            dense_gamma = int(arg[13:])

if __name__ == '__main__':
    options(sys.argv[1:])

    pidfile_filename = './RUNNING'
    if "PIPEX_WORK" in os.environ:
        pidfile_filename = './work/RUNNING'
    with open(pidfile_filename, 'w', encoding='utf-8') as f:
        f.write(str(os.getpid()))

    print(">>> Start time dot_extraction =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)

    qptiff_files = [f for f in os.listdir(data_folder) if f.endswith(".qptiff")]
    adata = anndata.read_h5ad(rf"{data_folder}/analysis/downstream/anndata.h5ad")
    labels = np.load(os.path.join(data_folder, 'analysis', 'segmentation_data.npy'))
    df = pd.read_csv(os.path.join(data_folder, 'analysis', 'downstream', 'cell_data_norm.csv'))

    for marker in spot_markers:
        print(f"Processing marker: {marker}", flush=True)
        marker_file = os.path.join(data_folder, f"{marker}.tif")

        if not os.path.isfile(marker_file) and len(qptiff_files) > 0:
            with TiffFile(os.path.join(data_folder, qptiff_files[0])) as tif:
                for page in tif.series[0].pages:
                    biomarker = ElementTree.fromstring(page.description).find('Biomarker').text
                    if biomarker == marker:
                        with TiffWriter(marker_file, bigtiff=False) as out_tif:
                            out_tif.write(page.asarray())

        if not os.path.isfile(marker_file):
            print(f"Marker {marker} not found, skipping.")
            continue

        img = imread(marker_file).astype(np.uint16)
        spots, _ = detection.detect_spots(images=img, return_threshold=True, voxel_size=(voxel_size, voxel_size), spot_radius=(spot_radius, spot_radius))
        try:
            spots_post_decomposition, _, _ = detection.decompose_dense(image=img, spots=spots, voxel_size=(voxel_size, voxel_size), spot_radius=(spot_radius, spot_radius), alpha=dense_alpha, beta=dense_beta, gamma=dense_gamma)
        except Exception as e:
            print(f"Error during decomposition for marker {marker}: {e}", flush=True)
            spots_post_decomposition = spots    
        spots_post_clustering, clusters = detection.detect_clusters(spots=spots_post_decomposition, voxel_size=(voxel_size, voxel_size), radius=cluster_radius, nb_min_spots=cluster_nb_min_spots)

        spots_df = pd.DataFrame(spots_post_clustering, columns=["y", "x", "cluster"])
        clusters_df = pd.DataFrame(clusters, columns=["y", "x", "size", "cluster"])
        spots_df.to_csv(os.path.join(data_folder, "analysis", "downstream", f"{marker}_spots.csv"), index=False)
        clusters_df.to_csv(os.path.join(data_folder, "analysis", "downstream", f"{marker}_clusters.csv"), index=False)

        adata.uns[f"{marker}_dots_spots"] = spots_df
        adata.uns[f"{marker}_dots_clusters"] = clusters_df

        mapping_df = pd.DataFrame({'cell_id': df['cell_id']})
        mapping_df.set_index('cell_id', inplace=True)
        mapping_df[f"{marker}_dots_spots_count"] = 0
        mapping_df[f"{marker}_dots_clusters_count"] = 0
        mapping_df[f"{marker}_dots_clusters_sum"] = 0

        for _, row in spots_df.iterrows():
            x, y = int(row['x']), int(row['y'])
            label = labels[y, x]
            if label in mapping_df.index:
                mapping_df.loc[label, f"{marker}_dots_spots_count"] += 1
        for _, row in clusters_df.iterrows():
            x, y, size = int(row['x']), int(row['y']), int(row['size'])
            label = labels[y, x]
            if label in mapping_df.index:
                mapping_df.loc[label, f"{marker}_dots_clusters_count"] += 1
                mapping_df.loc[label, f"{marker}_dots_clusters_sum"] += size

        df = df.drop(columns=[col for col in df.columns if col.startswith(f"{marker}_dots_")], errors='ignore')
        df = df.merge(mapping_df, how='inner', left_on='cell_id', right_index=True)
        df[f"{marker}_dots_spots_density"] = df[f"{marker}_dots_spots_count"] / df['size']

        adata.obs[f"{marker}_dots_spots_count"] = df[f"{marker}_dots_spots_count"].values.astype(int)
        adata.obs[f"{marker}_dots_spots_density"] = df[f"{marker}_dots_spots_density"].values.astype(float)
        adata.obs[f"{marker}_dots_clusters_count"] = df[f"{marker}_dots_clusters_count"].values.astype(int)
        adata.obs[f"{marker}_dots_clusters_sum"] = df[f"{marker}_dots_clusters_sum"].values.astype(int)

    df.to_csv(os.path.join(data_folder, 'analysis', 'downstream', 'cell_data_norm.csv'), index=False)
    adata.write_h5ad(rf"{data_folder}/analysis/downstream/anndata.h5ad")

    print(">>> End time dot_extraction =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)
