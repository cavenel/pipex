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


def options(argv):
    for arg in argv:
        if arg.startswith('-help'):
            print('Usage: \n\t-data=<optional /path/to/images/folder, defaults to /home/pipex/data> : example -> -data=/lab/projectX/images', flush=True)
            sys.exit()
        elif arg.startswith('-data='):
            global data_folder
            data_folder = arg[6:]


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
    if len(qptiff_files) > 0 or os.path.isfile(os.path.join(data_folder, "DO.tif")):
        if (not os.path.isfile(os.path.join(data_folder, "DO.tif"))):
            with TiffFile(os.path.join(data_folder, qptiff_files[0])) as tif:
                for page in tif.series[0].pages:
                    biomarker = ElementTree.fromstring(page.description).find('Biomarker').text
                    if biomarker == "DO":
                        if not os.path.isfile(os.path.join(data_folder, f'{biomarker}.tif')):
                            with TiffWriter(os.path.join(data_folder, f'{biomarker}.tif'), bigtiff=False) as tif:
                                tif.write(page.asarray())

        # read DO image:
        img = imread(os.path.join(data_folder, "DO.tif"))
        img_name = "DO"
        # Convert to int32
        img = img.astype(np.uint16)
        spots, threshold = detection.detect_spots(images=img, return_threshold=True, voxel_size=(103, 103), spot_radius=(150, 150))

        spots_post_decomposition, dense_regions, reference_spot = detection.decompose_dense(image=img, spots=spots, voxel_size=(103, 103), spot_radius=(150, 150), alpha=0.7, beta=1, gamma=5)

        spots_post_clustering, clusters = detection.detect_clusters(spots=spots_post_decomposition, voxel_size=(103, 103), radius=350, nb_min_spots=4)

        path = rf"{data_folder}/analysis/downstream/{img_name}_spots.csv"
        spots_df = pd.DataFrame(spots_post_clustering, columns=["y", "x", "cluster"])
        spots_df.to_csv(path, index=False)
        path = rf"{data_folder}/analysis/downstream/{img_name}_clusters.csv"
        clusters_df = pd.DataFrame(clusters, columns=["y", "x", "size", "cluster"])
        clusters_df.to_csv(path, index=False)

        # Add to adata object "F:\projects\SBPDA23\small\analysis\downstream\anndata_TissUUmaps.h5ad"
        adata = anndata.read_h5ad(rf"{data_folder}/analysis/downstream/anndata.h5ad")
        adata.uns["PLA_spots"] = spots_df
        adata.uns["PLA_clusters"] = clusters_df

        # Load segmentation data into numpy array format
        labels = np.load(os.path.join(data_folder, 'analysis', 'segmentation_data.npy'))
        df = pd.read_csv(os.path.join(data_folder, 'analysis', 'cell_data.csv'))
        try:
            df.drop(columns=['PLA_count', 'PLA_density'], inplace=True)
        except:
            pass
        # Create a dataframe that will map each label to its count and size
        mapping_df = pd.DataFrame(df['cell_id'])
        mapping_df['PLA_count'] = 0
        mapping_df.set_index('cell_id', inplace=True)

        # Update the count
        for i in spots_df.index:
            x = int(spots_df.loc[i, 'x'])
            y = int(spots_df.loc[i, 'y'])
            label = labels[y, x]
            if label in mapping_df.index:
                mapping_df.loc[label, 'PLA_count'] += 1

        # Merge the count and size data back to the main dataframe
        df = df.merge(mapping_df, how='inner', left_on='cell_id', right_index=True)
        
        # Calculate PLA density
        df['PLA_density'] = df['PLA_count'] / df['size']

        # Save df to file
        df.to_csv(os.path.join(data_folder, 'analysis', 'cell_data.csv'), index=False)
        
        adata.obs['PLA_count'] = df['PLA_count'].values.astype(int)
        adata.obs['PLA_density'] = df['PLA_density'].values.astype(float)

        adata.write_h5ad(rf"{data_folder}/analysis/downstream/anndata.h5ad")
    
    adata = anndata.read_h5ad(rf"{data_folder}/analysis/downstream/anndata.h5ad")

    # Compute the distribution of the intensity of the markers
    
    values_PD1 = adata.X[:, np.where(adata.var.index == "PD-1")].flatten().copy()
    values_PDL1 = adata.X[:, np.where(adata.var.index == "PD-L1")].flatten().copy()
    values_minPD1PDL1 = np.minimum(values_PD1, values_PDL1)
    values_timesPD1PDL1 = values_PD1 * values_PDL1

    markers = {
        "PD-1": values_PD1,
        "PD-L1": values_PDL1,
        "Min of PD-1 and PD-L1": values_minPD1PDL1,
        "PD-1 x PD-L1": values_timesPD1PDL1
    }
    fig, ax = plt.subplots(4, 2, figsize=(20, 20))
    for marker_index, (marker_name, marker) in enumerate(markers.items()):
        #Compute distribution of marker weigthed by PLA_count values:
        
        # Get the PLA count values
        counts = adata.obs["PLA_count"].to_numpy().copy()
        # Compute the weighted distribution
        weighted_values = np.repeat(marker, counts)
        # Compute the distribution
        distribution = np.histogram(weighted_values, bins=20, range=(0, 20), density=True)
        # Plot the distribution as bars
        ax[marker_index, 0].bar(distribution[1][:-1], distribution[0], width=distribution[1][1]-distribution[1][0])
        ax[marker_index, 0].set_xlim(-1, 20)
        ax[marker_index, 0].set_xlabel("Weighted intensity")
        ax[marker_index, 0].set_ylabel("Density")
        ax[marker_index, 0].set_title(f"Weighted distribution of {marker_name} intensity")
        
        # Plot distribution of marker for all cells:

        distribution = np.histogram(marker, bins=20, range=(0, 20), density=True)

        ax[marker_index, 1].bar(distribution[1][:-1], distribution[0], width=distribution[1][1]-distribution[1][0])
        ax[marker_index, 1].set_xlim(-1, 20)
        ax[marker_index, 1].set_xlabel("Intensity")
        ax[marker_index, 1].set_ylabel("Density")
        ax[marker_index, 1].set_title(f"Distribution of {marker_name} intensity")
    plt.tight_layout()
    plt.savefig(os.path.join(data_folder, 'analysis', 'downstream', 'intensity_distribution.png'))

    print(">>> End time dot_extraction =", datetime.datetime.now().strftime("%H:%M:%S"), flush=True)
