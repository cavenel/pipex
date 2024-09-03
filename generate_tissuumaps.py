import scanpy as sc
import scipy
import os
import json
import datetime
import sys
from skimage.measure import approximate_polygon
import numpy as np
import struct
from tifffile import TiffFile, TiffWriter
from xml.etree import ElementTree

data_folder = os.environ.get('PIPEX_DATA')
include_marker_images = "no"
include_geojson = "no"
compress_geojson = "no"
include_html = "no"

def exporting_tissuumaps ():
    # Check if the required files are present
    global include_geojson
    if include_marker_images == "no" and include_geojson == "yes":
        print(">>> Impossible to display geojson without a background image", flush=True)
        include_geojson = "no"
    if include_geojson == "yes":
        if not os.path.exists(os.path.join(data_folder, 'analysis/cell_segmentation_geo.json')):
            print(">>> Impossible to display geojson without a cell segmentation file", flush=True)
            include_geojson = "no"
    
    adata = sc.read_h5ad(os.path.join(data_folder, 'analysis/downstream/anndata.h5ad'))

    # Make sure that the X matrix is in the compressed sparse column (CSC) format (required by TissUUmaps)
    adata.X = scipy.sparse.csc_matrix(adata.X)

    # Add image layers and cell segmentation geoJSON file to the AnnData object:
    if include_marker_images == "yes":
        markers = adata.var_names
    elif include_marker_images == "no":
        markers = []
    else:
        markers = include_marker_images.split(",")
    
    qptiff_files = [f for f in os.listdir(data_folder) if f.endswith(".qptiff")]
    if len(qptiff_files) > 0:
        with TiffFile(os.path.join(data_folder, qptiff_files[0])) as tif:
            for page in tif.series[0].pages:
                biomarker = ElementTree.fromstring(page.description).find('Biomarker').text
                print (biomarker, markers)
                if biomarker in markers:
                    with TiffWriter(os.path.join(data_folder, f'{biomarker}.tif'), bigtiff=False) as tif:
                        tif.write(page.asarray())
                    #with open(os.path.join(data_folder, f'{biomarker}.tif'), 'wb') as f:
                    #    f.write(page.asarray())

    if include_geojson == "yes" and include_marker_images != "no":
        if compress_geojson == "yes":
            import geobuf
            # load os.path.join(data_folder, 'analysis/cell_segmentation_geo.json') as json
            with open(os.path.join(data_folder, 'analysis/cell_segmentation_geo.json'), 'r') as f:
                geojson_data = json.load(f)
                # Remove measurements from geojson_data[i].properties.measurements
                for f in geojson_data:
                    if "measurements" in f["properties"]:
                        del f["properties"]["measurements"]
                # Approximate region boundaries to reduce file size:
                geojson_data_approx = []
                for f in geojson_data:
                    f_approx = f.copy()
                    for i in range(len(f["geometry"]["coordinates"])):
                            f_approx["geometry"]["coordinates"][i] = approximate_polygon(np.array(f["geometry"]["coordinates"][i]), tolerance=0.75).tolist()
                    geojson_data_approx.append(f_approx)
                geojson_data_approx = {
                    "type":"FeatureCollection",
                    "features":geojson_data_approx
                }
                # Save as a compressed protobuf file in the geobuf format:
                pbf = geobuf.encode(geojson_data_approx, 3) # TissUUmaps uses a precision of 3 decimal
                # We need to add a tag to indicate the precision
                data = geobuf.geobuf_pb2.Data()
                data.ParseFromString(pbf)
                data.precision = 3
                pbf = data.SerializeToString()

                cell_segmentation_path = "../cell_segmentation_geo.pbf"
                with open(os.path.join(data_folder, 'analysis/cell_segmentation_geo.pbf'), 'wb') as f:
                    f.write(pbf)
        else:
            cell_segmentation_path = "../cell_segmentation_geo.json"
        regionFiles = [
            {
                "name": "Cell segmentation",
                "path": cell_segmentation_path,
                "title": "Cell segmentation",
                "settings":[
                    {
                        "module": "regionUtils",
                        "function": "_regionStrokeWidth",
                        "value": "0.5"
                    },
                    {
                        "module": "regionUtils",
                        "function": "_regionStrokeAdaptOnZoom",
                        "value": True
                    },
                    {
                        "module": "glUtils",
                        "function": "_regionShowOnTop",
                        "value": False
                    }
                ]
            }
        ]
    else:
        regionFiles = []
    #check if adata.uns["PLA_spots"] exists
    if "PLA_spots" in adata.uns:
        markerFiles = [
            {
                "title": "Detected PLA spots",
                "name": "Detected PLA spots",
                "autoLoad": False,
                "hideSettings": True,
                "uid": "RawSpots",
                "expectedHeader": {
                    "X": "/uns/PLA_spots/x",
                    "Y": "/uns/PLA_spots/y",
                    "gb_col": "",
                    "gb_name": "",
                    "cb_cmap": "",
                    "cb_col": "",
                    "cb_gr_dict": "{\"All\":\"#ffff00\"}",
                    "scale_col": "",
                    "scale_factor": 0.1,
                    "coord_factor": 1,
                    "pie_col": "",
                    "pie_dict": "",
                    "shape_col": "",
                    "shape_fixed": "disc",
                    "shape_gr_dict": "",
                    "edges_col": "/obsp/spatial_connectivities;join",
                    "collectionItem_col": "",
                    "collectionItem_fixed": "0",
                    "opacity_col": "",
                    "opacity": 1,
                    "stroke_width": 5.5,
                    "sortby_col": "",
                    "z_order": 1,
                    "tooltip_fmt": ""
                },
                "expectedRadios": {
                    "cb_col": False,
                    "cb_gr": True,
                    "cb_gr_rand": False,
                    "cb_gr_dict": True,
                    "cb_gr_key": False,
                    "pie_check": False,
                    "scale_check": False,
                    "shape_col": False,
                    "shape_gr": False,
                    "shape_gr_rand": True,
                    "shape_gr_dict": False,
                    "shape_fixed": True,
                    "opacity_check": False,
                    "_no_outline": False,
                    "no_fill": True,
                    "collectionItem_col": False,
                    "collectionItem_fixed": True,
                    "sortby_check": False,
                    "sortby_desc_check": False,
                    "edges_check": False
                },
                "path": "data/files/anndata_TissUUmaps.h5ad",
                "fromButton": 3
                }
        ]
    else:
        markerFiles = []
    adata.uns["tmap"] = json.dumps({
        "markerFiles": markerFiles, 
        "layers": [
            {
                "name": f"{marker}",
                "tileSource": f"../../{marker}.tif.dzi"
            }
            for marker in markers
        ],
        "regionFiles": regionFiles,
        "plugins": ["Feature_Space","InteractionQC","Spot_Inspector"],
        "settings": [
            {
                "module": "pluginUtils",
                "function": "startPlugin",
                "value": ["Spot_Inspector",
                [
                    {"name": "_layer_format", "value":"{layout-row6}"},
                    {"name": "_cmap", "value":"undefined"},
                ],False]
            },
            {
                "module": "pluginUtils",
                "function": "startPlugin",
                "value": ["InteractionQC",[],False]
            },
            {
                "module": "pluginUtils",
                "function": "startPlugin",
                "value": ["Feature_Space",[],False]
            }
        ],
    })
    adata.write_h5ad(os.path.join(data_folder, 'analysis/downstream/anndata_TissUUmaps.h5ad'))

    if include_html == "yes":
        import tissuumaps

        state = tissuumaps.read_h5ad.h5ad_to_tmap("", os.path.join(data_folder, 'analysis/downstream/anndata_TissUUmaps.h5ad'))
        tissuumaps.views.exportToStatic(
            json.dumps(state), 
            os.path.join(data_folder, 'analysis/downstream/TissUUmaps_webexport/'),
            os.path.join(data_folder, 'analysis/downstream')
        )
        from urllib.request import urlretrieve
        for plugin in ["Feature_Space","InteractionQC","Spot_Inspector"]:
            url = f"https://tissuumaps.github.io/TissUUmaps/plugins/latest/{plugin}.js"
            filename = os.path.join(data_folder, 'analysis/downstream/TissUUmaps_webexport/plugins/', f"{plugin}.js")
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            urlretrieve(url, filename)

#Function to handle the command line parameters passed
def options(argv):
    if (len(argv) == 0):
        print('export_tissuumaps.py arguments:\n\t-data=<optional /path/to/images/folder, defaults to /home/pipex/data> : example -> -data=/lab/projectX/images\n\t-include_marker_images=<yes or no or list of present specific markers to display as image layers> : example -> -include_marker_images=DAPI,SST,GORASP2\n\t-include_geojson=<yes or no to include cell segmentation as regions> : example -> -include_geojson=yes\n\t-compress_geojson=<yes or no to compress geojson regions into pbf> : example -> -compress_geojson=yes\n\t-include_html=<yes or no to export html page for sharing the TissUUmaps project on the web> : example -> -include_marker_images=yes', flush=True)
        sys.exit()
    else:
        for arg in argv:
            if arg.startswith('-help'):
                print('export_tissuumaps.py arguments:\n\t-data=<optional /path/to/images/folder, defaults to /home/pipex/data> : example -> -data=/lab/projectX/images\n\t-include_marker_images=<yes or no or list of present specific markers to display as image layers> : example -> -include_marker_images=DAPI,SST,GORASP2\n\t-include_geojson=<yes or no to include cell segmentation as regions> : example -> -include_geojson=yes\n\t-compress_geojson=<yes or no to compress geojson regions into pbf> : example -> -compress_geojson=yes\n\t-include_html=<yes or no to export html page for sharing the TissUUmaps project on the web> : example -> -include_marker_images=yes', flush=True)
                sys.exit()
            elif arg.startswith('-data='):
                global data_folder
                data_folder = arg[6:]
            elif arg.startswith('-include_marker_images='):
                global include_marker_images
                include_marker_images = arg[23:]
            elif arg.startswith('-include_geojson='):
                global include_geojson
                include_geojson = arg[17:]
            elif arg.startswith('-compress_geojson='):
                global compress_geojson
                compress_geojson = arg[18:]
            elif arg.startswith('-include_html='):
                global include_html
                include_html = arg[14:]

if __name__ =='__main__':
    options(sys.argv[1:])

    pidfile_filename = './RUNNING'
    if "PIPEX_WORK" in os.environ:
        pidfile_filename = './work/RUNNING'
    with open(pidfile_filename, 'w', encoding='utf-8') as f:
        f.write(str(os.getpid()))
        f.close()
    with open(data_folder + '/log_settings_tissuumaps.txt', 'w+', encoding='utf-8') as f:
        f.write(">>> Start time tissuumaps = " + datetime.datetime.now().strftime(" %H:%M:%S_%d/%m/%Y") + "\n")
        f.write(' '.join(sys.argv))
        f.close()

    print(">>> Start time exporting tissuumaps =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)

    exporting_tissuumaps()

    print(">>> End time exporting tissuumaps =", datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"), flush=True)