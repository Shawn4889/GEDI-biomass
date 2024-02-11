# Author: Xiaoxuan Li
import pandas as pd
import numpy as np
import shutil
import os
import h5py

# SA setup

studysite = ['Agincourt', 'DNyala', 'Welverdiendt',
             'Venetia', 'Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh']

beam0 = ["BEAM0000", "BEAM0001", "BEAM0010", "BEAM0011"]
beam1 = ["BEAM0101", "BEAM0110", "BEAM1000", "BEAM1011"]
GEDI_loc = r"E:\GEDI\GEDI_archive/"
GEDI_L4A_dir = GEDI_loc + "GEDI04A\SA_2023/"

CSV_dir = GEDI_loc + "SA.csv"
Result_dir = r"E:\GEDI\Result\Basic_GEDI/"

try:
    os.makedirs(Result_dir)
except:
    print("Folder already exists")

for Location in studysite:
    print(Location)
    # result directory
    coordinate_output = Result_dir + Location + "_GEDI_biomass.csv"
    coordinate_initial = []
    df_coordinate_initial = pd.DataFrame(coordinate_initial)
    df_geo = pd.read_csv(CSV_dir, sep=",")
    lat_max = df_geo.loc[df_geo['Location'] == Location, 'ul_lat'].iloc[0]
    lon_min = df_geo.loc[df_geo['Location'] == Location, 'ul_lon'].iloc[0]
    lat_min = df_geo.loc[df_geo['Location'] == Location, 'lr_lat'].iloc[0]
    lon_max = df_geo.loc[df_geo['Location'] == Location, 'lr_lon'].iloc[0]
    HDF_folder = os.listdir(GEDI_L4A_dir)
    for i in HDF_folder:
        # read HDF file
        file = GEDI_L4A_dir + i
        print(file)
        f_L4A = h5py.File(file, 'r')
        # beam list
        beam = list(f_L4A.keys())
        beam.remove("METADATA")
        beam.remove("ANCILLARY")
        for b in beam:
            print(b)
            # beam sel
            if b in beam0:
                Beam = b + "_0"
            else:
                Beam = b + "_1"
            print(Location + "_" + Beam)
            df_lat = pd.DataFrame(np.array(f_L4A[b + '/lat_lowestmode'][:]), columns=['lat'])
            df_long = pd.DataFrame(np.array(f_L4A[b + '/lon_lowestmode'][:]), columns=['long'])
            df_quality_l2A = pd.DataFrame(np.array(f_L4A[b + '/l2_quality_flag'][:]), columns=['ql2'])
            df_quality_l4A = pd.DataFrame(np.array(f_L4A[b + '/l4_quality_flag'][:]), columns=['ql4'])
            df_Coor = pd.concat([df_lat, df_long,df_quality_l2A,df_quality_l4A], axis=1)
            df_interaction = df_Coor[(df_Coor["lat"] > lat_min) &
                                     (df_Coor["lat"] < lat_max) &
                                     (df_Coor["long"] > lon_min) &
                                     (df_Coor["long"] < lon_max) &
                                     (df_Coor["ql2"] == 1) &
                                     (df_Coor["ql4"] == 1)]
            index_result = np.asarray(df_interaction.index.values)
            print(index_result)
            if not df_interaction.empty:
                # selected algor
                df_sa = pd.DataFrame(np.array(f_L4A[b + '/selected_algorithm'][:]), columns=['sel_al_biom'])
                # Geolocation & Quality data & AGBD
                df_shot = pd.DataFrame(np.array(f_L4A[b + '/shot_number'][:]), columns=['shot_number'])
                df_agbd = pd.DataFrame(np.array(f_L4A[b + '/agbd'][:]), columns=['AGBD'])
                df_agbd1 = pd.DataFrame(np.array(f_L4A[b + '/agbd_prediction/agbd_a1'][:]), columns=['AGBD_1'])
                df_agbd2 = pd.DataFrame(np.array(f_L4A[b + '/agbd_prediction/agbd_a2'][:]), columns=['AGBD_2'])
                df_agbd3 = pd.DataFrame(np.array(f_L4A[b + '/agbd_prediction/agbd_a3'][:]), columns=['AGBD_3'])
                df_agbd4 = pd.DataFrame(np.array(f_L4A[b + '/agbd_prediction/agbd_a4'][:]), columns=['AGBD_4'])
                df_agbd5 = pd.DataFrame(np.array(f_L4A[b + '/agbd_prediction/agbd_a5'][:]), columns=['AGBD_5'])
                df_agbd6 = pd.DataFrame(np.array(f_L4A[b + '/agbd_prediction/agbd_a6'][:]), columns=['AGBD_6'])
                df_agbd10 = pd.DataFrame(np.array(f_L4A[b + '/agbd_prediction/agbd_a10'][:]), columns=['AGBD_10'])
                # Merge all variables into one dataframe
                df = pd.concat([df_shot, df_agbd, df_sa, df_quality_l4A, df_agbd1,  df_agbd2,
                                df_agbd3, df_agbd4, df_agbd5, df_agbd6, df_agbd10], axis=1)
                df = df.loc[index_result, :]
                df["Beam"] = Beam
                df_coordinate_initial = pd.concat([df_coordinate_initial, df])
            else:
                print("no interaction between the GEDI track and bounding box")
    if not df_coordinate_initial.empty:
        df_coordinate_initial['shot_number'] = 'Shot_' + df_coordinate_initial['shot_number'].astype(str)
        df_coordinate_initial.to_csv(coordinate_output, index=None)
    else:
        print("No GEDI shots locate within the area: " + Location)
    if len(os.listdir(Result_dir)) == 0:
        shutil.rmtree(Result_dir)
