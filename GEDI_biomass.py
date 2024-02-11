# Author: Xiaoxuan Li
import os
import arcpy
from arcpy.sa import *
import pandas as pd
import numpy as np
arcpy.env.overwriteOutput = True



def merge_biomass_GEDI():
    Result_dir = r"E:\GEDI\Result\Basic_GEDI/"
    out_csv = r"E:\GEDI\Result\Result\All_biomass.csv"
    initial = []
    df = pd.DataFrame(initial)
    csv_folder = os.listdir(Result_dir)
    for i in csv_folder:
        if "biomass" in i:
            print("Merging..." + i)
            in_csv = Result_dir + i
            df_item = pd.read_csv(in_csv)
            df = pd.concat([df, df_item])
    df.to_csv(out_csv, index=None)
    rh_csv = r"E:\GEDI\Result\Result\All_PCM.csv"
    df_rh = pd.read_csv(rh_csv)
    df_merge = pd.merge(df_rh, df, left_on='shot_number', right_on='shot_number', how='left')
    out_final = r"E:\GEDI\Result\Result\All_rh_biom.csv"
    print("Joining...")
    df_merge.to_csv(out_final, index=None)


def biom_gedi_ras_output():
    studysite = ['DNyala', 'Venetia','Justicia', 'Ireagh',
                 'Limpopo1', 'Limpopo2', 'Limpopo3']
    for site in studysite:
        print(site)
        shp = r"E:\GEDI\Result\SHP/" + site + "_c_buffer.shp"
        ras_chm = r"E:\ALS_archive\LiDAR_CHM_Mosaic/" + site + ".tif"
        out_dir = r"E:\Biomass\CSIR\ras/" + site + "/"
        arcpy.MakeFeatureLayer_management(shp, "fLayer")
        print([f.name for f in arcpy.ListFields("fLayer")])
        with arcpy.da.SearchCursor("fLayer", ["Unnamed__1",'SHAPE@']) as cursor:
            for row in cursor:
                print(row[0])
                query = """ "Unnamed__1" = {0}""".format(row[0])
                arcpy.SelectLayerByAttribute_management("fLayer", 'NEW_SELECTION', query)
                outraster = out_dir + "ras_" + str(row[0]) + ".tif"
                extent = str(row[1].extent.XMin) + " " + \
                         str(row[1].extent.YMin) + " " + \
                         str(row[1].extent.XMax) + " " + \
                         str(row[1].extent.YMax)
                print("Output..." + outraster)
                arcpy.Clip_management(ras_chm, extent, outraster,
                                      "fLayer", 0, "ClippingGeometry", "NO_MAINTAIN_EXTENT")
        arcpy.Delete_management('fLayer')


#biom_gedi_ras_output()


def int_zonal_biom():
    studysite = ['Agincourt','DNyala', 'Venetia','Justicia', 'Ireagh',
                 'Limpopo1', 'Limpopo2', 'Limpopo3']
    for site in studysite:
        print(site)
        arcpy.env.workspace = r"E:\Biomass\CSIR\ras/" + site
        ras_filenames = arcpy.ListRasters('*.tif*')
        list = []
        list_file = []
        for ras_filename in ras_filenames:
            print(ras_filename)
            if "ras_" in ras_filename:
                    list_file.append(ras_filename)
                    ras = r"E:\Biomass\CSIR\ras/" + site + "/" + ras_filename
                    zon_ras = r"E:\Biomass\CSIR\lidR_crown/" + site + "/" + ras_filename
                    arc_dir = r"D:\temp\ArcGIS_project\GEDI\GEDI.gdb/"
                    out = r"E:\Biomass\CSIR\result/" + site + "/"
                    table = arc_dir + "table_height_" + ras_filename.split(".")[0]
                    outInt1 = Int(zon_ras)
                    ZonalStatisticsAsTable(outInt1, "Value", ras, table, statistics_type="MAXIMUM")
                    xlsx = out + ras_filename.split("_")[1] + ".xlsx"
                    arcpy.TableToExcel_conversion(table, xlsx)
                    df = pd.read_excel(xlsx)
                    df = df[df.COUNT > 8]
                    df['biomass'] = pow(10, (0.0207 +
                                             0.301 * np.log10(df['AREA']) +
                                             0.196 * pow(np.log10(df['AREA']), 2) +
                                             1.74 * np.log10(df['MAX']))) / 1000
                    print(sum(df['biomass']))
                    list.append(sum(df['biomass']))

        df = pd.DataFrame(list)
        df_list = pd.DataFrame(list_file)
        df = pd.concat([df,df_list],axis=1)
        df.to_csv(r'E:\Biomass\CSIR\result/' + site + '/biom.csv')


def concat_GEDI_ALS_biom():
    #studysite = ['Agincourt','DNyala', 'Venetia','Justicia', 'Ireagh',
    #             'Welverdiendt','Limpopo1', 'Limpopo2', 'Limpopo3']
    studysite = ["Welverdiendt"]
    for site in studysite:
        print(site)

        GEDI_biom = r"E:\GEDI\Result\Result/"+ site + "_c.csv"
        out = r"E:\Biomass\CSIR\result/" + site + "_biom_combo.csv"
        ALS_tree_biom = r"E:\Biomass\CSIR\result/" + site + "/biom.csv"
        ALS_shurb_biom = r"E:\Biomass\CSIR\result/" + site + "/biom_shurb.csv"
        df_GEDI_biom = pd.read_csv(GEDI_biom)

        df_ALS_tree_biom = pd.read_csv(ALS_tree_biom, skiprows=1, usecols=range(1,3), names=["Tree", "Name1"])
        df_ALS_shurb_biom = pd.read_csv(ALS_shurb_biom, skiprows=range(0,2), usecols=range(1,3), names=["Shurb", "Name2"])
        df_GEDI_biom.rename(columns={'Unnamed: 0.1': 'Name'}, inplace=True)

        df_GEDI_biom['Name'] = "ras_" + df_GEDI_biom['Name'].astype(str) + ".tif"
        df_ALS_shurb_biom['Name2'] = df_ALS_shurb_biom['Name2'] + ".tif"
        df = pd.concat([df_ALS_tree_biom,df_ALS_shurb_biom],axis=1)
        df_merge = pd.merge(df, df_GEDI_biom, left_on='Name1', right_on='Name', how='left')
        df_merge.drop_duplicates(keep='first', inplace=True)
        df_merge['AGB'] =  df_merge['AGBD'] * 0.049
        df_merge['Combo'] = df_merge['Tree'] + df_merge['Shurb']
        df_merge['Diff'] = df_merge['Combo'] - df_merge['AGB']
        df_merge.to_csv(out,index=None)


def int_zonal_biom_CSIR():
    arcpy.env.workspace = r"E:\Biomass\CSIR\ras\CSIR/"
    ras_filenames = arcpy.ListRasters('*_h.tif*')
    list = []
    list_file = []
    for ras_filename in ras_filenames:
        try:
            print(ras_filename)
            list_file.append(ras_filename)
            ras = r"E:\Biomass\CSIR\ras\CSIR/" + ras_filename
            zon_ras = r"E:\Biomass\CSIR\lidR_crown\CSIR/" + ras_filename
            arc_dir = r"D:\temp\ArcGIS_project\GEDI\GEDI.gdb/"
            out = r"E:\Biomass\CSIR\result\CSIR/"
            table = arc_dir + "table_height_" + ras_filename.split(".")[0]
            outInt1 = Int(zon_ras)
            ZonalStatisticsAsTable(outInt1, "Value", ras, table, statistics_type="MAXIMUM")
            xlsx = out + ras_filename.split("_")[1] + ".xlsx"
            arcpy.TableToExcel_conversion(table, xlsx)
            df = pd.read_excel(xlsx)
            df = df[df.COUNT > 8]
            df['biomass'] = pow(10, (0.0207 +
                                     0.301 * np.log10(df['AREA']) +
                                     0.196 * pow(np.log10(df['AREA']), 2) +
                                     1.74 * np.log10(df['MAX']))) / 1000
            print(sum(df['biomass']))
            list.append(sum(df['biomass']))
        except:
            print("No file found..." + ras_filename)
    df = pd.DataFrame(list)
    df_list = pd.DataFrame(list_file)
    df = pd.concat([df,df_list],axis=1)
    df.to_csv(r'E:\Biomass\CSIR\result\CSIR\biom_tree_CSIR.csv')


def height_clip_25m():
    dir_shp = r"E:\Biomass\CSIR\shp/CSIR_25m_Ven65.shp"
    dir_out = r"E:\Biomass\CSIR\ras\CSIR_25m/"
    dir_ras_h = r"E:\ALS_archive\LiDAR_CHM_Mosaic/"
    arcpy.MakeFeatureLayer_management(dir_shp, 'CSIR_25')
    Proj_dir = r'E:\GEDI\Boundingbox/'
    with arcpy.da.UpdateCursor('CSIR_25', ['FID']) as cursor:
        for row in cursor:
            print(row[0])
            query = """ "FID" = %s""" % row[0]
            arcpy.SelectLayerByAttribute_management('CSIR_25', 'NEW_SELECTION', query)
            site = "Venetia"
            #Proj = Proj_dir + "UTM35S.prj"
            out_shp_p = r"E:\Biomass\CSIR_5_p_25.shp"
            arcpy.Project_management('CSIR_25', out_shp_p, Proj)
            arcpy.MakeFeatureLayer_management(out_shp_p, "CSIR_proj_25")
            ras_chm = dir_ras_h + site + ".tif"
            outraster_h = dir_out + site + "_Ven_" + str(row[0]) + "_h.tif"
            desc = arcpy.Describe("CSIR_proj_25")
            extent = str(desc.extent.XMin) + " " + \
                     str(desc.extent.YMin) + " " + \
                     str(desc.extent.XMax) + " " + \
                     str(desc.extent.YMax)
            print("Output..." + outraster_h)
            arcpy.Clip_management(ras_chm, extent, outraster_h, "CSIR_proj_25", -999, "ClippingGeometry", "NO_MAINTAIN_EXTENT")



def int_zonal_biom_CSIR_25m():
    arcpy.env.workspace = r"E:\Biomass\CSIR\ras\CSIR_25m/"
    ras_filenames = arcpy.ListRasters('*_Ven_*')
    list = []
    list_file = []
    for ras_filename in ras_filenames:
        print(ras_filename)
        try:
            list_file.append(ras_filename)
            ras = r"E:\Biomass\CSIR\ras\CSIR_25m/" + ras_filename
            zon_ras = r"E:\Biomass\CSIR\lidR_crown\CSIR_25m/" + ras_filename
            arc_dir = r"D:\temp\ArcGIS_project\GEDI\GEDI.gdb/"
            out = r"E:\Biomass\CSIR\result\CSIR_25m/"
            table = arc_dir + "table_height_" + ras_filename.split(".")[0]
            outInt1 = Int(zon_ras)
            ZonalStatisticsAsTable(outInt1, "Value", ras, table, statistics_type="MAXIMUM")
            xlsx = out + ras_filename.split("_")[1] + ".xlsx"
            arcpy.TableToExcel_conversion(table, xlsx)
            df = pd.read_excel(xlsx)
            df = df[df.COUNT > 8]
            df['biomass'] = pow(10, (0.0207 +
                                     0.301 * np.log10(df['AREA']) +
                                     0.196 * pow(np.log10(df['AREA']), 2) +
                                     1.74 * np.log10(df['MAX']))) / 1000
            print(sum(df['biomass']))
            list.append(sum(df['biomass']))
        except:
            print("No file found..." + ras_filename)
    df = pd.DataFrame(list)
    df_list = pd.DataFrame(list_file)
    df = pd.concat([df,df_list],axis=1)
    df.to_csv(r'E:\Biomass\CSIR\result\CSIR_25m\biom_CSIR_25m_Ven.csv')



def concat_GEDI_ALS_biom_CSIR():
    GEDI_biom = r"E:\GEDI\Result\Result\All.csv"
    ALS_tree_biom = r"E:\Biomass\CSIR\result\biom_CSIR.csv"
    ALS_shurb_biom = r"E:\Biomass\CSIR\result\biom_shurb_CSIR.csv"
    df_GEDI_biom = pd.read_csv(GEDI_biom, usecols=['ql2a', 'ql2b', 'ql4',
                                                  'site', 'orbit', 'Beam',
                                                  'Unnamed: 0.1', 'AGBD','status'])

    df_ALS_tree_biom = pd.read_csv(ALS_tree_biom, skiprows=1, usecols=range(1,3), names=["Tree", "Name1"])
    df_ALS_shurb_biom = pd.read_csv(ALS_shurb_biom, skiprows=range(0,2), usecols=range(1,3), names=["Shurb", "Name2"])
    df_GEDI_biom.rename(columns={'Unnamed: 0.1': 'Name'}, inplace=True)

    df_GEDI_biom['Name'] = "ras_" + df_GEDI_biom['Name'].astype(str) + ".tif"
    df_ALS_shurb_biom['Name2'] = df_ALS_shurb_biom['Name2'] + ".tif"
    print(df_ALS_tree_biom)
    print(df_ALS_shurb_biom)
    print(df_GEDI_biom)
    df = pd.concat([df_ALS_tree_biom,df_ALS_shurb_biom],axis=1)
    df_merge = pd.merge(df, df_GEDI_biom, left_on='Name1', right_on='Name', how='left')
    df_merge.drop_duplicates(keep='first', inplace=True)
    df_merge['AGB'] =  df_merge['AGBD'] * 0.049
    df_merge['Combo'] = df_merge['Tree'] + df_merge['Shurb']
    df_merge['Diff'] = df_merge['Combo'] - df_merge['AGB']
    df_merge.to_csv(r'E:\Biomass\CSIR\result\biom_CSIR_combo.csv',index=None)



def concat_ALS_biom_CSIR():
    ALS_tree_biom = r"E:\Biomass\CSIR\result\CSIR\biom_tree_CSIR.csv"
    ALS_shurb_biom = r"E:\Biomass\CSIR\result\CSIR\biom_shurb_CSIR.csv"

    df_ALS_tree_biom = pd.read_csv(ALS_tree_biom, skiprows=1, usecols=range(1,3), names=["Tree", "Name1"])
    df_ALS_shurb_biom = pd.read_csv(ALS_shurb_biom, skiprows=range(0,2), usecols=range(1,3), names=["Shurb", "Name2"])

    df_ALS_shurb_biom['Name2'] = df_ALS_shurb_biom['Name2'] + ".tif"
    print(df_ALS_tree_biom)
    print(df_ALS_shurb_biom)

    df_merge = pd.concat([df_ALS_tree_biom,df_ALS_shurb_biom],axis=1)
    df_merge['Combo'] = df_merge['Tree'] + df_merge['Shurb']
    df_merge.to_csv(r'E:\Biomass\CSIR\result\CSIR\biom_CSIR_combo.csv',index=None)


def concat_ALS_biom_CSIR_25m():
    ALS_tree_biom = r"E:\Biomass\CSIR\result\CSIR_25m\biom_CSIR_25m.csv"
    ALS_shurb_biom = r"E:\Biomass\CSIR\result\CSIR_25m\biom_shurb_CSIR_25m.csv"

    df_ALS_tree_biom = pd.read_csv(ALS_tree_biom, skiprows=1, usecols=range(1,3), names=["Tree", "Name1"])
    df_ALS_shurb_biom = pd.read_csv(ALS_shurb_biom, skiprows=range(0,2), usecols=range(1,3), names=["Shurb", "Name2"])

    df_ALS_shurb_biom['Name2'] = df_ALS_shurb_biom['Name2'] + ".tif"
    print(df_ALS_tree_biom)
    print(df_ALS_shurb_biom)

    df_merge = pd.concat([df_ALS_tree_biom,df_ALS_shurb_biom],axis=1)
    df_merge['Combo'] = df_merge['Tree'] + df_merge['Shurb']
    df_merge.to_csv(r'E:\Biomass\CSIR\result\CSIR_25m\biom_CSIR_25m_combo.csv',index=None)


def merge_biom_pos():
    studysite = ['Agincourt','DNyala', 'Venetia','Justicia', 'Ireagh',
                 'Welverdiendt','Limpopo1', 'Limpopo2', 'Limpopo3']
    for site in studysite:
        print(site)
        dir_biom = r"E:\Biomass\CSIR\individual_tree\Biom_clean.csv"
        dir_pos = r"E:\Biomass\CSIR\individual_tree\Justicia_merge.xlsx"
        out = r"E:\Biomass\CSIR\individual_tree\Justicia.csv"
        df_biom = pd.read_csv(dir_biom)
        df_pos = pd.read_excel(dir_pos)
        print(df_biom)
        print(df_pos)
        df_merge = pd.merge(df_pos, df_biom, left_on='combo', right_on='combo', how='left')
        df_merge.to_csv(out,index=None)


def merge_biom_result():
    initial = []
    initial_df = pd.DataFrame(initial)
    dir = r"E:\Biomass\CSIR\result/"
    biom_csvs = os.listdir(dir)
    for file in biom_csvs:
        if "csv" in file:
            file_df = pd.read_csv(dir + file)
            initial_df = pd.concat([initial_df,file_df])
    initial_df.to_csv(r"E:\Biomass\CSIR\result\All.csv")


def GEDI_L2B_merge():
    initial=[]
    df_initial = pd.DataFrame(initial)
    dir = r"E:\GEDI\Result\Basic_GEDI/"
    out_dir = r"E:\GEDI\Result\Basic_GEDI\All_L2B.csv"
    lists = os.listdir(dir)
    for file in lists:
        if "_L2B.csv" in file:
            df = pd.read_csv(dir+file)
            df_initial = pd.concat([df_initial,df])
    df_initial.to_csv(out_dir)


def merge_biom_l2b():
    dir1 = r"E:\Biomass\CSIR\result\All.csv"
    dir2 = r"E:\GEDI\Result\Basic_GEDI\All_L2B.csv"
    out = r"E:\Biomass\CSIR\result\All_biom_l2b.csv"
    df_biom = pd.read_csv(dir1)
    df_l2b = pd.read_csv(dir2)
    df_merge = pd.merge(df_biom, df_l2b, left_on='shot_number', right_on='shot_number', how='left')
    df_merge.drop_duplicates(keep='first', inplace=True)
    df_merge.to_csv(out,index=None)


def change_combo():
    dir1 = r"E:\Biomass\CSIR\result\All_biom_l2b.csv"
    alldf = pd.read_csv(dir1)
    df = alldf.copy()
    for i in range(0, len(df)):
        if df.iloc[i]['Tree'] < df.iloc[i]['Shurb']:
            df.at[i, 'Combo'] = df.iloc[i]['Tree'] + df.iloc[i]['Shurb']
        else:
            df.at[i, 'Combo'] = df.iloc[i]['Tree']
    out = r"E:\Biomass\CSIR\result\All_biom_l2b_2.csv"
    df.to_csv(out, index=None)


def biom_csv_shp():
    Result_dir= r"E:\Biomass\CSIR\result/"
    merge_files = os.listdir(Result_dir)
    Proj_dir = r'E:\GEDI\Boundingbox/'
    GDB_dir = r"D:\temp\ArcGIS_project\GEDI\GEDI.gdb/"
    SHP_dir = r"E:\Biomass\CSIR\shp/"
    for file in merge_files:
        if "_biom_combo" in file:
            if "Addo" in file or "Venetia" in file or "DNyala" in file:
                Proj = Proj_dir + "UTM35S.prj"
            elif "Agulhas" in file:
                Proj = Proj_dir + "UTM34S.prj"
            else:
                Proj = Proj_dir + "UTM36S.prj"
            out_shp = GDB_dir + str(file).split(".")[0]
            out_buffer = SHP_dir + str(file).split(".")[0]
            print("Converting csv to buffered shp: " + out_buffer)
            arcpy.management.XYTableToPoint(Result_dir + file, out_shp,
                                            "longitude", "latitude", '', Proj)
            arcpy.Buffer_analysis(out_shp, out_buffer, "12.5 Meters")


def zonal_max_mean():
    shp_dir = r"E:\Biomass\CSIR\shp/"
    chm_dir = r"E:\ALS_archive\LiDAR_CHM_Mosaic/"
    GDB_dir = r"D:\temp\ArcGIS_project\GEDI\GEDI.gdb/"
    zonal_dir = r"E:\Biomass\CSIR\result/"
    arcpy.env.workspace = shp_dir
    fcs = arcpy.ListFeatureClasses()
    for fc in fcs:
        print("Zonal stats... " + fc)
        table = GDB_dir + "table_" + fc.split(".")[0]
        chm = chm_dir + fc.split(".")[0].split("_")[0] + ".tif"
        ZonalStatisticsAsTable(fc, "shot_numbe", chm, table, '')
        xlsx = zonal_dir + fc.split(".")[0] + "_zonal.xlsx"
        arcpy.TableToExcel_conversion(table, xlsx)
        print("Join..." + table)
        csv1 = zonal_dir + fc.split(".")[0] + ".csv"
        out = zonal_dir + fc.split(".")[0] + "_z.csv"
        df_csv1 = pd.read_csv(csv1)
        df_csv2 = pd.read_excel(xlsx)
        df_csv2 = df_csv2[["shot_numbe", "MAX", "MEAN"]]
        df_merge = pd.merge(df_csv1, df_csv2, left_on='shot_number', right_on='shot_numbe', how='left')
        df_merge.drop_duplicates(keep='first', inplace=True)
        df_merge.to_csv(out, index=None)


def merge_biomass_Z():
    Result_dir = r"E:\Biomass\CSIR\result/"
    out_csv = Result_dir + r"All_z.csv"
    initial = []
    df = pd.DataFrame(initial)
    csv_folder = os.listdir(Result_dir)
    for i in csv_folder:
        if "_biom_combo_z.csv" in i:
            print("Merging..." + i)
            in_csv = Result_dir + i
            df_item = pd.read_csv(in_csv)
            df = pd.concat([df, df_item])
    df.to_csv(out_csv, index=None)


def merge_biom_z_l2b():
    dir1 = r"E:\Biomass\CSIR\result\All_z.csv"
    dir2 = r"E:\GEDI\Result\Basic_GEDI\All_L2B.csv"
    out = r"E:\Biomass\CSIR\result\All_biom_z_L2B.csv"
    df_biom = pd.read_csv(dir1)
    df_l2b = pd.read_csv(dir2)
    df_merge = pd.merge(df_biom, df_l2b, left_on='shot_number', right_on='shot_number', how='left')
    df_merge.drop_duplicates(keep='first', inplace=True)
    df_merge.to_csv(out,index=None)


def zonal_max_p98():
    shp_dir = r"E:\Biomass\CSIR\shp/"
    chm_dir = r"E:\ALS_archive\LiDAR_CHM_Mosaic/"
    GDB_dir = r"D:\temp\ArcGIS_project\GEDI\GEDI.gdb/"
    zonal_dir = r"E:\Biomass\CSIR\result/"
    arcpy.env.workspace = shp_dir
    fcs = arcpy.ListFeatureClasses()
    for fc in fcs:
        print("Zonal stats... " + fc)
        table = GDB_dir + "P98"
        chm = chm_dir + fc.split(".")[0].split("_")[0] + ".tif"
        ZonalStatisticsAsTable(fc, "shot_numbe", chm, table, '', 'PERCENTILE', '', "98")
        xlsx = zonal_dir + fc.split(".")[0] + "_p98.xlsx"
        arcpy.TableToExcel_conversion(table, xlsx)
        print("Join..." + table)
        csv1 = zonal_dir + fc.split(".")[0] + ".csv"
        out = zonal_dir + fc.split(".")[0] + "_p98.csv"
        df_csv1 = pd.read_csv(csv1)
        df_csv2 = pd.read_excel(xlsx)
        df_csv2 = df_csv2[["shot_numbe", "PCT98"]]
        df_merge = pd.merge(df_csv1, df_csv2, left_on='shot_number', right_on='shot_numbe', how='left')
        df_merge.drop_duplicates(keep='first', inplace=True)
        df_merge.to_csv(out, index=None)


def merge_biomass_p98():
    Result_dir = r"E:\Biomass\CSIR\result/"
    out_csv = Result_dir + r"All_p98.csv"
    initial = []
    df = pd.DataFrame(initial)
    csv_folder = os.listdir(Result_dir)
    for i in csv_folder:
        if "_biom_combo_p98.csv" in i:
            print("Merging..." + i)
            in_csv = Result_dir + i
            df_item = pd.read_csv(in_csv)
            df = pd.concat([df, df_item])
    df.to_csv(out_csv, index=None)


def join_p98():
    dir1 = r"E:\GEDI\Result\backup\Result\All.csv"
    dir2 = r"E:\Biomass\CSIR\result\All_p98.csv"
    out = r"E:\GEDI\Result\backup\Result\All_p98.csv"
    df_ori = pd.read_csv(dir1)
    df_p98 = pd.read_csv(dir2)
    df_p98 = df_p98[["shot_number", "PCT98"]]
    df_merge = pd.merge(df_ori,df_p98, left_on='shot_number', right_on='shot_number', how='left')
    df_merge.drop_duplicates(keep='first', inplace=True)
    df_merge.to_csv(out,index=None)


def remove_dup():
    dir = r"E:\Biomass\CSIR\individual_tree\Biomass_trees.csv"
    out = r"E:\Biomass\CSIR\individual_tree\Biomass_trees_2.csv"
    df = pd.read_csv(dir)
    df = df.sort_values('Biomass', ascending=False).drop_duplicates(['ROI', 'Site', 'SP', 'GPS ']).sort_index()
    df.to_csv(out,index=None)


def extract_ras_pnt():
    inPointFeatures = r"D:\temp\ArcGIS_project\Biomass\Biomass.gdb\Ireagh_1"
    inRaster = r"E:\Biomass\CSIR\AGBD_CHM_CC_2012\2012_L7_veght_Grt0-5m.tif"
    outPointFeatures = r"D:\temp\ArcGIS_project\Biomass\Biomass.gdb\Ireagh_1_18_12"
    ExtractValuesToPoints(inPointFeatures, inRaster, outPointFeatures, "", "VALUE_ONLY")
    print("Complete...Ireagh_1_18_12")
    inPointFeatures = r"D:\temp\ArcGIS_project\Biomass\Biomass.gdb\Justicia_1"
    inRaster = r"E:\Biomass\CSIR\AGBD_CHM_CC_2012\2012_L4L5_Veght_Grt0-5m.tif"
    outPointFeatures = r"D:\temp\ArcGIS_project\Biomass\Biomass.gdb\Justicia_1_18_12"
    ExtractValuesToPoints(inPointFeatures, inRaster, outPointFeatures, "", "VALUE_ONLY")
    print("Complete...Justicia_1_18_12")
    inPointFeatures = r"D:\temp\ArcGIS_project\Biomass\Biomass.gdb\Welverdiendt_1"
    inRaster = r"E:\Biomass\CSIR\AGBD_CHM_CC_2012\2012_andoverwelverd_veght_Grt0-5m.tif"
    outPointFeatures = r"D:\temp\ArcGIS_project\Biomass\Biomass.gdb\Welverdiendt_1_18_12"
    ExtractValuesToPoints(inPointFeatures, inRaster, outPointFeatures, "", "VALUE_ONLY")
    print("Complete...Welverdiendt_1_18_12")


def output_1m_table():
    inTable = r"D:\temp\ArcGIS_project\Biomass\Biomass.gdb\Ireagh_1_18_12"
    outTable = r"E:\Biomass\CSIR\ALS_1\Ireagh.csv"
    expression = "RASTERVALU IS NOT NULL And grid_code > 0.5"
    arcpy.conversion.ExportTable(inTable, outTable, expression, "NOT_USE_ALIAS")
    print("Complete...Ireagh")
    inTable = r"D:\temp\ArcGIS_project\Biomass\Biomass.gdb\Justicia_1_18_12"
    outTable = r"E:\Biomass\CSIR\ALS_1\Justicia.csv"
    expression = "RASTERVALU IS NOT NULL And grid_code > 0.5"
    arcpy.conversion.ExportTable(inTable, outTable, expression, "NOT_USE_ALIAS")
    print("Complete...Justicia")
    inTable = r"D:\temp\ArcGIS_project\Biomass\Biomass.gdb\Welverdiendt_1_18_12"
    outTable = r"E:\Biomass\CSIR\ALS_1\Welverdiendt.csv"
    expression = "RASTERVALU IS NOT NULL And grid_code > 0.5"
    arcpy.conversion.ExportTable(inTable, outTable, expression, "NOT_USE_ALIAS")
    print("Complete...Welverdiendt")


def zonal_2012_2018():
    arc_dir = r"D:\temp\ArcGIS_project\Biomass\Biomass.gdb/"
    chm_AGBD_2012 = Raster(r"E:\Biomass\CSIR\AGBD_CHM_CC_2012\2012_VegVars_Stack_AGB_CC_H_25m.tif\Band_1")
    chm_H_2012 = Raster(r"E:\Biomass\CSIR\AGBD_CHM_CC_2012\2012_VegVars_Stack_AGB_CC_H_25m.tif\Band_3")
    studysite = ['Ireagh']
    for site in studysite:
        inpolygon = r"E:\Biomass\CSIR\Result_0606\shp_25/" + site + "_25_Sel.shp"
        #2012 AGBD
        print("Input ras: " + site)
        table = arc_dir + site + "_2012_AGBD"
        print("Output zonal stats..." + table)
        ZonalStatisticsAsTable(inpolygon, "FID", chm_AGBD_2012, table, "", "MEAN")
        #2012 H
        table = arc_dir + site + "_2012_H"
        print("Output zonal stats..." + table)
        ZonalStatisticsAsTable(inpolygon, "FID", chm_H_2012, table, "", "MEAN")
        #2018 H
        chm_H_2018 = Raster(r"E:\ALS_archive\LiDAR_CHM_Mosaic/"+site+".tif")
        table = arc_dir + site + "_2018_H"
        print("Output zonal stats..." + table)
        ZonalStatisticsAsTable(inpolygon, "FID", chm_H_2018, table, "", "MEAN")


def zonal_2012_2018_table_output():
    arc_dir = r"D:\temp\ArcGIS_project\Biomass\Biomass.gdb/"
    output_dir = r"E:\Biomass\CSIR\Result_0606\zonal_table/"
    studysite = ['Agincourt', 'Welverdiendt', 'Justicia', 'Ireagh']
    for site in studysite:
        #index
        inTable = "E:\Biomass\CSIR\Result_0606\shp_25/" + site + "_25_Sel.shp"
        outTable = output_dir + site + "_index.csv"
        print("Output zonal stats..." + outTable)
        arcpy.conversion.ExportTable(inTable, outTable, "", "NOT_USE_ALIAS")
        #2018 H
        inTable = arc_dir + site + "_2018_H"
        outTable = output_dir + site + "_2018_H.csv"
        print("Output zonal stats..." + outTable)
        arcpy.conversion.ExportTable(inTable, outTable, "", "NOT_USE_ALIAS")

        #2012 H
        inTable = arc_dir + site + "_2012_H"
        outTable = output_dir + site + "_2012_H.csv"
        print("Output zonal stats..." + outTable)
        arcpy.conversion.ExportTable(inTable, outTable, "", "NOT_USE_ALIAS")
        #2012 AGBD
        inTable = arc_dir + site + "_2012_AGBD"
        outTable = output_dir + site + "_2012_AGBD.csv"
        print("Output zonal stats..." + outTable)
        arcpy.conversion.ExportTable(inTable, outTable, "", "NOT_USE_ALIAS")





def AGB_point():
    #PALSAR2_ScanSAR_HV_mtf_5_db_2014-09-06.tif
    #PALSAR2_ScanSAR_HV_mtf_5_db_2014-11-29.tif
    #PALSAR2_ScanSAR_HV_mtf_5_db_2015-01-10.tif
    #PALSAR2_ScanSAR_HV_mtf_5_db_2015-02-21.tif
    #PALSAR2_ScanSAR_HV_mtf_5_db_2015-04-04.tif
    suf = "PALSAR2_ScanSAR_HV_mtf_5_db_2015-04-04.tif"
    dir = r"E:\Biomass\CSIR\AGBD_CHM_CC_2012\Resample_25/"
    shp = r"E:\Biomass\CSIR\Result_0614\shp\2015-04-04/"
    sar = r"E:\ScanSAR\ScanSAR\single/" + suf
    out = r"E:\Biomass\CSIR\Result_0614\SAR_Biomass_14_15\2015-04-04/"
    studysite = ['Agincourt', 'Welverdiendt', 'Justicia', 'Ireagh', 'L1', 'L2', 'Athol']
    for site in studysite:
        dir_agb = dir + site + "_biomass.tif"
        dir_shp = shp + site + ".shp"
        out_shp = shp + site + "_SAR.shp"
        outTable = out + site + ".csv"
        print(site)
        print("Raster to point...")
        arcpy.RasterToPoint_conversion(dir_agb, dir_shp)
        print("Extract SAR to ALS AGB shp...")
        ExtractValuesToPoints(dir_shp, sar, out_shp, "", "VALUE_ONLY")
        print("Shp to table...")
        arcpy.conversion.ExportTable(out_shp, outTable, "", "NOT_USE_ALIAS")


def AGB_point_12_18():
    arc_dir = r"D:\temp\ArcGIS_project\Biomass_new\Biomass_new.gdb/"
    output = r"E:\Biomass\CSIR\Result_0615/"
    studysite = ['Agincourt', 'Welverdiendt', 'Justicia', 'Ireagh']
    for site in studysite:
        years = ['2012', '2018']
        for year in years:
            print(site)
            print(year)
            inFeatures = arc_dir + site + "_c_25_p"
            agb_ras = arc_dir + site + "_" + str(year) + "_SAS"
            outFeatureClass = arc_dir + site + "_c_25_p_" + str(year)
            outTable = output + site + "_" + year + "_AGB.csv"
            ExtractValuesToPoints(inFeatures, agb_ras, outFeatureClass, "", "VALUE_ONLY")
            arcpy.conversion.ExportTable(outFeatureClass, outTable, "", "NOT_USE_ALIAS")



def SAR_point_14_15():
    arc_dir = r"D:\temp\ArcGIS_project\Biomass_new\Biomass_new.gdb/"
    output = r"E:\Biomass\CSIR\Result_0615/"
    studysite = ['Agincourt', 'Welverdiendt', 'Justicia', 'Ireagh']
    for site in studysite:
        sars = ['PALSAR2_ScanSAR_HV_mtf_5_db_2014-09-06', 'PALSAR2_ScanSAR_HV_mtf_5_db_2014-11-29',
                'PALSAR2_ScanSAR_HV_mtf_5_db_2015-01-10','PALSAR2_ScanSAR_HV_mtf_5_db_2015-02-21',
                'PALSAR2_ScanSAR_HV_mtf_5_db_2015-04-04']
        for sar in sars:
            print(site)
            inFeatures = arc_dir + site + "_c_25_p"
            suf = sar.split('_')[6].split('-')[0] + sar.split('_')[6].split('-')[1] + "_SAS"
            sar_ras = r"E:\ScanSAR\ScanSAR\single/" + sar + ".tif"
            outFeatureClass = arc_dir + site + "_c_25_p_SAR_" + suf
            outTable = output + site + "_SAR_" + suf + ".csv"
            ExtractValuesToPoints(inFeatures, sar_ras, outFeatureClass, "", "VALUE_ONLY")
            arcpy.conversion.ExportTable(outFeatureClass, outTable, "", "NOT_USE_ALIAS")




def zonalstats_CHM_75m():
    studysite = ['Agincourt', 'Welverdiendt', 'Justicia', 'Ireagh']
    for site in studysite:
        arc_dir = r"D:\temp\ArcGIS_project\ScanSAR\ScanSAR.gdb/"
        inpolygon = r"E:\ScanSAR\ScanSAR\Grid/" + site + "_75m.shp"
        chm = r"E:\ALS_archive\LiDAR_CHM_Mosaic/" + site + ".tif"
        print("Input ras: " + chm)
        table = arc_dir + "CHM_2012_" + site + "_75m"
        print("Output zonal stats..." + table)
        ZonalStatisticsAsTable(inpolygon, "FID", chm, table, "", "MEAN")
        xlsx = r"E:\ScanSAR\ScanSAR\site_75m/" + site + "/" + "CHM_2012_" + site + "_75m" + ".xlsx"
        arcpy.TableToExcel_conversion(table, xlsx)


def biom_csv_ind_0619():
    csv_dir = r"E:\Biomass\CSIR\Result_0618\2018/csv/"
    Proj_dir = r"E:\GEDI\Boundingbox/"
    files = os.listdir(csv_dir)
    result_dir = r"E:\Biomass\CSIR\Result_0618\2018\shp/"
    for file in files:
        print(file)
        if "All" not in file:
            if "Addo" in file or "Venetia" in file or "DNyala" in file:
                Proj = Proj_dir + "UTM35S.prj"
            elif "Agulhas" in file:
                Proj = Proj_dir + "UTM34S.prj"
            else:
                Proj = Proj_dir + "UTM36S.prj"
            out_shp = result_dir + str(file).split(".")[0] + ".shp"
            out_buffer = result_dir + str(file).split(".")[0] + "_25m.shp"
            print("Converting csv to buffered shp: " + out_buffer)
            arcpy.management.XYTableToPoint(csv_dir + file, out_shp,
                                            "longitude", "latitude", '', Proj)
            arcpy.Buffer_analysis(out_shp, out_buffer, "12.5 Meters")



def biom_csv_clip_0619():
    arcpy.env.workspace = r"E:\Biomass\CSIR\Result_0618\2018\shp/"
    clip_dir = r"D:\temp\ArcGIS_project\ScanSAR_New\ScanSAR_New.gdb/"
    fc_list = arcpy.ListFeatureClasses()
    for fc in fc_list:
        if "25m" in fc:
            print(fc)
            out = r"E:\Biomass\CSIR\Result_0618\2018\shp/"
            arcpy.MakeFeatureLayer_management(clip_dir + fc.split("_")[0] + "_clip", 'clip')
            arcpy.MakeFeatureLayer_management(out + fc, fc.split(".")[0] + "_clip")
            arcpy.SelectLayerByLocation_management(fc.split(".")[0] + "_clip", 'INTERSECT', 'clip')
            arcpy.FeatureClassToShapefile_conversion(fc.split(".")[0] + "_clip", out)



def biom_csv_HCC_0619():
    arcpy.env.workspace = r"E:\Biomass\CSIR\Result_0618\2018\shp/"
    dir_chm = r"E:\ALS_archive\LiDAR_CHM_Mosaic05/"
    dir_cc = r"E:\ALS_archive\LiDAR_CHM_CC15/"
    dir_out = r"D:\temp\ArcGIS_project\Biomass_new\Biomass_new.gdb/"
    fc_list = arcpy.ListFeatureClasses()
    for fc in fc_list:
        if "_25m_clip.shp" in fc:
            print(fc)
            chm = dir_chm + fc.split("_")[0] + "05.tif"
            #cc = dir_cc + fc.split("_")[0] + ".tif"
            table_chm = dir_out + fc.split("_")[0] + "_H"
            #table_cc = dir_out + fc.split("_")[0] + "_CC"
            ZonalStatisticsAsTable(fc, "shot_numbe", chm, table_chm, statistics_type="MEAN")
            #ZonalStatisticsAsTable(fc, "shot_numbe", cc, table_cc, statistics_type="SUM")
            outTable_H = r"E:\Biomass\CSIR\Result_0618\2018\csv/" + fc.split("_")[0] + "_H_0.csv"
            #outTable_CC = r"E:\Biomass\CSIR\Result_0618\2018\csv/" + fc.split("_")[0] + "_CC.csv"
            arcpy.conversion.ExportTable(table_chm, outTable_H, "", "NOT_USE_ALIAS")
            #arcpy.conversion.ExportTable(table_cc, outTable_CC, "", "NOT_USE_ALIAS")


def biom_csv_merge_0619():
    df_initial_list = []
    df_initial = pd.DataFrame(df_initial_list)
    sites = ['Agincourt', 'DNyala', 'Welverdiendt', 'Venetia',
                 'Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh']
    dir = r"E:\Biomass\CSIR\Result_0618\2018\csv/"
    for site in sites:
        dir_GEDI = dir + site + ".csv"
        dir_H = dir + site + "_H_0.csv"
        dir_CC =  dir + site + "_CC.csv"
        df_GEDI = pd.read_csv(dir_GEDI)
        df_H = pd.read_csv(dir_H)
        df_CC = pd.read_csv(dir_CC)
        dir_H_sel = df_H.iloc[:, [1, 5]]
        dir_CC_sel = df_CC.iloc[:, [1, 5]]
        df_merge1 = pd.merge(df_GEDI, dir_H_sel, left_on='shot_number', right_on='shot_numbe', how='left')
        df_merge2 = pd.merge(df_merge1, dir_CC_sel, left_on='shot_number', right_on='shot_numbe', how='left')
        df_initial = pd.concat([df_initial, df_merge2])
    output = dir + "All_06192023_0.csv"
    df_initial.to_csv(output, index=None)



def CHM_CC_2018_25():
    sites = ['Agincourt', 'DNyala', 'Welverdiendt', 'Venetia',
                 'Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh']
    dir_shp_25 = r"E:\Biomass\CSIR\Result_0618\2018\shp/"
    dir_chm = r"E:\ALS_archive\LiDAR_CHM_Mosaic_05/"
    dir_CC = r"E:\ALS_archive\LiDAR_CHM_CC15/"
    dir_output_CHM = r"E:\Biomass\CSIR\Result_0618\2018\chm_25/"
    dir_output_CC = r"E:\Biomass\CSIR\Result_0618\2018\cc_25/"
    for site in sites:
        shp_25 = dir_shp_25 + site + "_c_25.shp"
        chm_1 = dir_chm + site + "_05.tif"
        CC_1 = dir_CC + site + ".tif"
        OutRas_CHM = ZonalStatistics(shp_25, "FID", chm_1, statistics_type="MEAN")
        OutRas_CHM.save(dir_output_CHM + site + ".tif")
        OutRas_CC = ZonalStatistics(shp_25, "FID", CC_1, statistics_type="SUM")
        OutRas_CC.save(dir_output_CC + site + ".tif")


def AGB_2018_25():
    sites = ['Agincourt', 'DNyala', 'Welverdiendt', 'Venetia',
             'Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh']
    dir_CHM = r"E:\Biomass\CSIR\Result_0618\2018\chm_25_05/"
    dir_CHM0 = r"E:\Biomass\CSIR\Result_0618\2018\chm_25/"
    dir_CC = r"E:\Biomass\CSIR\Result_0618\2018\cc_25/"
    dir_output_AGB = r"E:\Biomass\CSIR\Result_0618\2018\AGB/"
    for site in sites:
        print(site)
        CHM_ras = Raster(dir_CHM + site + ".tif")
        CHM0_ras = Raster(dir_CHM0 + site + ".tif")
        CC_ras = Raster(dir_CC + site + ".tif")
        #HCC -- linear
        #AGB_ras_linear = 10.35 * CHM_ras * CC_ras/490 - 5.9236
        #AGB_ras_linear.save(dir_output_AGB + "AGB_HCC_05_linear/" + site + ".tif")
        #HCC -- polynomial
        #AGB_ras_poly = 0.7694 * (CHM_ras * CC_ras/490) ** 2 + 3.9 * (CHM_ras * CC_ras/490) + 1.7568
        #AGB_ras_poly.save(dir_output_AGB + "AGB_HCC_05_poly/" + site + ".tif")
        #H -- >0.5m CHM
        #AGB_ras_05 = 0.3502 * CHM_ras ** 2.5405
        #AGB_ras_05.save(dir_output_AGB + "AGB_H_05/" + site + ".tif")
        #H original
        AGB_ras = 2.0997 * CHM0_ras ** 1.6593
        AGB_ras.save(dir_output_AGB + "AGB_H/" + site + ".tif")



def AGB_2018_25_0802():
    sites = ['Agincourt', 'Welverdiendt','Justicia', 'Ireagh']
    dir_CHM = r"E:\Biomass\CSIR\Result_0618\2018\chm_25_05/"
    dir_CC = r"E:\Biomass\CSIR\Result_0618\2018\cc_25/"
    dir_output_AGB = r"E:\Biomass\CSIR\Result_0802\SAS_HCC_2018/"
    for site in sites:
        print(site)
        CHM_ras = Raster(dir_CHM + site + ".tif")
        CC_ras = Raster(dir_CC + site + ".tif")
        #HCC -- linear
        AGB_ras = 9.0665 * CHM_ras * CC_ras/490
        AGB_ras.save(dir_output_AGB + site + ".tif")



def MCH_CC_25():
    arc_dir = r"D:\temp\ArcGIS_project\Biomass_new\Biomass_new.gdb/"
    CC_dir = r"E:\ALS_archive\LiDAR_CHM_CC15/"
    MCH_dir = r"E:\ALS_archive\LiDAR_CHM_Mosaic_05/"
    sites1 = ['Agincourt', 'DNyala', 'Welverdiendt', 'Venetia',
                 'Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh']
    sites = ['Justicia', 'Ireagh']
    for site in sites:
        print(site)
        Polygon_clip = arc_dir + site + "_c_25"
        Ras_CC = CC_dir + site + ".tif"
        Ras_MCH = MCH_dir + site + "_05.tif"
        table_CC = arc_dir + site + "_CC_25"
        table_MCH = arc_dir + site + "_MCH_25"

        ZonalStatisticsAsTable(Polygon_clip, "OID", Ras_CC, table_CC, "", "SUM")
        outTable_CC = r"E:\Biomass\CSIR\Result_01292024/" + site + "_CC.csv"
        arcpy.conversion.ExportTable(table_CC, outTable_CC, "", "NOT_USE_ALIAS")
        ZonalStatisticsAsTable(Polygon_clip, "OID", Ras_MCH, table_MCH, "", "MEAN")
        outTable_MCH = r"E:\Biomass\CSIR\Result_01292024/" + site + "_MCH.csv"
        arcpy.conversion.ExportTable(table_MCH, outTable_MCH, "", "NOT_USE_ALIAS")


def MCH_CC_25_I2():
    arc_dir = r"D:\temp\ArcGIS_project\Biomass_new\Biomass_new.gdb/"
    poly_dir = r"C:\Users\Shawn\Desktop\GEDI_0130/"
    CC_dir = r"E:\ALS_archive\LiDAR_CHM_CC15/"
    MCH_dir = r"E:\ALS_archive\LiDAR_CHM_Mosaic_05/"
    sites = ['Agincourt', 'DNyala', 'Welverdiendt', 'Venetia',
                 'Limpopo1', 'Limpopo2', 'Limpopo3', 'Justicia', 'Ireagh']
    for site in sites:
        print(site)
        Polygon_clip = poly_dir + "I2_100m.shp"
        Ras_CC = CC_dir + site + ".tif"
        Ras_MCH = MCH_dir + site + "_05.tif"
        table_CC = arc_dir + site + "_CC_I2_100m"
        table_MCH = arc_dir + site + "_MCH_I2_100m"
        ZonalStatisticsAsTable(Polygon_clip, "FID", Ras_CC, table_CC, "", "SUM")
        outTable_CC = r"C:\Users\Shawn\Desktop\GEDI_0130\I2/" + site + "_I2_100m_CC.csv"
        arcpy.conversion.ExportTable(table_CC, outTable_CC, "", "NOT_USE_ALIAS")
        ZonalStatisticsAsTable(Polygon_clip, "FID", Ras_MCH, table_MCH, "", "MEAN")
        outTable_MCH = r"C:\Users\Shawn\Desktop\GEDI_0130\I2/" + site + "_I2_100m_MCH.csv"
        arcpy.conversion.ExportTable(table_MCH, outTable_MCH, "", "NOT_USE_ALIAS")


