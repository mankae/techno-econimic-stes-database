import pandas as pd
import numpy as np

def data_import(file_path):
    '''Imports data from an Excel file and returns a DataFrame

    Excel file has to contain a sheet called 'data' with the following columns:
    - time: Date and Time
    - Q_storage (NEEDED FOR VALIDATION. FOR SIMULATION ONLY THE INITIAL VALUE IS NEEDED): Energy content of the storage (MWh)
    - Q_charge: Energy charged into the storage (MWh)
    - Q_discharge: Energy discharged from the storage (MWh)
    - T_air: Air temperature at the storage location (°C)
    - T_soil: Soil temperature at the storage location (°C)
    - T_XX (NEEDED FOR VALIDATION. FOR SIMULATION TEMPERATURE MAP IS USED): Temperature at different heights (increasing from bottom to top XX = 01, 02, ...) in the storage (°C)'''

    # data import
    df_data = pd.read_excel(file_path, sheet_name='data')
    df_data.set_index(df_data.columns[0], inplace=True)
    df_data.index = pd.to_datetime(df_data.index)

    return df_data

def extract_storage_temperature(data):
    '''Extracts the storage temperature data from the DataFrame and returns it as a NumPy array'''

    first_T = data.columns.get_loc('T_01')
    T_storage = data.iloc[:, first_T:].to_numpy()

    return T_storage

# def PTES_geometry_import(file_path):
#     '''Imports geometry data from an Excel file and returns the parameters
    
#     Excel file has to contain a sheet called 'geometry' with the following columns:
#     - h: Height of the truncated pyramid (m)
#     - a: Length of the upper rectangle (m)
#     - b: Width of the upper rectangle (m)
#     - c: Length of the lower rectangle (m)
#     - d: Width of the lower rectangle (m)'''

#     # geometry import
#     df_geometry = pd.read_excel(file_path, sheet_name='geometry')
#     h = df_geometry['h'][0]
#     a = df_geometry['a'][0]
#     b = df_geometry['b'][0]
#     c = df_geometry['c'][0]
#     d = df_geometry['d'][0]

#     return h, a, b, c, d