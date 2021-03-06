# -*- coding: utf-8 -*-
import xarray as xr
import numpy as np
import pandas as pd
import pathlib


class BaselineDatabase(object):
    def __init__(self):
        self.line_table = pd.DataFrame(columns=['site','install', 'line', 'instrument_id', 'comment'])
        self.instrument_table = pd.DataFrame(columns = ['instrument_id','type_id', 'sn', 'config'])

    def add2line_table(self, site,install_datetime, line_id, instrument_id, comment = ''):#20200205
        install_datetime = pd.to_datetime(install_datetime)
        new_line_table_entry = {'site':site,'install': install_datetime, 'line': line_id, 'instrument_id': instrument_id, 'comment': comment}
        self.line_table = self.line_table.append(new_line_table_entry, ignore_index=True)
        return
    
    def addnewinstrument(self, instrument_id, type_id, sn, config_id):
        self.instrument_table =  self.instrument_table.append({'instrument_id': instrument_id,'type_id':type_id, 'sn': sn, 'config':config_id}, ignore_index=True)
        return 
    
    def get_instrument(self, site, line, date):
    #     site = 'mlo'
    #     line = 121
#         date = df_all.index[0]
        lt_site_line = self.line_table[np.logical_and(self.line_table.site == 'mlo', self.line_table.line == line)]
        previous_installs = lt_site_line[lt_site_line.install <= date]
        if previous_installs.shape[0] == 0:
            raise IndexError(f'Instrument not installed (line:{line}, site: {site}, date: {date}')
        lt_found = previous_installs.iloc[-1]

        instrument_found = self.instrument_table[self.instrument_table.instrument_id == lt_found.instrument_id].iloc[0]
        return instrument_found
    
database = BaselineDatabase()
conf_1= {'A': 368, 'B': 1050, 'C': 610, 'D': 778}
conf_2= {'A': 412, 'B': 500, 'C': 675, 'D': 862}
database.addnewinstrument(1,1,1032,conf_2)
database.addnewinstrument(2,1,1046,conf_1)

installdate = '20131125'
database.add2line_table('mlo', installdate, 121, 2)
database.add2line_table('mlo', installdate, 221, 1)

database.add2line_table('mlo', '20200205', 121, 1)
database.add2line_table('mlo', '20200205', 221, 2)

database.add2line_table('mlo', '20200620', 121, 2)
database.add2line_table('mlo', '20200620', 221, 1)


def get_lines_from_station_header(path = '/nfs/grad/gradobs/documentation/station_headers/MLO_header.xlsx', line_ids = [121, 221]):
    path2header = pathlib.Path(path)

    df = pd.read_excel(path2header)

    col_names = {}
    lines = []
    for line_id in line_ids:
        idx = (df['Unnamed: 1'] == line_id).argmax()
        header = df.iloc[idx-1].dropna().values[1:]
        col_names[line_id] = header
        lines.append(dict(line_id = line_id, column_labels = header))
        
    return lines

def read_file(path2raw, lines = None, 
              # collabels = ['lineID', 'Year', 'DOY', 'HHMM', 'A', 'B', 'C', 'D','temp'],
              collabels = ['lineID', 'Year', 'DOY', 'HHMM'],
              database = None,
              site = None
              ):
    """
    

    Parameters
    ----------
    path2raw : TYPE
        DESCRIPTION.
    lines : list, optional
        List of lines to consider (e.g. [121, 221] for sp02 at MLO). The default is None -> all.
    collabels : TYPE, optional
        DESCRIPTION. The default is ['lineID', 'Year', 'DOY', 'HHMM'].
    database : TYPE, optional
        DESCRIPTION. The default is None.
    site : str, optional
        DESCRIPTION. The default is None. If None the site is infered from the
        file path. Set if the path is not the standard path

    Returns
    -------
    out_list : TYPE
        DESCRIPTION.

    """
    # the particular way I am reading here allows for later implementation of
    # reading old data from Longenecker. And also allows to read other raw files
    
    # if isinstance(instruement, type(None)):
    #     collabels = 
    # collabels = lines[0]['column_labels']
    out = {}
    collabels = np.array(collabels)
    
    df_all = pd.read_csv(path2raw, header=None
    #             names = False
               )
    # out['df_all_01'] = df_all.copy()
    colsis = df_all.columns.values
    colsis = [int(col) for col in colsis]
    
    # todo: assigne collumn labels accoreding to instrument info
    # if 0:
    colsis[:collabels.shape[0]] = collabels
    df_all.columns = colsis   
    # out['df_all_02'] = df_all.copy()
    # df_all = pd.read_csv(path2raw, names=lines[0]['column_labels'])

    # make datetime index
    df_all['HHMM'] = df_all.apply(lambda row: f'{int(row.HHMM):04d}', axis=1)
    df_all.index = df_all.apply(lambda row: pd.to_datetime(f'{int(row.Year)}0101') + pd.to_timedelta(row.DOY, 'd') + pd.to_timedelta(int(row.HHMM[:2]), 'h') + pd.to_timedelta(int(row.HHMM[2:]), 'm'), axis=1)
    df_all.index.name = 'datetime'
    
    data_list = []
    df_inst_temp = pd.DataFrame()
    df_inst_channels = pd.DataFrame()
    out['df_all'] = df_all.copy()
    # return out
    out_list = []
    date = df_all.index[0]
    # print(df_all.lineID.unique())
    for lid in df_all.lineID.unique():
        if isinstance(lines, list):
            if lid not in lines:
                print(f'{lid} not in lines ({lines})')
                continue
            
        df_lid = df_all[df_all.lineID == lid].copy()
        instrument = database.get_instrument(site, lid, date)
        df_lid = df_lid.drop(['lineID', 'Year','DOY', 'HHMM'], axis=1)
        df_lid.columns = ['A', 'B', 'C', 'D', 'temp']
    
        # replace invalid values with nan
        df_lid[df_lid == -99999] = np.nan
        df_lid[df_lid == -7999.0] = np.nan
    
        # seperate photo detector readings from temp 
        df_temp = df_lid.temp
        df_voltages = df_lid.reindex(['A', 'B', 'C', 'D'], axis= 1)
    
        df_voltages.columns.name = 'channel'
    
        # create dataset
        ds = xr.Dataset()
        ds['raw_data'] = df_voltages
        ds['internal_temperature'] = df_temp
        ser = pd.Series(instrument.config)
        ser.index.name = 'channel'
        ds['channle_wavelengths'] = ser
    
        ds['line_id'] = lid
        sn = instrument['sn']
        ds['serial_no'] = sn
    #     ds_by_instrument[f'sp02_{lid}_{sn}'] = ds
        out_list.append(ds)
        
    return out_list

#     for line in lines:
#         lid = line['line_id']
#         dft = df_all[df_all.lineID == lid].copy()
#         dft = dft.dropna(axis = 1)

#         # replace placeholder with correct column labels
#         dft.columns = line['column_labels']
#         line['df'] = dft.copy()

#         # clean up the channel voltages
#         df_channels = dft.drop(['lineID', 'Year', 'DOY', 'HHMM', 'SPO2 internal temp [degC]'], axis=1)
#         channels = [int(col.split(' ')[2]) for col in df_channels.columns]
#         df_channels.columns = channels
# #         df_channels.columns.name = f'wavelength_lid{lid}'
#         df_channels[df_channels == -99999] = np.nan
#         df_channels[df_channels == -7999.0] = np.nan
#         data_list.append(df_channels.copy())
#         # clean up temp
#         temp = dft['SPO2 internal temp [degC]'].copy()
#         temp[temp == -99999] = np.nan
#         temp[temp == -7999.0] = np.nan
#         df_inst_temp[lid] = temp
# #         print(len(channels))
# #         print(channels)
#         df_inst_channels[lid] = channels
# #         line['raw_data'] = df_channels
# #         ds[f'rawdata_line_id_{lid}'] = df_channels
# #         ds[f'instrument_temperature_line_id_{lid}'] = temp
# #     ds['line_ids'] = lines

#     ds = xr.Dataset()
#     data = pd.concat(data_list, axis=1).sort_index(axis=1)
#     data.columns.name = 'channel_wavelength'
#     ds['raw_data'] = data
#     df_inst_temp.columns.name = 'line_id'
#     ds['instrument_temperatures'] = df_inst_temp
#     df_inst_channels.columns.name = 'line_id'
#     df_inst_channels.index = [chr(ord('A') + i) for i in df_inst_channels.index]
#     df_inst_channels.index.name = 'channel'
#     ds['instrument_channels'] = df_inst_channels
#     return ds

def convert_raw2nc(path2rawfolder = '/nfs/grad/gradobs/raw/mlo/2020/', path2netcdf = '/mnt/telg/data/baseline/mlo/2020/', 
                   # database = None,
                   start_date = '2020-02-06',
                   pattern = '*sp02.*',
                   sernos = [1032, 1046],
                   overwrite = False, 
                   verbose = False, 
                   test = False):
    """
    

    Parameters
    ----------
    path2rawfolder : TYPE, optional
        DESCRIPTION. The default is '/nfs/grad/gradobs/raw/mlo/2020/'.
    path2netcdf : TYPE, optional
        DESCRIPTION. The default is '/mnt/telg/data/baseline/mlo/2020/'.
    # database : TYPE, optional
        DESCRIPTION. The default is None.
    start_date : TYPE, optional
        DESCRIPTION. The default is '2020-02-06'.
    pattern : str,  optional
        Only files with this pattern are considered. In newer raw data 
        versions this would be '*sp02.*'. In older ones: 'MLOD*'
    sernos : TYPE, optional
        DESCRIPTION. The default is [1032, 1046].
    overwrite : TYPE, optional
        DESCRIPTION. The default is False.
    verbose : TYPE, optional
        DESCRIPTION. The default is False.
    test : TYPE, optional
        If True only one file is processed. The default is False.

    Returns
    -------
    None.

    """
    # lines = get_lines_from_station_header()
    path2rawfolder = pathlib.Path(path2rawfolder)
    path2netcdf = pathlib.Path(path2netcdf)
    path2netcdf.mkdir(exist_ok=True)

    file_list = list(path2rawfolder.glob(pattern))
    # print(len(file_list))

    # file_contents = []
    return file_list
    
    df_in = pd.DataFrame(file_list, columns=['path_in'])
    df_in.index = df_in.path_in.apply(lambda x: pd.to_datetime(x.name.split('.')[2]))
    df_in.sort_index(inplace=True)
    df_in = df_in.truncate(before=start_date)
    
    df_out = pd.DataFrame(columns=['path_out'])
    
    for sn in sernos:
        for idx, row in df_in.iterrows():
            fnnc = row.path_in.name.replace('.dat','.nc')
            fnnc = fnnc.replace('-sp02', '.sp02')
            fnns = fnnc.split('.')
            fnns = fnns[:3] + [f'sn{str(sn)}'] + fnns[3:]
            fnnc = '.'.join(fnns)
            path2netcdf_file = path2netcdf.joinpath(fnnc)
    
            df_add = pd.DataFrame({'path_in': row.path_in, 'path_out':path2netcdf_file}, index = [idx]
        #                            ignore_index=True
                                  )
    
            df_out = df_out.append(df_add)
    
    df_out['exists'] = df_out.path_out.apply(lambda x: x.is_file())
    df_work = df_out[~df_out.exists]
    work_array = df_work.path_in.unique()

    print(f'No of files that need to be processed: {len(work_array)}')

    # exists = 0
    # new = 0
    for e, file in enumerate(work_array):
        # if e == 3: break
        # ds = read_file(file, lines)
        try:
            dslist = read_file(file, database = database, site = None)
        except IndexError:
            print('Instrument not installed ... skip', end = '...')
            continue
        ### generate output file name
        # processing
        for ds in dslist:
            fnnc = file.name.replace('.dat','.nc')
            fnnc = fnnc.replace('-sp02', '.sp02')
            fnns = fnnc.split('.')
            fnns = fnns[:3] + [f'sn{str(ds.serial_no.values)}'] + fnns[3:]
            fnnc = '.'.join(fnns)
            path2netcdf_file = path2netcdf.joinpath(fnnc)

            # ds = read_file(file, lines)
    
    
    
            # save to file
            ds.to_netcdf(path2netcdf_file)
            # new += 1
        if test:
            break
    # out = dict(processed = new,
    #            skipped = exists,
    #            last_ds_list = dslist)
    df_out['exists'] = df_out.path_out.apply(lambda x: x.is_file())
    df_work = df_out[~df_out.exists]
    work_array = df_work.path_in.unique()
    
    assert(df_work.shape[0] == 0)
    return 


