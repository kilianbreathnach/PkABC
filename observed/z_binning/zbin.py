'''


Code to analyze redshift binning of galaxy catalog 

Author(s): ChangHoon Hahn 


'''
import numpy as np 
import os.path
import multiprocessing as mp

# --- Local ---
from Spectrum import data as spec_data
from Spectrum import fft as spec_fft
from Spectrum import spec as spec_spec

def zslice_data(DorR, cat_corr, **kwargs): 
    ''' Construct redshift binned galaxy/random catalog 

    Parameters
    ----------
    DorR : 'data' or 'random'
    cat_corr : Catalog and Correction dictionary. Specifies the catalog and correction parameters

    Notes
    -----
    * cat_corr example {'catalog': {'name': 'nseries', 'n_mock': 1}, 'correction': {'name': 'zbin1of5'}}
    * For the sake of simplicity, if 'of5' is specified in correction name the code will build 1,2,3,4,5 of 5 not 
    only the specified correction 
    * Example :
        cat_corr = {'catalog': {'name': 'nseries', 'n_mock': 1},
                'correction': {'name': 'zbin1of5'}}
        zslice_data('data', cat_corr)

    '''
    cat = cat_corr['catalog']       # catalog dictionary
    corr = cat_corr['correction']   # correction dictionary

    # catalog redshift limits
    if cat['name'].lower() in ('nseries'): 
        cat_zmin, cat_zmax = 0.43, 0.70
    else: 
        raise NotImplementedError("Only Nseries catalog has been implemented so far") 

    # import galaxy/random catalog to slice into redshift bins 
    def_cat_corr = {'catalog': cat, 'correction': {'name': 'default'}} 
    data = spec_data.Data(DorR, def_cat_corr, **kwargs) 
    data.Read(**kwargs)
    print 'Reading ', data.file_name 

    # slice catalog into redshifts 
    # interpret correction name (e.g. 1of5)
    corr_name = corr['name'].lower()
    n_slice = int(corr_name.split('of')[-1])
    slice_corr_names = [ 
            ''.join([ (corr_name.split('of')[0])[:-1], str(i_slice+1), 'of', str(n_slice) ]) 
            for i_slice in range(n_slice)
            ]
    
    # slice redshift limits 
    slice_zlow = np.arange(cat_zmin, cat_zmax, (cat_zmax - cat_zmin)/np.float(n_slice)) 
    slice_zhigh = np.append( np.arange(cat_zmin, cat_zmax, (cat_zmax - cat_zmin)/np.float(n_slice))[1:], cat_zmax )

    for i_z in range(n_slice): 
        # For redshift slices 
        z_slice = np.where( (data.z >= slice_zlow[i_z]) & (data.z < slice_zhigh[i_z]) ) 
    
        slice_cat_corr = {'catalog': cat, 'correction': {'name': slice_corr_names[i_z]}}
        slice_data = spec_data.Data(DorR, slice_cat_corr, **kwargs)
        for col in data.columns: 
            setattr(slice_data, col, getattr(data, col)[z_slice])
        print 'Writing ', slice_data.file_name
        slice_data.Write()

    return 

def Pk_zslice(mock_name, i_mock, n_slices, **kwargs):
    ''' Calculate P(k) for redshift sliced survey catalogs  

    Parameters
    ----------
    mock_name : name of mock catalog 
    i_mock : mock index (e.q. 1-84 for Nseries mocks)
    slice_str : string specifying the redshift slicing (e.g. "1of5")

    Notes
    -----
    * Only implemented for Nseries for now 
    * Multiprocessing function available to calculate P(k) simultaneously   
    * Calls spec from Spectrum module 

    '''
    if mock_name.lower() != 'nseries': 
        raise NotImplementedError("Only Nseries implemented")
    
    slice_corr_name = ''.join(['zbin1of', str(n_slices)])

    cat_corr = {
            'catalog': {'name': mock_name, 'n_mock': i_mock}, 
            'correction': {'name': slice_corr_name} 
            }
    
    # Check that file exists 
    R_data = spec_data.Data('random', cat_corr, **kwargs)
    D_data = spec_data.Data('data', cat_corr, **kwargs)
    if not os.path.isfile(R_data.file_name): 
        print ' '
        print 'Writing ', R_data.file_name
        print ' '
        zslice_data('random', cat_corr, **kwargs)
    if not os.path.isfile(D_data.file_name): 
        print ' '
        print 'Writing ', D_data.file_name
        print ' '
        zslice_data('data', cat_corr, **kwargs)

    spec = spec_spec.Spec('power', cat_corr, **kwargs) 
    spec.calculate()

    return None 

def Pk_zslice_mp(mock_name, n_mocks, n_slices, Nthreads=5, **kwargs):
    ''' Multiprocessing wrapper for Pk_zslice function which 
    calculates the P(k) for specified redshift slices 

    Paramaters
    ----------
    n_mocks : number of mocks to calculate

    '''

    pool = mp.Pool(processes=Nthreads)
    mapfn = pool.map

    arglist = [ [mock_name, i_mock, n_slices, kwargs] for i_mock in range(1, n_mocks+1) ] 

    mapfn(Pk_zslice, [arg for arg in arglist])

    pool.close()
    pool.terminate()
    pool.join()

    return

if __name__=='__main__': 
    Pk_zslice('nseries', 1, 5, Ngrid=360)
