"""
When switching projects:
    1. Copy `config.template.py` under `config.py` to avoid adding it to git.
    2. Modify the paths, names, and variables defined below. Pay extra caution
       to how the priors might get affected.
    3. Adapt `opencube.py`: it serves as a "middle-man" for all spectral cube
       I/O operations (data, uncertainties, [...and also spectral model?])
    4. Change the spectral model in the main sampler (`innocent_script.py`)
"""

import os

# [Project Settings]
sampler_script_file = 'innocent_script.py'  # Takes cmd arguments: NPEAKS, Y, X
name_id = 'H-MM1'
# If called from shell with no arguments, innocent_script.py will
# revert to these parameters:
default_yx = (80, 60) # Y, X position in the cube to extract a spectrum from
default_npeaks = 2

# [MultiNest configuration]
n_live_points = 400
sampling_efficiency = 0.8

# [Settings taken for the next sampler run]
lines = ['nh311', 'nh322']
line_names = ['oneone', 'twotwo']
npars = 6

# [Directory Paths]
proj_dir = os.path.expanduser('~/H-MM1_bayesian/')
chain_dir = os.path.join(proj_dir, 'nested-sampling/') # NOTE: needs ending "/"
logs_dir = os.path.join(proj_dir, 'nested-sampling/logs/')
cube_storage_dir = os.path.join(proj_dir, 'nested-sampling/cubexarr')

# [File Paths]
file_Ks = os.path.join(chain_dir, '{0}-Ks.fits'.format(name_id))
file_Zs = os.path.join(chain_dir, '{0}-Zs.fits'.format(name_id))
file_mle_formatter = os.path.join(chain_dir,
                                  '{}-mle'.format(name_id)+'-x{}.fits')
file_mle_x1 = file_mle_formatter.format(1) # that's prorbaly not the best way
file_mle_x2 = file_mle_formatter.format(2) # that's prorbaly not the best way
file_mle_x3 = file_mle_formatter.format(3) # that's prorbaly not the best way

# [File Paths: GAS DR1 fits files]
gasdata_dir = os.path.join(proj_dir, 'data')
file_nh311_dr1 = os.path.join(gasdata_dir, 'H-MM1_NH3_11_bls.fits')
file_nh322_dr1 = os.path.join(gasdata_dir, 'H-MM1_NH3_22_bls.fits')
file_rms_nh311_dr1 = os.path.join(gasdata_dir, 'H-MM1_rms_11.fits')
file_rms_nh322_dr1 = os.path.join(gasdata_dir, 'H-MM1_rms_22.fits')
file_sig_dr1 = os.path.join(gasdata_dir, 'sample_sig.fits')
file_esig_dr1 = os.path.join(gasdata_dir, 'sample_sig.fits')

# [Kwargs for saving data cube and xarr info (for lazy loading)]
cube_save_kwargs = dict(
    target_dir=cube_storage_dir,
    target_xarr='{}-xarr.npy'.format(name_id),
    target_xarrkwargs='{}-xarrkwargs.p'.format(name_id),
    target_cubefile='{}-data.npy'.format(name_id),
    target_errfile='{}-errors.npy'.format(name_id),
    target_header='{}-header.p'.format(name_id),
    mmap_mode='r')

# [Priors and where they came from ]
# - Temperatures: without a better knowledge, setting the priors to be as
#                 wide as we expect them to be - T_CMB to 25 (upper bound
#                 from the fact that we probably can't constrain it much
#                 further from the ammonia level population)
tk_prior = [5, 40]
tex_prior = [3, 15]
# - Line widths: from the thermal line width at 10 K
sigma_prior = [0.06, 0.8]
# - Total ammonia column: ranging from dex 12 to 15, for a typical NH3
#   abundance of 1e-8 we're tracing H2 densities of 1e20 to 1e23
ntot_prior = [13, 15]
# - Vlsr - from prior knowledge of cloud kinematics
xoff_prior = [3, 4.7]
dv_min, dv_max = 0.1, 0.5 # min/max separation of the velocity components
dxoff_prior = [dv_min, dv_max]
o2p_prior = [0.05, 1] # doesn't matter - will get fixed at 0.5 anyway

def get_priors_xoff_wrapped(npeaks):
    """
    Generates a list with prior edges for arbitrary spectral multiplicity

    What is xoff_wrapped, you ask? It's a workaround for the component
    switching via reparametrisation hack described in here:
    https://github.com/vlas-sokolov/pyspecnest/issues/1
    """
    priors_xoff_transformed = ((
              [tk_prior, tex_prior, ntot_prior, sigma_prior, dxoff_prior, o2p_prior])[::-1] * (npeaks - 1)
            + [tk_prior, tex_prior, ntot_prior, sigma_prior, xoff_prior, o2p_prior][::-1])

    return priors_xoff_transformed

