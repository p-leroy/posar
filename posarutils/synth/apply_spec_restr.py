import numpy as np

def apply_spec_restr( o_SPEC, SPEC, CFG, RD, GEOM, PRM):

  if PRM.proc_Dphi < 0 or PRM.proc_Dphi > np.pi:
    PRM.proc_Dphi = np.pi

  PRM.sin_phi_max = np.sin( PRM.proc_Dphi / 2 )

  if ( PRM.check == 1 ):

    #######
    # RANGE

    # Limit bandwidth
    k = np.linspace( -o_SPEC.ksys_rg / 2, o_SPEC.ksys_rg / 2, CFG.Nf )
    ind = find( np.abs(k) > SPEC.kres_rg / 2 )
    RD[ :,ind ] = 0
    GEOM.rg_res = 2 * pi / SPEC.kres_rg

    if GEOM.gr_sr:         # Slant range
      GEOM.d_gr = 2 * pi / SPEC.ks_rg
      
    else:                   # Ground range : worst case (far range)
      r_max = np.sqrt( CFG.H**2 + GEOM.gr_pos_max**2 )
      GEOM.d_gr = GEOM.gr_pos_max - np.sqrt( (r_max - 2 * np.pi / o_SPEC.ks_rg)**2 - CFG.H**2 )
      #clear r_max

    #########
    # AZIMUTH
    if(~PRM.force_res):
      PRM.sin_phi_max = 0.5 * SPEC.kres_az / CFG.kc
    # 2*k*sin_phi_max : ambiguity free image bandwidth

    GEOM.az_res = 2 * np.pi / min( SPEC.Delta_k_az, SPEC.kres_az ) # Effective resolution
    GEOM.d_az   = 2 * np.pi / SPEC.ks_az                           # Effective ambiguity free sampling

  return (GEOM, PRM)
