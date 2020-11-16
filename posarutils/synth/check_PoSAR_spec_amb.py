from posarutils.synth.SPECt import SPECt
import numpy as np

def check_PoSAR_spec_amb( CFG, ant_pos, GEOM ):

  # Re-calculate intermediate MEAN variables

  Y_Tx = np.mean( ant_pos.Tx_y )
  Y_Rx = np.mean( ant_pos.Rx_y )
  Z_Tx = np.mean( ant_pos.Tx_z )
  Z_Rx = np.mean( ant_pos.Rx_z )

  if( GEOM.gr_sr == 1 ): # slant range

    rmin_Tx = np.sqrt( (Y_Tx - GEOM.gr_pos_min)**2 + Z_Tx**2 )
    rmin_Rx = np.sqrt( (Y_Rx - GEOM.gr_pos_min)**2 + Z_Rx**2 )
    rmin = ( rmin_Tx + rmin_Rx ) / 2
    rmax_Tx = np.sqrt( ( Y_Tx - GEOM.gr_pos_max )**2 + Z_Tx**2 )
    rmax_Rx = np.sqrt( ( Y_Rx - GEOM.gr_pos_max )**2 + Z_Rx**2 )
    rmax = ( rmax_Tx + rmax_Rx ) / 2

    r = np.arange(rmin, rmax, GEOM.d_gr)

  else: # ground range

    rgr = np.arange(GEOM.gr_pos_min, GEOM.gr_pos_max, GEOM.d_gr)
    rmin_Tx = np.sqrt( np.amin( ( Y_Tx - rgr )**2 ) + Z_Tx**2 )
    rmin_Rx = np.sqrt( np.amin( ( Y_Rx - rgr )**2 ) + Z_Rx**2 )
    r = ( np.sqrt( ( Y_Tx - rgr )**2 + Z_Tx**2 ) + np.sqrt( ( Y_Rx - rgr )**2 + Z_Rx**2 ) ) / 2

  eps = np.finfo(dtype=float).eps

  GEOM.rg_res = ( GEOM.rg_res == 0 ) * eps + GEOM.rg_res

  #######
  # Range
  #######

  Delta_k_rg = 4 * np.pi * ( CFG.Fmax - CFG.Fmin ) / CFG.c  # Signal range bandwidth
  ksys_rg    = 4 * np.pi * ( CFG.Fmax - CFG.Fmin ) / CFG.c  # System range bandwidth
  kres_rg    = 2 * np.pi / GEOM.rg_res                  # Desired range bandwidth
  ks_rg      = 2 * np.pi / np.amax( np.diff( r ) )      # Image minimum range sampling rate


  # Original spectral properties
  o_SPEC = SPECt()

  o_SPEC.Delta_k_rg  = Delta_k_rg
  o_SPEC.ksys_rg     = ksys_rg
  o_SPEC.kres_rg     = kres_rg
  o_SPEC.ks_rg       = ks_rg

  # Check user provided values

  if(ksys_rg < Delta_k_rg):
    print("1,\n Error ! \n")
    print("System rg sampling < Signal rg bandwidth : nothing to do")

  if(kres_rg > Delta_k_rg):
    print("\n Warning ! \n Required range bandwidth superior to the signal one");
    print("--> Range resolution decreased to the signal best resolution");
    kres_rg = Delta_k_rg

  if(ks_rg < kres_rg):
    print("\n Warning ! \n Required range sampling unsufficient w.r.t. required range bandwidth");
    print("--> Range sampling interval decreased to resolution");
    ks_rg = kres_rg

  if(kres_rg < Delta_k_rg):
    print("\n Warning ! \n Required range bandwidth inferior to the signal one");
    print("--> Signal range bandwith has to be limited");
    Delta_k_rg = kres_rg;

  # Modified spectral properties

  SPEC = SPECt()

  SPEC.Delta_k_rg  = Delta_k_rg
  SPEC.ksys_rg     = ksys_rg
  SPEC.kres_rg     = kres_rg
  SPEC.ks_rg       = ks_rg

  ##########
  ## Azimuth
  ##########

  GEOM.az_res = ( GEOM.az_res == 0 ) * eps + GEOM.az_res

# Assumptions :
# The aperture length has been correctely limited before processing

# Rigourous computation k_az = k sin(phi) = k (x-xa) / np.sqrt( (x-xa)^2 + r^2 )
# Compute k_az in worst case : image lower corners and maximum frequency

  k_az0 = 4 * np.pi * CFG.Fmax / CFG.c
  
# first lower corner
# Old monostatic configuration
# Delta_x = GEOM.az_pos_min - [ min(ant_pos.x), max(ant_pos.x) ];
# k_az = k_az0 * Delta_x ./ np.sqrt( Delta_x.^2 + min(r)^2 );

  Delta_x_Tx = GEOM.az_pos_min - np.array( [ np.amin(ant_pos.Tx_x), np.amax(ant_pos.Tx_x) ] )
  k_az_Tx = k_az0 * Delta_x_Tx / np.sqrt( Delta_x_Tx**2 + rmin_Tx**2 )
  Delta_x_Rx = GEOM.az_pos_min - np.array( [ np.amin(ant_pos.Rx_x ), np.amax(ant_pos.Rx_x) ] )
  k_az_Rx = k_az0 * Delta_x_Rx / np.sqrt( Delta_x_Rx**2 + rmin_Rx**2 )
  k_az = ( k_az_Tx + k_az_Rx ) / 2

  print("(1) k_az " + str(k_az) )

# second lower corner
# Old monostatic configuration
# Delta_x = GEOM.az_pos_max - [ min(ant_pos.x), max(ant_pos.x) ];
# k_az = [ k_az; k_az0 * Delta_x ./ np.sqrt( Delta_x.^2 + min(r)^2 ) ];

  Delta_x_Tx = GEOM.az_pos_max - np.array( [ np.amin(ant_pos.Tx_x), np.amax(ant_pos.Tx_x) ] )
  k_az_Tx = k_az0 * Delta_x_Tx / np.sqrt( Delta_x_Tx**2 + rmin_Tx**2 )
  Delta_x_Rx = GEOM.az_pos_max - np.array( [ np.amin(ant_pos.Rx_x), np.amax(ant_pos.Rx_x) ] )
  k_az_Rx = k_az0 * Delta_x_Rx / np.sqrt( Delta_x_Rx**2 + rmin_Rx**2 )
  k_az = np.stack( (k_az, ( k_az_Tx + k_az_Rx ) / 2) )

  print("(2) k_az " + str(k_az) )

  Delta_k_az = np.amax( np.abs( np.diff( k_az, axis = 1 ) ) )         # Image pixel max azimuth bandwidth
  ksys_az = 2 * np.pi / np.mean( np.diff( ant_pos.mean_x ) )          # System azimuth bandwidth
  kres_az = 2 * np.pi /  GEOM.az_res                                  # Desired azimuth bandwidth
  ks_az = 2 * np.pi /  GEOM.d_az                                      # Image azimuth sampling rate
  Delta_k_az_max  = np.abs( np.amax( k_az[:] )- np.amin( k_az[:] ) )  # Image total azimuth bandwidth
  k_az_mean = np.sum( k_az, 1 ) / 2                                   # Central azimuth frequencies

  print("(3) k_az_mean " + str(k_az_mean) )
  print("(3) Delta_k_az " + str(Delta_k_az) )

  o_SPEC.Delta_k_az = Delta_k_az
  o_SPEC.ksys_az = ksys_az
  o_SPEC.kres_az = kres_az
  o_SPEC.ks_az = ks_az
  o_SPEC.Delta_k_az_max = Delta_k_az_max

  if(ksys_az<Delta_k_az_max):
    print("1,\n Error ! \n")
    #error('System az sampling < Image azimuth bandwidth : reduce image absolute dimensions');

  if(ksys_az>Delta_k_az_max):
    print("\n Warning ! \n System az sampling > Image azimuth bandwidth ");
    print("--> Raw Data could be down-sampled in the azimuth direction");

  if(kres_az>Delta_k_az):
    print("\n Warning ! \n Required azimuth bandwidth superior to the signal one");
    print("--> Azimuth resolution decreased to the signal best resolution");
    kres_az = Delta_k_az

  if(kres_az<Delta_k_az):
    print("\n Warning ! \n Required azimuth bandwidth inferior to the signal one ");
    print("--> Signal azimuth bandwidth has to be limited (angular integration)");
    Delta_k_az = kres_az


  k_az_max = np.amax( np.abs( k_az_mean ) ) + Delta_k_az / 2
  if( ks_az < 2 * k_az_max ):
    print("\n Warning ! \n Required azimuth sampling unsufficient w.r.t. required azimuth bandwidth");
    print("--> Azimuth sampling interval decreased")
    ks_az = 2 * k_az_max

# Effective spectral properties
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SPEC.Delta_k_az = Delta_k_az
  SPEC.ksys_az = ksys_az
  SPEC.kres_az = kres_az
  SPEC.ks_az = ks_az
  SPEC.Delta_k_az_max = Delta_k_az_max


  print("\nSummary");
  print("-------");

  print("\nRange characteristics")
  print("---------------------");
  print("range resolution : system " + str( 2 * np.pi / o_SPEC.ksys_rg )
    + ", required " + str( 2 * np.pi / o_SPEC.kres_rg )
    + ", effective " + str( 2 * np.pi / SPEC.kres_rg) )
  print("range sampling   : system " + str( 2 * np.pi / o_SPEC.ksys_rg )
    + ", required " + str( 2 * np.pi / o_SPEC.ks_rg )
    + ", effective " + str( 2 * np.pi / SPEC.ks_rg ) )

  print("\nAzimuth characteristics")
  print("-----------------------");
  print("azimuth resolution : system " + str( 2 * np.pi / o_SPEC.Delta_k_az)
    +", required " + str( 2 * np.pi / o_SPEC.kres_az )
    +", effective " + str( 2 * np.pi / SPEC.kres_az ) )
  print("azimuth sampling   : system " + str( 2 * np.pi / o_SPEC.ksys_az )
    +", required " + str( 2 * np.pi / o_SPEC.ks_az )
    +", effective " + str( 2 * np.pi / SPEC.ks_az ) )

  print("\n2*max azimuth frequency : " + str( 2 * k_az_max )
    + ", effective azimuth sampling: " + str( SPEC.ks_az ) )

  return (SPEC, o_SPEC)
