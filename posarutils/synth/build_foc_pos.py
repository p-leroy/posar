from posarutils.synth.ImPosT import ImPosT
import numpy as np

def build_foc_pos( GEOM, ant_pos ):

# $$$ Reference track to be done later
    Y_Tx = np.mean( ant_pos.Tx_y )
    Y_Rx = np.mean( ant_pos.Rx_y )
    Z_Tx = np.mean( ant_pos.Tx_z )
    Z_Rx = np.mean( ant_pos.Rx_z )

    #############
    # SLANT RANGE
    #############

    if( GEOM.gr_sr == 1):
        
        rmin = (
            np.sqrt( (Y_Tx - GEOM.gr_pos_min)**2 + (Z_Tx - GEOM.el_pos)**2 ) +
            np.sqrt( (Y_Rx - GEOM.gr_pos_min)**2 + (Z_Rx - GEOM.el_pos)**2 )
            ) / 2

        rmax = (
            np.sqrt( (Y_Tx - GEOM.gr_pos_max)**2 + (Z_Tx - GEOM.el_pos)**2 ) +
            np.sqrt( (Y_Rx - GEOM.gr_pos_max)**2 + (Z_Rx - GEOM.el_pos)**2 )
            ) / 2

        r = np.arange(rmin, rmax, GEOM.d_gr)
        
  #%% Non valid for non null Y_Tx or Y_Rx
  #%%
  #vec_gr_pos = np.sqrt(
  #              r.^2 - (Z_Tx^2 + Z_Rx^2)/2
  #                   + ( (Z_Tx^2 - Z_Rx^2).^2 ) ./ (16*r.^2)
  #                   )

  #%% Equations are too long, numerical search instead

        x1 = np.ones( r.size ) * min( Y_Tx, Y_Rx )
        x2 = r + max( Y_Tx, Y_Rx )
        while( sum( abs( x1 - x2 ) > 1e-8 ) != 0 ):
            xtry = ( x1 + x2 ) / 2
            rtry = (
                np.sqrt( (Y_Tx - xtry)**2 + (Z_Tx - GEOM.el_pos)**2 ) + 
                np.sqrt( (Y_Rx - xtry)**2 + (Z_Rx - GEOM.el_pos)**2 )
                ) / 2
            comp = rtry > r
            x2 = x2 * (1-comp) + xtry * comp
            x1 = x1 * comp + xtry * (1-comp)

        vec_gr_pos = ( x1 + x2 ) / 2

        vec_rg_pos = r

    ##############
    # GROUND RANGE
    ##############

    else:
        vec_gr_pos = np.arange(GEOM.gr_pos_min, GEOM.gr_pos_max + GEOM.d_gr, GEOM.d_gr )
        vec_rg_pos = (
            np.sqrt( (Y_Tx - vec_gr_pos)**2 + (Z_Tx - GEOM.el_pos)**2 ) + 
            np.sqrt( (Y_Rx - vec_gr_pos)**2 + (Z_Rx - GEOM.el_pos)**2 )
        ) / 2

    vec_az_pos = np.arange(GEOM.az_pos_min, GEOM.az_pos_max + GEOM.d_az, GEOM.d_az )
                         
    im_pos = ImPosT()

    print("vec_az_pos " + str(vec_az_pos.size) + " ** vec_gr_pos " + str(vec_gr_pos.size ) )

    im_pos.gr = np.ones( (vec_az_pos.size, 1) )            * vec_gr_pos.reshape( (1, vec_gr_pos.size) )
    im_pos.az = vec_az_pos.reshape( (vec_az_pos.size, 1) ) * np.ones( (1, vec_gr_pos.size) )
    im_pos.el = np.ones( im_pos.gr.shape ) * GEOM.el_pos
    im_pos.vec_rg = vec_rg_pos
    im_pos.vec_gr = vec_gr_pos
    im_pos.vec_az = vec_az_pos

    im_pos.gr_alt = vec_gr_pos.reshape( (1, vec_gr_pos.size) )
    im_pos.az_alt = vec_az_pos.reshape( (vec_az_pos.size, 1) )
    im_pos.el_alt = np.ones( im_pos.gr_alt.shape ) * GEOM.el_pos

    return im_pos
