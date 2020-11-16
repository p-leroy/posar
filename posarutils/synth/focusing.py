import numpy as np
import numexpr as ne
import time

def d_im_calc( d_im, im_pos, ant_pos, pos ):

  # Tx
  ant_pos_Tx_x = ant_pos.Tx_x[ pos ]
  ant_pos_Tx_y = ant_pos.Tx_y[ pos ]
  ant_pos_Tx_z = ant_pos.Tx_z[ pos ]
  # Rx
  ant_pos_Rx_x = ant_pos.Rx_x[ pos ]
  ant_pos_Rx_y = ant_pos.Rx_y[ pos ]
  ant_pos_Rx_z = ant_pos.Rx_z[ pos ]

  el = im_pos.el
  az = im_pos.az
  gr = im_pos.gr

  ne.evaluate("( sqrt( (el - ant_pos_Tx_z)**2 + (az - ant_pos_Tx_x)**2 + (gr - ant_pos_Tx_y)**2 )"
    + "+ sqrt( (el - ant_pos_Rx_z)**2 + (az - ant_pos_Rx_x)**2 + (gr - ant_pos_Rx_y)**2 ) ) / 2", out = d_im)

def d_im_calc_ple( d_im, im_pos, ant_pos, pos ):

  xa = ant_pos.mean_x[pos]
  ya = ant_pos.mean_y[pos]
  za = ant_pos.mean_z[pos]

  el = im_pos.el
  az = im_pos.az
  gr = im_pos.gr

  ne.evaluate("sqrt( (az - xa)**2 + (gr - ya)**2 + (el - za)**2 )", out = d_im)

def sin_phi_calc( sin_phi, im_pos, ant_pos, d_im, pos ):

  #sin_phi = (im_pos.az - (ant_pos.Tx_x[ pos ] + ant_pos.Rx_x[ pos ]) ) ./ d_im

  ant_pos_Tx_x = ant_pos.Tx_x[ pos ]
  ant_pos_Rx_x = ant_pos.Rx_x[ pos ]

  az = im_pos.az

  sin_phi[:,:] = ne.evaluate("(az - (ant_pos_Tx_x + ant_pos_Rx_x) / 2 ) / d_im")

def focusing( CFG, RD, im_pos, ant_pos, PRM, focusedImage ):

  ne.set_num_threads(8)

  #May be useless, but ...
  d_rg = CFG.c / ( CFG.Fmax - CFG.Fmin ) / 2
  # d = nr;

  #Up-sampled distance
  up_d = -CFG.d_shift + np.arange( 0, CFG.n_up_smp, 1 ) * d_rg / CFG.up_smp
  # up_d=nr
  d_min = np.amin( up_d )
  d_max = np.amax( up_d )

  nbPointsImage = im_pos.el.size

  d_im    = np.zeros( im_pos.el.shape, dtype = np.float64 )
  yi      = np.zeros( im_pos.el.shape, dtype = np.complex128 )
  sin_phi = np.zeros( im_pos.el.shape, dtype = np.float64 )

  ind = np.zeros(nbPointsImage, dtype = int)
  ind_alt = np.zeros( im_pos.el.shape, dtype = bool )

  up_vec  = np.zeros( CFG.n_up_smp, dtype = np.complex128 )
  # Up sampled raw data vector (one azimuth position) in frequency domain
  vec_RD  = np.zeros( CFG.n_up_smp, dtype = np.complex128  )

  vec_ind = int( np.ceil( ( CFG.Nf + 1 ) / 2 ) )

  t = time.time()

  sin_phi_max = PRM.sin_phi_max

  threshold = 0

  for pos in range(CFG.Npos):
  #for pos in range(55700, 55701):

      aux = pos / CFG.Npos * 100
      if aux >= threshold:
        print( str(int(aux) ) + "% ..")
        threshold += 10

      #upsample raw data
      vec_RD[ 0 : vec_ind ] = RD[ pos, 0 : vec_ind ]
      vec_RD[ - ( CFG.Nf - vec_ind ) : ] = RD[ pos, vec_ind : ]
      # vec_RD contains RD(pos,:) and zeros: up-sampling

      up_vec = np.fft.ifft( vec_RD ) # Transformation to time domain

      #Compute image (focusing positions) delays at azimuth position xa
      d_im_calc_ple( d_im, im_pos, ant_pos, pos )

      # Compute azimuth look angle : sin( phi ) = ( x - xa ) / d
      #sin_phi_calc( sin_phi, im_pos, ant_pos, d_im, pos )
      ant_pos_Tx_x = ant_pos.Tx_x[ pos ]
      ant_pos_Rx_x = ant_pos.Rx_x[ pos ]
      az = im_pos.az
      ne.evaluate("(az - (ant_pos_Tx_x + ant_pos_Rx_x) / 2 ) / d_im", out = sin_phi)

      # Interpolate and focus raw data at image delay values
      # Lookfor effectively observed image positions & within the image bandwidth

      #ind = np.where( (d_im >= d_min) & (d_im <= d_max) & (np.abs( sin_phi ) < PRM.sin_phi_max) )
      ne.evaluate( "(d_im >= d_min) & (d_im <= d_max) & (abs( sin_phi ) < sin_phi_max)", out = ind_alt )

      # Interpolate raw data at focusing positions in time domain
      yi.fill( 0 )
      #interp_fix1d( up_d, up_vec, d_im, yi, ind_alt )
      yi[ ind_alt ] = np.interp( d_im[ind_alt], up_d, up_vec )

      # update image
      #focusedImage[ind_alt] += yi[ind_alt] * np.exp( 1j * CFG.kc * d_im[ ind_alt ] )
      kc = CFG.kc
      ne.evaluate("focusedImage + yi * exp( 1j * kc * d_im )", out = focusedImage)

  elapsed = time.time() - t
  print("execution time = " + str(elapsed))

  return d_im, sin_phi, ind_alt, vec_RD, yi, up_vec