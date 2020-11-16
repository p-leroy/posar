import numpy as np
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np
# We now need to fix a datatype for our arrays. I've used the variable
# DTYPE for this, which is assigned to the usual NumPy runtime
# type info object.
DTYPE = np.float
# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef np.int_t DTYPE_t

def d_im_calc( np.ndarray[double, ndim=2] d_im, im_pos, ant_pos, pos ):

  nbArrows, nbColumns = im_pos.el.shape
  indices = nbArrows * nbColumns

  im_pos_el = 0.
  im_pos_az = 0.
  im_pos_gr = 0.

  # Tx
  ant_pos_Tx_x = ant_pos.Tx_x[ pos ]
  ant_pos_Tx_y = ant_pos.Tx_y[ pos ]
  ant_pos_Tx_z = ant_pos.Tx_z[ pos ]
  # Rx
  ant_pos_Rx_x = ant_pos.Rx_x[ pos ]
  ant_pos_Rx_y = ant_pos.Rx_y[ pos ]
  ant_pos_Rx_z = ant_pos.Rx_z[ pos ]

  for k in range(indices):
    im_pos_el = im_pos.el.flat[ k ]
    im_pos_az = im_pos.az.flat[ k ]
    im_pos_gr = im_pos.gr.flat[ k ]
    #d_im.flat[ k ] = (
    #  np.sqrt(
    #    (im_pos_el - ant_pos_Tx_z)**2
    #    + (im_pos_az - ant_pos_Tx_x)**2
    #    + (im_pos_gr - ant_pos_Tx_y)**2
    #    )
    #  + np.sqrt(
    #    (im_pos_el - ant_pos_Rx_z)**2
    #    + (im_pos_az - ant_pos_Rx_x)**2
    #    + (im_pos_gr - ant_pos_Rx_y)**2
    #    )
    #) / 2
