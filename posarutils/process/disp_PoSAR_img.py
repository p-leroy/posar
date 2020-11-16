import matplotlib.pyplot as plt
import numpy as np

class OPTt(object):
  """docstring for OPTt"""
  def __init__(self):
    super(OPTt, self).__init__()
    
    self.db = 0
    self.mod = 0
    self.med_dyn = 0
    self.cmap = "jet"
    self.title = "default title"

  def print(self):
    print("db " + str(self.db))
    print("mod " + str(self.mod))
    print("med_dyn " + str(self.med_dyn))
    print("cmap " + self.cmap )
    print("title " + self.title )

def disp_PoSAR_img( focusedImage, opt, az_pos, rg_gr_pos, im_extent, cmap = 'jet' ):

  opt.print()

  disp_im = focusedImage

  eps = np.finfo(dtype=float).eps

  if( opt.mod == 1 ):
    disp_im = np.abs( disp_im )
  elif ( opt.db == 1 ):
    disp_im = 20 * np.log10( np.abs( disp_im ) + eps )

  if( opt.med_dyn > 0 ):
    med = np.median( disp_im )
    disp_im = np.minimum( np.maximum( disp_im, med - opt.med_dyn / 2 ), med + opt.med_dyn / 2 )

  plt.figure()
  #imagesc( rg_gr_pos, az_pos, disp_im )

  im_max_db = 20 * np.log10( np.amax( np.abs( focusedImage ) ) )
  if im_extent == [0,0,0,0]:
    plt.imshow( disp_im, cmap=cmap )
  else:
    plt.imshow( disp_im, extent=im_extent, cmap=cmap )

  plt.colorbar()

  plt.title(opt.title)
  ax = plt.gca()
  ax.invert_xaxis()
  ax.xaxis.tick_top()
  ax.yaxis.tick_right()

  #colormap( opt.cmap )
  #if ( isfield( opt, 'title' ) )
  #  title( opt.title, 'interpreter', 'none' )

def processImg(img, mod, db, med_dyn, box):
    if mod:
        disp_im = np.abs( img )
    elif db:
        disp_im = 20 * np.log10( np.abs( img ) + eps )

    if( med_dyn > 0 ):
        med = np.median( disp_im )
        disp_im = np.minimum( 
            np.maximum( disp_im, med - med_dyn / 2 ), 
            med + med_dyn / 2 )
        
    return box_filter( disp_im , 2 )
