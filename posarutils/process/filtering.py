import numpy as np
import scipy.signal as signal
import scipy.ndimage as ndimage

def filter2( H, X, mode='same' ):
    # Matlab: filter2( H, X, shape ) is equivalent to conv2( X, rot90( H, 2 ), shape )
    return signal.convolve2d( X, np.rot90( H, 2 ), mode=mode )


def filter2_alt( H, X ):
    # Matlab: filter2( H, X, shape ) is equivalent to conv2( X, rot90( H, 2 ), shape )
    return ndimage.convolve( X, np.rot90( H, 2 ) )

def box_filter( im, Nsmooth ):
    
    Nsmooth = Nsmooth + 1 - Nsmooth % 2 # Make Nsmooth an odd number
    
    filt = np.ones( (Nsmooth, Nsmooth) ) / Nsmooth**2

    si = im.shape
    filt_im = np.zeros( si )
    
    print( "im.shape = {}".format( im.shape ) )

    if len( si ) > 2:
        warning('Careful conventions have changed')

    case = len( si )
        
    if case == 2:
        filt_im = filter2( filt, im )
        
    elif case == 3 :
        for loop in range( si[ 2 ] ):
                filt_im[ :, :, loop ] = filter2( filt, im[ :, :, loop ] )
                
    elif case == 4:
        for loop1 in range( si[ 2 ] ):
            for loop2 in range( si[ 3 ] ): 
                filt_im[ :, :, loop1, loop2 ] = filter2( filt, im[ :, :, loop1, loop2 ] )
                
    else:
        print( 'Dimension error: size(im)= {}'.format( si.shape ) )
        
    return filt_im
