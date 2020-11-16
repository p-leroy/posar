import numpy as np
import numexpr as ne

import numpy as np

def interp_fix1d( x, y, xi, yi, ind ):

	x_first = x[0]
	dx = x[1]-x[0]

	if ( dx < 0 ):
		error('x_last < x_first')
 
	xi_rel = ( xi[ind] - x_first ) / dx
	ind1 = np.floor( xi_rel ).astype(int)
	ind2 = np.ceil( xi_rel ).astype(int)

	x_ind1 = x[ ind1 ]
	y_ind1 = y[ ind1 ]
	y_ind2 = y[ ind2 ]
	xi_ind = xi[ ind ]

	yi[ind] = ne.evaluate( "y_ind1 + ( xi_ind - x_ind1 ) * ( y_ind2 - y_ind1 ) / dx" )

def interp_fix1d_alt( x, y, xi ):

	x_first = x[0]
	dx = x[1]-x[0]

	if ( dx < 0 ):
		error('x_last < x_first')
 
	xi_rel = ( xi - x_first ) / dx
	ind1 = np.floor( xi_rel ).astype(int)
	ind2 = np.ceil( xi_rel ).astype(int)

	x_ind1 = x[ ind1 ]
	y_ind1 = y[ ind1 ]
	y_ind2 = y[ ind2 ]

	yi = ne.evaluate( "y_ind1 + ( xi - x_ind1 ) * ( y_ind2 - y_ind1 ) / dx" )

	return yi

#function [yi]=interp_fix1d(x,y,xi)

#x = x(:);
#y = y(:);
#x_first = x(1);
#dx = x(2)-x(1);

#[l,c] = size(xi);
#xi = xi(:);

#if(dx<0)
# error('x_last<x_first');
#end 
#xi_rel = (xi-x_first)/dx;
#ind1 = floor(xi_rel)+1;
#ind2 = ceil(xi_rel)+1;

#yi = y(ind1) + (xi-x(ind1)).*(y(ind2)-y(ind1))/dx;

#yi = reshape(yi,l,c);
