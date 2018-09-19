/*
 * wrapper.hpp
 *
 *  Created on: Oct 20, 2017
 *      Author: dmarce1
 */



#define NDIM 3
#define XDIM 2
#define YDIM 1
#define ZDIM 0



extern "C" {

void set_cpp_grid_( int* icell, double* x, double* y, double *z, int* nx, int* ny, int* nz, int* ncell );

void set_cpp_gas_(double* temp, double* eraddens, double* ur, double* rho,
		double* bcoef, double* decaygamma, double* vol, double* mass,
		double* ye, double* natom, double* nelec, double* natom1fr,
		double* matsrc, double* edep, int* nelem, int* niso);

void update_gas_cpp_(double dt);


int test();


}
