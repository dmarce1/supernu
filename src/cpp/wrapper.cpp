/*
 * wrapper.cpp
 *
 *  Created on: Oct 20, 2017
 *      Author: dmarce1
 */

#include "wrapper.hpp"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <array>
#include <vector>
#include <valarray>
#include <memory>

struct grid_t {
	std::array<std::vector<double>, NDIM> X;
	std::array<int, NDIM> N;
	std::vector<std::vector<std::vector<int>>>index;
	int ncell;
	grid_t(int* icell, double* x, double* y, double *z, int* nx,
			int* ny, int* nz, int* nc) {
		ncell = *nc;
		N[XDIM] = *nx;
		N[YDIM] = *ny;
		N[ZDIM] = *nz;
		for (int d = 0; d != NDIM; ++d) {
			X[d].resize(N[d]);
		}
		for (int i = 0; i != *nx; ++i) {
			X[XDIM][i] = x[i];
		}
		for (int i = 0; i != *nx; ++i) {
			X[YDIM][i] = y[i];
		}
		for (int i = 0; i != *nx; ++i) {
			X[ZDIM][i] = z[i];
		}
		index.resize(N[XDIM]);
		for( int i = 0; i != N[XDIM]; ++i) {
			index[i].resize(N[YDIM]);
			for( int j = 0; j != N[YDIM]; ++j) {
				index[i][j].resize(N[ZDIM]);
				for( int k = 0; k != N[ZDIM]; ++k) {
					const int fi = i + N[XDIM] * (j + N[YDIM] * k);
					index[i][j][k] = icell[fi];
				}
			}
		}
	}
};

struct gas_t {
	double* temp;
	double* eraddens;
	double* ur;
	double* rho;
	double* bcoef;
	double* decaygamma;
	double* vol;
	double* mass;
	double* ye;
	double* natom;
	double* nelec;
	double* natom_fr;
	double* niso_fr;
	double* matsrc;
	double* edep;
	int iso_cnt;
	int ele_cnt;
};

std::shared_ptr<grid_t> grid;
std::shared_ptr<gas_t> gas;

const double cons_a = 7.5657e-15;
const double cons_kb = 1.380658e-16;

extern "C" {

void set_cpp_grid_(int* icell, double* x, double* y, double *z, int* nx,
		int* ny, int* nz, int* nc) {
	grid = std::make_shared < grid_t > (icell, x, y, z, nx, ny, nz, nc);
}

}

double minmod(double a, double b) {
	return (std::copysign(0.5, a) + std::copysign(0.5, b))
			* std::min(std::abs(a), std::abs(b));
}

#define NF 5
#define BW 2

class array3d: public std::valarray<double> {
private:
	std::valarray<double>& v;
	std::array<int, NDIM> N;
	std::array<int, NDIM> dN;
	int sz;
	int index(int i, int j, int k) const {
		return i + N[XDIM] * (j + N[YDIM] * k);
	}
public:
	array3d() :
			v(*this) {
	}
	array3d(array3d&& other): v(*this) {
		*this = std::move(other.v);
	}
	array3d(const array3d& other): v(*this) {
		*this = other.v;
	}
	array3d(std::valarray<double>&& other): v(*this) {
		*this = std::move(other);
	}
	array3d(const std::valarray<double>& other): v(*this) {
		*this = other;
	}
	array3d(int nx, int ny, int nz) :
	v(*this), N( { {nx, ny, nz}}) {
		sz = nx * ny * nz;
		v.resize(sz);
		dN[XDIM] = 1;
		dN[YDIM] = nx;
		dN[ZDIM] = ny;
	}
	array3d& operator=(array3d&& other) {
		v = std::move(other.v);
		return *this;
	}
	array3d& operator=(const array3d& other) {
		v = other.v;
		return *this;
	}
	array3d& operator=(std::valarray<double>&& other) {
		v = std::move(other);
		return *this;
	}
	array3d& operator=(const std::valarray<double>& other) {
		v = other;
		return *this;
	}
	double& operator()(int i, int j, int k) {
		return v[index(i,j,k)];
	}
	double operator()(int i, int j, int k) const {
		return v[index(i,j,k)];
	}
	void reconstruct(array3d& p, array3d& m, int dim) const {

		for( int i = BW; i != sz - BW; ++i) {
			const double sp = v[i+dN[dim]] - v[i];
			const double sm = v[i] - v[i-dN[dim]];
			const double s = minmod(sp,sm);
			const double d = 0.5 * s;
			p[i] = v[i] - d;
			m[i + dN[dim]] = v[i] + d;
		}

	}
	array3d inv() const {
		array3d a;
		a.resize(sz);
		for( int i = 0; i != sz; ++i) {
			a[i] = 1.0 / v[i];
		}
		return a;
	}

	operator std::valarray<double>&() {
		return v;
	}

	operator const std::valarray<double>&()const {
		return v;
	}

	void zero() {
		for( int i = 0; i != sz; ++i) {
			v[i] = 0.0;
		}
	}

	array3d shift_p( int dim ) {
		array3d a;
		a.v.resize(sz);
		for( int i = 0; i != sz - dN[dim]; ++i) {
			a.v[i] = v[i + dN[dim]];
		}
		return a;
	}

	array3d shift_m( int dim ) {
		array3d a;
		a.v.resize(sz);
		for( int i = dN[dim]; i != sz; ++i) {
			a.v[i] = v[i - dN[dim]];
		}
		return a;
	}

	friend array3d max( const array3d& a, double b );
	friend array3d min( const array3d& a, double b );
	friend array3d max( const array3d& a, const array3d& b );
	friend array3d min( const array3d& a, const array3d& b );

	friend array3d max( const std::valarray<double>& a, double b );
	friend array3d min( const std::valarray<double>& a, double b );
	friend array3d max( const std::valarray<double>& a, const array3d& b );
	friend array3d min( const std::valarray<double>& a, const array3d& b );

};

array3d max(const array3d& a, double b) {
	array3d c;
	c.v.resize(a.sz);
	for (int i = 0; i != a.sz; ++i) {
		c.v[i] = std::max(a.v[i], b);
	}
	return c;
}

array3d max(const array3d& a, const array3d& b) {
	array3d c;
	c.v.resize(a.sz);
	for (int i = 0; i != a.sz; ++i) {
		c.v[i] = std::max(a.v[i], b.v[i]);
	}
	return c;
}

array3d min(const array3d& a, double b) {
	array3d c;
	c.v.resize(a.sz);
	for (int i = 0; i != a.sz; ++i) {
		c.v[i] = std::min(a.v[i], b);
	}
	return c;
}

array3d min(const array3d& a, const array3d& b) {
	array3d c;
	c.v.resize(a.sz);
	for (int i = 0; i != a.sz; ++i) {
		c.v[i] = std::min(a.v[i], b.v[i]);
	}
	return c;
}

array3d max(const std::valarray<double>& a, double b) {
	array3d c;
	c.v.resize(a.size());
	for (int i = 0; i != a.size(); ++i) {
		c.v[i] = std::max(a[i], b);
	}
	return c;
}

array3d min(const std::valarray<double>& a, double b) {
	array3d c;
	c.v.resize(a.size());
	for (int i = 0; i != a.size(); ++i) {
		c.v[i] = std::min(a[i], b);
	}
	return c;
}

struct hydro_vars_t {
	std::array<array3d, NDIM> X;
	std::array<array3d, NF> U;
	std::array<array3d, NF> dU;
	array3d& rho;
	array3d& sx;
	array3d& sy;
	array3d& sz;
	array3d& egas;
	array3d* s = &(U[1]);
	hydro_vars_t() :
			rho(U[0]), sx(U[1]), sy(U[2]), sz(U[3]), egas(U[4]), s(&(U[1])) {
		for (int f = 0; f != NF; ++f) {
			U[f] = array3d(nx, ny, nz);
			dU[f] = array3d(nx, ny, nz);
		}
	}

	void update(double dt) {

		const double gamma = 5.0 / 3.0;
		std::array<array3d, NDIM> v;
		array3d& vx = v[XDIM];
		array3d& vy = v[YDIM];
		array3d& vz = v[ZDIM];
		array3d rho_p(nx, ny, nz);
		array3d v_p[3];
		array3d& vx_p(v_p[XDIM]);
		array3d& vy_p(v_p[YDIM]);
		array3d& vz_p(v_p[ZDIM]);
		array3d ei_p(nx, ny, nz);
		array3d rho_m(nx, ny, nz);
		array3d v_m[3];
		array3d& vx_m(v_m[XDIM]);
		array3d& vy_m(v_m[YDIM]);
		array3d& vz_m(v_m[ZDIM]);
		array3d ei_m(nx, ny, nz);
		array3d s_f[3];
		array3d& sx(s_f[XDIM]);
		array3d& sy(s_f[YDIM]);
		array3d& sz(s_f[ZDIM]);
		array3d rho_f(nx, ny, nz);
		array3d sx_f(nx, ny, nz);
		array3d sy_f(nx, ny, nz);
		array3d sz_f(nx, ny, nz);
		array3d egas_f(nx, ny, nz);
		array3d sx_m(nx, ny, nz);
		array3d sy_m(nx, ny, nz);
		array3d sz_m(nx, ny, nz);
		array3d egas_m(nx, ny, nz);
		array3d sx_p(nx, ny, nz);
		array3d sy_p(nx, ny, nz);
		array3d sz_p(nx, ny, nz);
		array3d egas_p(nx, ny, nz);
		array3d U0[NF];
		array3d* U_f[] = { &rho_f, &sx_f, &sy_f, &sz_f, &egas_f };
		array3d* U_p[] = { &rho_p, &sx_p, &sy_p, &sz_p, &egas_p };
		array3d* U_m[] = { &rho_m, &sx_m, &sy_m, &sz_m, &egas_m };

		for (int f = 0; f != NF; ++f) {
			for (int i = BW; i != nx - BW; i++) {
				for (int j = BW; j != ny - BW; j++) {
					U[f](i, j, 0) = U[f](i, j, 2);
					U[f](i, j, 1) = U[f](i, j, 2);
					U[f](i, j, nz - 1) = U[f](i, j, nz - 3);
					U[f](i, j, nz - 2) = U[f](i, j, nz - 3);
				}
			}
			for (int i = BW; i != nx - BW; i++) {
				for (int j = BW; j != nz - BW; j++) {
					U[f](i, 0, j) = U[f](i, 2, j);
					U[f](i, 1, j) = U[f](i, 2, j);
					U[f](i, ny - 1, j) = U[f](i, ny - 3, j);
					U[f](i, ny - 2, j) = U[f](i, ny - 3, j);
				}
			}
			for (int i = BW; i != ny - BW; i++) {
				for (int j = BW; j != nz - BW; j++) {
					U[f](0, i, j) = U[f](2, i, j);
					U[f](1, i, j) = U[f](2, i, j);
					U[f](nx - 1, i, j) = U[f](nx - 3, i, j);
					U[f](nx - 2, i, j) = U[f](nx - 3, i, j);
				}
			}
		}

		for (int f = 0; f != NF; ++f) {
			U0[f] = U[f];
		}

		for (int rk = 0; rk != 2; ++rk) {

			const array3d rhoinv = rho.inv();
			const array3d ek(rhoinv * (sx * sx + sy * sy + sz * sz) * 0.5);
			const array3d ei(max(egas - ek, 0.0) * rhoinv);

			for (int d = 0; d != NDIM; ++d) {
				v[d] = s[d] * rhoinv;
			}

			for (int d = 0; d != NDIM; ++d) {
				rho.reconstruct(rho_p, rho_m, d);
				ei.reconstruct(ei_p, ei_m, d);
				for (int d1 = 0; d1 != NDIM; ++d1) {
					v[d1].reconstruct(v_p[d1], v_m[d1], d);
				}

				array3d p_p((gamma - 1.0) * rho_p * ei_p);
				array3d p_m((gamma - 1.0) * rho_m * ei_m);

				auto cs_p = sqrt(gamma * p_p / rho_p);
				auto cs_m = sqrt(gamma * p_m / rho_m);

				array3d a_p(abs(v_p[d]) + cs_p);
				array3d a_m(abs(v_m[d]) + cs_m);

				array3d a(max(a_p, a_m));

				auto ek_m = 0.5 * (vx_m * vx_m + vy_m * vy_m + vz_m * vz_m);
				auto ek_p = 0.5 * (vx_p * vx_p + vy_p * vy_p + vz_p * vz_p);

				sx_m = rho_m * vx_m;
				sy_m = rho_m * vy_m;
				sz_m = rho_m * vz_m;
				egas_m = rho_m * (ei_m + ek_m);

				sx_p = rho_p * vx_p;
				sy_p = rho_p * vy_p;
				sz_p = rho_p * vz_p;
				egas_p = rho_p * (ei_p + ek_p);

				for (int f = 0; f != NF; ++f) {
					*(U_f[f]) = 0.5 * (*(U_p[f]) * v_p[d] + *(U_m[f]) * v_m[d]);
					*(U_f[f]) -= 0.5 * a * (*(U_p[f]) - *(U_m[f]));
				}

				s_f[d] += 0.5 * (p_p + p_m);
				egas_f += 0.5 * (v_p[d] * p_p + v_m[d] * p_m);

				for (int f = 0; f != NF; ++f) {
					const auto fp = U_f[f]->shift_p(d);
					const auto& fm = *(U_f[f]);
					const auto xp = X[d].shift_p(d);
					const auto& xm = X[d];
					if (d == 0) {
						dU[f].zero();
					}
					dU[f] -= (fp - fm) / (xp - xm);
				}
			}
			double beta;
			if (rk == 0) {
				beta = 1.0;
			} else {
				beta = 0.5;
			}
			for (int f = 0; f != NF; ++f) {
				U[f] += (dt * beta) * dU[f] + (1.0 - beta) * (U0[f] - U[f]);
			}
		}

	}

	int nx;
	int ny;
	int nz;
};

hydro_vars_t vars;

extern "C" {

void set_cpp_gas_(double* temp, double* eraddens, double* ur, double* rho,
		double* bcoef, double* decaygamma, double* vol, double* mass,
		double* ye, double* natom, double* nelec, double* natom1fr,
		double* matsrc, double* edep, int* nelem, int* niso) {
	gas = std::make_shared<gas_t>();
	gas->temp = temp;
	gas->eraddens = eraddens;
	gas->ur = ur;
	gas->rho = rho;
	gas->bcoef = bcoef;
	gas->decaygamma = decaygamma;
	gas->vol = vol;
	gas->mass = mass;
	gas->ye = ye;
	gas->natom = natom;
	gas->nelec = nelec;
	gas->natom_fr = natom1fr + *niso + 1;
	gas->niso_fr = natom1fr;
	gas->matsrc = matsrc;
	gas->edep = edep;
	gas->iso_cnt = *niso;
	gas->ele_cnt = *nelem;

	vars.nx = grid->N[XDIM] + 2 * BW;
	vars.ny = grid->N[YDIM] + 2 * BW;
	vars.nz = grid->N[ZDIM] + 2 * BW;

	for (int i = BW; i != vars.nx - BW; ++i) {
		for (int j = BW; j != vars.ny - BW; ++j) {
			for (int k = BW; k != vars.nz - BW; ++k) {
				const int index = grid->index[i + BW][j + BW][k + BW];
				const double vinv = 1.0 / gas->vol[index];
				const double rho = gas->mass[index] * vinv;
				vars.rho(i, j, k) = rho;
				vars.sx(i, j, k) = grid->X[XDIM][i] * rho;
				vars.sy(i, j, k) = grid->X[YDIM][j] * rho;
				vars.sz(i, j, k) = grid->X[ZDIM][k] * rho;
				const double eint = cons_kb * gas->natom[index]
						* gas->temp[index] * vinv;
				const double ek = 0.5
						* (std::pow(vars.sx(i, j, k), 2)
								+ std::pow(vars.sy(i, j, k), 2)
								+ std::pow(vars.sz(i, j, k), 2)) / rho;
				vars.egas(i, j, k) = ek + eint;
				vars.X[XDIM](i, j, k) = grid->X[XDIM][i + BW];
				vars.X[YDIM](i, j, k) = grid->X[YDIM][j + BW];
				vars.X[ZDIM](i, j, k) = grid->X[ZDIM][k + BW];
			}
		}
	}

}

void update_gas_cpp_(double dt) {
}

}

