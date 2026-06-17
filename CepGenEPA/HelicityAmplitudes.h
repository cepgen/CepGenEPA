// Computes different helicity amplitudes as defined in
// Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787

#ifndef ggMatrixElements_HelicityAmplitudes_h
#define ggMatrixElements_HelicityAmplitudes_h

#include <complex>

int limits(double sred, double tred);

std::complex<double> Mxxxx_fermion(double x, double y);
std::complex<double> Mpppp_fermion(double sred, double tred, int exclude_loops);
std::complex<double> Mpmmp_fermion(double sred, double tred, int exclude_loops);
std::complex<double> Mpmpm_fermion(double sred, double tred, int exclude_loops);
std::complex<double> Mpppm_fermion(double sred, double tred, int exclude_loops);
std::complex<double> Mppmm_fermion(double sred, double tred, int exclude_loops);

std::complex<double> Mxxxx_vector(double x, double y);
std::complex<double> Mpppp_vector(double sred, double tred, int exclude_loops);
std::complex<double> Mpmmp_vector(double sred, double tred, int exclude_loops);
std::complex<double> Mpmpm_vector(double sred, double tred, int exclude_loops);
std::complex<double> Mpppm_vector(double sred, double tred, int exclude_loops);
std::complex<double> Mppmm_vector(double sred, double tred, int exclude_loops);

std::complex<double> Mxxxx_spin0even(double x, double y, double m, double f0, double w_const, double a2);
std::complex<double> Mpppp_spin0even(double x, double y, double m, double f0, double w_const, double a2);
std::complex<double> Mpmmp_spin0even(double x, double y, double m, double f0, double w_const, double a2);
std::complex<double> Mpmpm_spin0even(double x, double y, double m, double f0, double w_const, double a2);
std::complex<double> Mppmm_spin0even(double x, double y, double m, double f0, double w_const, double a2);
std::complex<double> Mpppm_spin0even(double x, double y, double m, double f0, double w_const, double a2);

std::complex<double> Mpppp_eft(double zeta1, double zeta2, double s, double t);
std::complex<double> Mpmmp_eft(double zeta1, double zeta2, double s, double t);
std::complex<double> Mpmpm_eft(double zeta1, double zeta2, double s, double t);
std::complex<double> Mppmm_eft(double zeta1, double zeta2, double s, double t);
std::complex<double> Mpppm_eft(double zeta1, double zeta2, double s, double t);

#endif
