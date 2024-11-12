#include <CepGen/Physics/Constants.h>
#include <CepGen/Physics/PDG.h>

#include <array>
#include <cmath>
#include <stdexcept>

#include "CepGenEPA/HelicityAmplitudes.h"
#include "CepGenEPA/MatrixElements.h"

namespace sm_aaaa {
  double prefac_W = -1.;
  const std::array<double, 9> SM_weight = {1, 1, 1, 16. / 27., 16. / 27., 16. / 27., 1. / 27., 1. / 27., 1. / 27.};
  std::array<double, 9> SM_masses;
  bool initialised = false;

  void me_SM(void (*me)(double, double, double *, double *, int),
             double s,
             double t,
             double *re,
             double *im,
             bool exclude_loops) {
    if (!initialised) {
      prefac_W = 0.25 / std::pow(cepgen::PDG::get().mass(23 /*W*/), 2);
      SM_masses = {cepgen::PDG::get().mass(11),
                   cepgen::PDG::get().mass(13),
                   cepgen::PDG::get().mass(15),
                   cepgen::PDG::get().mass(2),
                   cepgen::PDG::get().mass(4),
                   cepgen::PDG::get().mass(6),
                   cepgen::PDG::get().mass(1),
                   cepgen::PDG::get().mass(3),
                   cepgen::PDG::get().mass(5)};
      initialised = true;
    }
    // This routine computes the complex SM amplitude
    // The first argument can be any of the helicity amplitudes Mpppp,Mppmm,Mpmpm,Mpmmp,Mpppm

    // SM fermion content: (e,mu,tau,u,c,t,d,s,b)
    // SM_weight equals (number of colors) * (el. charge)^4
    // SM masses in GeV

    double d_re;
    double d_im;

    *re = 0;
    *im = 0;

    for (int i = 0; i <= 8; i++) {
      const auto prefac_f = 1. / (4 * SM_masses.at(i) * SM_masses.at(i));
      me(s * prefac_f, t * prefac_f, &d_re, &d_im, exclude_loops);
      *re += d_re * SM_weight.at(i);
      *im += d_im * SM_weight.at(i);
    }

    // Add also the W contribution
    if (me == Mpppp_fermion)
      Mpppp_vector(s * prefac_W, t * prefac_W, &d_re, &d_im, exclude_loops);
    else if (me == Mppmm_fermion)
      Mppmm_vector(s * prefac_W, t * prefac_W, &d_re, &d_im, exclude_loops);
    else if (me == Mpmpm_fermion)
      Mpmpm_vector(s * prefac_W, t * prefac_W, &d_re, &d_im, exclude_loops);
    else if (me == Mpmmp_fermion)
      Mpmmp_vector(s * prefac_W, t * prefac_W, &d_re, &d_im, exclude_loops);
    else if (me == Mpppm_fermion)
      Mpppm_vector(s * prefac_W, t * prefac_W, &d_re, &d_im, exclude_loops);

    *re += d_re;
    *im += d_im;

    *re *= 8 * cepgen::constants::ALPHA_EM * cepgen::constants::ALPHA_EM;
    *im *= 8 * cepgen::constants::ALPHA_EM * cepgen::constants::ALPHA_EM;

    // the factor of 8 is needed because of the conventions in
    // Costantini, DeTollis, Pistoni
  }

  // compute the SM squared matrix element, including leptons, quarks and the W boson
  double sqme(double s, double t, bool exclude_loops) {
    double re;
    double im;
    double value = 0;

    if (s < 0 || t > 0 || t < -s)
      throw std::runtime_error("Invalid domain. Valid range is s>=0 and -s<=t<=0");

    me_SM(Mpppm_fermion, s, t, &re, &im, exclude_loops);
    value += 4 * (re * re + im * im);

    me_SM(Mppmm_fermion, s, t, &re, &im, exclude_loops);
    value += re * re + im * im;

    me_SM(Mpppp_fermion, s, t, &re, &im, exclude_loops);
    value += re * re + im * im;

    me_SM(Mpmmp_fermion, s, t, &re, &im, exclude_loops);
    value += re * re + im * im;

    me_SM(Mpmpm_fermion, s, t, &re, &im, exclude_loops);
    value += re * re + im * im;

    return 0.5 * value;
  }

}  //namespace sm_aaaa
