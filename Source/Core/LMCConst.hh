#ifndef LMCCONST_HH_
#define LMCCONST_HH_

#include <cmath>

namespace locust {

/**
 * This class contains various fundamental constants.
 * Values are taken from PDG edition 2006, unless pointed out otherwise. The naming conventions are: normal name for SI units, a suffix _unit for something else.
 **/
class LMCConst
{
public:
    LMCConst() = delete;

    //mathematical numbers
    template<class XFloatT = double>
    constexpr static XFloatT Pi()
    {
        return 3.141592653589793238462643383279502884L;
    } //!< pi

    constexpr static double C()
    {
        return 299792458.0;
    } //!< c im m/s

    constexpr static double Q()
    {
        return 1.602176634E-19;  // updated from https://pdg.lbl.gov/2021
    } //!< elementary charge  in C(>0)

    constexpr static double Hbar()
    {
        return 1.054571817E-34; // updated from https://pdg.lbl.gov/2021
    }//!< hbar in J s-1

    constexpr static double HbarC_eV()
    {
    	return 197.3269804; // updated from https://pdg.lbl.gov/2021
    }//!<hbar c in m eV.

    constexpr static double kB()
    {
        return 1.380649E-23; // updated from https://pdg.lbl.gov/2021
    }//!< Boltzmann constant J/K

    constexpr static double kB_eV()
    {
        return 8.617333262E-5; // updated from https://pdg.lbl.gov/2021
    }//!< Boltzmann constant eV/K

    //EM coupling constants
    constexpr static double EpsNull()
    {
        return 8.8541878128E-12; // updated from https://pdg.lbl.gov/2021
    } //!< epsilon0, Constant of Newtons force.

    constexpr static double FourPiEps()
    {
        return 4. * Pi() * EpsNull();
    } //!< 4  pi  epsilon0, Constant of Newtons force.

    constexpr static double MuNull()
    {
        return 4.E-7 * Pi();
    }//!< permeability of free space

    //masses
    constexpr static double M_el_kg()
    {
        return 9.1093837015E-31; // updated from https://pdg.lbl.gov/2021
    } //!< electron mass in kg

    constexpr static double M_el_eV()
    {
        return 510.998950E3; // updated from https://pdg.lbl.gov/2021
    } //!< electron mass in ev/c^2

    constexpr static double a0()
    {
        return 5.29177210903E-11;
    } //!< bohr radius in m

    constexpr static double E_Rydberg()
    {
        return 13.605693122994; // updated from https://pdg.lbl.gov/2021
    } //!< hydrogen ionization energy in eV

};

}

#endif //LMCCONST_HH
