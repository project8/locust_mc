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
        return 1.60217653E-19;
    } //!< elementary charge  in C(>0)

    constexpr static double Hbar()
    {
        return 1.05457168E-34;
    }//!< hbar in J s-1

    constexpr static double HbarC_eV()
    {
        return 197.326968E-9;
    }//!<hbar c in m eV.

    constexpr static double kB()
    {
        return 1.3806505E-23;
    }//!< Boltzmann constant J/K

    constexpr static double kB_eV()
    {
        return 8.617343E-5;
    }//!< Boltzmann constant eV/K

    //EM coupling constants
    constexpr static double EpsNull()
    {
        return 8.854187817E-12;
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
        return 9.1093826E-31;
    } //!< electron mass in kg

    constexpr static double M_el_eV()
    {
        return 510.998918E3;
    } //!< electron mass in ev/c^2

    constexpr static double M_p_kg()
    {
        return 1.67262192369E-27;
    } //!< proton mass in kg

    constexpr static double a0()
    {
        return 5.29177210903E-11;
    } //!< bohr radius in m

    constexpr static double E_Rydberg()
    {
        return 13.605693009;
    } //!< hydrogen ionization energy in eV

};

}

#endif //LMCCONST_HH
