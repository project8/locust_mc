/*
 * LMCTrappedElectronGenerator.cc
 *
 *  Created on: Mar 4, 2015
 *      Author: plslocum after nsoblath
 */

#include "LMCTrappedElectronGenerator.hh"

#include "../Core/LMCLogger.hh"

using std::string;

namespace locust
{
    LMCLOGGER( lmclog, "TrappedElectronGenerator" );

    MT_REGISTER_GENERATOR(TrappedElectronGenerator, "trapped-electron");

    TrappedElectronGenerator::TrappedElectronGenerator( const std::string& aName ) :
            Generator( aName ),
            fDoGenerateFunc( &TrappedElectronGenerator::DoGenerateTime )
    {
        fRequiredSignalState = Signal::kTime;
    }

    TrappedElectronGenerator::~TrappedElectronGenerator()
    {
    }

    bool TrappedElectronGenerator::Configure( const ParamNode* aParam )
    {
        if( aParam == NULL) return true;

        if( aParam->Has( "domain" ) )
        {
            string domain = aParam->GetValue( "domain" );
            if( domain == "time" )
            {
                SetDomain( Signal::kTime );
                LMCDEBUG( lmclog, "Domain is equal to time.");
            }
            else if( domain == "freq" )
            {
                SetDomain( Signal::kFreq );
            }
            else
            {
                LMCERROR( lmclog, "Unable to use domain requested: <" << domain << ">" );
                return false;
            }
        }

        return true;
    }

    void TrappedElectronGenerator::Accept( GeneratorVisitor* aVisitor ) const
    {
        aVisitor->Visit( this );
        return;
    }

 
    Signal::State TrappedElectronGenerator::GetDomain() const
    {
        return fRequiredSignalState;
    }

    void TrappedElectronGenerator::SetDomain( Signal::State aDomain )
    {
        if( aDomain == fRequiredSignalState ) return;
        fRequiredSignalState = aDomain;  // pls changed == to =.
        if( fRequiredSignalState == Signal::kTime )
        {
            fDoGenerateFunc = &TrappedElectronGenerator::DoGenerateTime;
        }
        else if( fRequiredSignalState == Signal::kFreq )
        {
            fDoGenerateFunc = &TrappedElectronGenerator::DoGenerateFreq;
        }
        else
        {
            LMCWARN( lmclog, "Unknown domain requested: " << aDomain );
        }
        return;
    }


    bool TrappedElectronGenerator::DoGenerate( Signal* aSignal ) const
    {
        return (this->*fDoGenerateFunc)( aSignal );
    }


    bool TrappedElectronGenerator::DoGenerateTime( Signal* aSignal ) const
    {

        RunLengthCalculator *RunLengthCalculator1 = new RunLengthCalculator;
        double time = 0.;
        double StartingEnergy = 30.; // keV
        double StartingPitchAngle = 90.; // initial pitch angle in degrees 
        double TimeDependentEnergy = StartingEnergy;
        double LarmorPower = 0.;
        double TimeDependentAmplitude = 0.0; 
        double BBFreq = 0.;  // baseband frequency
        double ElectronStartTime = 0.000; // seconds
        double ElectronDuration = 0.001; // seconds

        double dt = RunLengthCalculator1->GetBinWidth(); // seconds
        double *StartPosition = StartElectron(StartingEnergy, StartingPitchAngle);
        double *position = StartPosition;
        double TimeDependentPitchAngle = StartingPitchAngle;
        double mu0 = GetMu(StartPosition, StartingEnergy, StartingPitchAngle);
        double CyclotronFrequency = CalculateCyclotronFrequency(CalculateGamma(TimeDependentEnergy), position);
        double ShiftedCyclotronFrequency = GetCyclotronFreqAntennaFrame(CyclotronFrequency, position[4]);
        double ShiftedCyclotronFrequencyAtShort = GetCyclotronFreqAntennaFrame(CyclotronFrequency, -position[4]);
        double phase1 = 2.*PI*(CENTER_TO_ANTENNA - position[2])/(C/ShiftedCyclotronFrequency);  // 2PI*L/lambda
        double phase2 = 2.*PI*(position[2] + 2.*CENTER_TO_SHORT + CENTER_TO_ANTENNA)/  // 2PI*L/lambda
          (C/ShiftedCyclotronFrequencyAtShort);
//printf("position[4] is %g\n", position[4]);
//printf("shifted cyc freq is %f\n", ShiftedCyclotronFrequencyAtShort);
printf("phase 1 is %g and phase 2 is %g radians, (phase1-phase2)/2PI is %g\n", phase1, phase2, (phase1-phase2)/2./PI);
 getchar();
        double LO_phase = 0.;
        double real_part1 = 0.;  // antenna
        double real_part2 = 0.;  // short
        double imaginary_part1 = 0.;
        double imaginary_part2 = 0.;

//printf("still outside loop\n");
//getchar();

        double TimeDependentMu = mu0;

        for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
        {
        time = (double)index/(RunLengthCalculator1->GetAcquisitionRate()*1.e6); // seconds
        if (time > ElectronStartTime && time < ElectronDuration)
        {
//        printf("taking a step\n\n");
        position = StepElectron(position, TimeDependentEnergy, TimeDependentMu, dt);  
        LarmorPower = CalculateLarmorPower(CalculateGamma(TimeDependentEnergy), GetBMag(position[0], position[1], position[2]));  // keV/s            
        CyclotronFrequency = CalculateCyclotronFrequency(CalculateGamma(TimeDependentEnergy), position);
        ShiftedCyclotronFrequency = GetCyclotronFreqAntennaFrame(CyclotronFrequency, position[4]);
        ShiftedCyclotronFrequencyAtShort = GetCyclotronFreqAntennaFrame(CyclotronFrequency, -position[4]);

        TimeDependentAmplitude = pow(2.e-15/2.,0.5);  // for debugging.
        phase1 += GetVoltagePhase(ShiftedCyclotronFrequency, dt);
        phase2 += GetVoltagePhaseFromShort(ShiftedCyclotronFrequencyAtShort, dt);  // reflecting short.
        LO_phase = -2.*PI*LO_FREQUENCY*time;

        real_part1 = cos(phase1)*cos(LO_phase) - sin(phase1)*sin(LO_phase);
        real_part2 = cos(phase2)*cos(LO_phase) - sin(phase2)*sin(LO_phase);
        imaginary_part1 = sin(phase1)*cos(LO_phase) + cos(phase1)*sin(LO_phase);
        imaginary_part2 = sin(phase2)*cos(LO_phase) + cos(phase2)*sin(LO_phase);

        if (abs((phase1+LO_phase)/2./PI/time)<100.e6)  // low pass filter
          {
          aSignal->SignalTime( index ) += TimeDependentAmplitude * (real_part1);
          aSignal->SignalTime( index ) += TimeDependentAmplitude * (real_part2);
          }

        }
        }



        delete RunLengthCalculator1;
        delete position;

        return true;
    }

    bool TrappedElectronGenerator::DoGenerateFreq( Signal* aSignal ) const
    {
        return true;
    }

    double TrappedElectronGenerator::GetCyclotronFreqAntennaFrame( double RFFreq, double Vparallel) const
    {
    double FPrime = 0.;
    double gamma = CalculateGamma(GetKineticEnergy(Vparallel));  // should replace C with C_w.
//    double gamma = 1.;
    FPrime = RFFreq * gamma * ( 1. - Vparallel/C );  // relativistic Doppler shift.
    return FPrime;
    }



    double TrappedElectronGenerator::GetVoltagePhase( double Freq, double dt ) const
    {    
    double phase = 0.;
    phase = 2.*PI*Freq*dt;
    return phase;
    }

    double TrappedElectronGenerator::GetVoltagePhaseFromShort( double Freq, double dt ) const
    {    
    double phase = 0.;
    phase = 2.*PI*Freq*dt;
    return phase;
    }





    double TrappedElectronGenerator::CalculateCyclotronFrequency( double Gamma, double *position ) const
    {
    double B = GetBMag(position[0], position[1], position[2]);
//    double B=GetBMag(0.,0.,0.); // constant B field.  
    double Frequency = 1.602e-19*B/(2.*PI*Gamma*9.11e-31);  // Hz
//    printf("B is %g\n", B);
    return Frequency;
    }

    double TrappedElectronGenerator::CalculateBasebandFrequency( double CyclotronFrequency ) const
    {
    double BasebandFrequency = (CyclotronFrequency - LO_FREQUENCY);  // Hz
    return BasebandFrequency;
    }


double TrappedElectronGenerator::GetMu(double *position, double KineticEnergy, double PitchAngle) const
{
double B = GetBMag(position[0], position[1], position[2]);
double Vperp = GetSpeed(KineticEnergy)*sin(PitchAngle*PI/180.);
double Eperp = GetKineticEnergy(Vperp);
double mu = Eperp/B;
printf("Eperp is %g\n", Eperp);
//printf("mu is %g and B is %g\n", mu, B);
//getchar();
return mu;
}



double TrappedElectronGenerator::CalculateLarmorPower( double gamma, double B) const  
{
double power = 1./(4.*PI*8.85e-12)*(2./3.)*pow(1.602e-19,4.)/pow(9.11e-31,2.)/3.e8*
   pow(B,2.)*(gamma*gamma-1.)*pow(sin(90.*PI/180.),2.)*(1./1.6027e-16);  // keV/s, theta = 90.
return power;
}



double *TrappedElectronGenerator::StartElectron(double KineticEnergy, double PitchAngle) const
{
double *position = new double[5];
position[0] = 0.;  // x in cm
position[1] = 0.;  // y in cm
position[2] = 0.;  // z in cm
double mu = GetMu(position, KineticEnergy, PitchAngle);
double Eperp = mu * GetBMag(position[0],position[1],position[2]);
double Eparallel = KineticEnergy - Eperp;
if (Eparallel < 0.) Eparallel = 0.;  

if (cos(PitchAngle*PI/180.) > 0.) position[3] = 1.;
else position[3] = -1.; 

double Vparallel = position[3]*GetSpeed(Eparallel);  // sign*speed = velocity

printf("vparallel is %g and Eparallel is %g\n", Vparallel, Eparallel);

position[4] = Vparallel;

return position;
}



bool TrappedElectronGenerator::StopElectron(double *position) const
{
delete position;
return true;
}

double *TrappedElectronGenerator::StepElectron(double *OldPosition, double KineticEnergy, double mu, double dt) const
{

double x = OldPosition[0];
double y = OldPosition[1];
double z = OldPosition[2];
double B0 = GetBMag(x,y,z);
double Eperp = mu * B0;
double Eparallel = KineticEnergy - Eperp;
double Vparallel = OldPosition[3]*GetSpeed(Eparallel);  // sign*speed = velocity
double dx = 0.;
double dy = 0.;
double dz = 0.;
double Eovershoot = 0.;
bool direction;

if (OldPosition[3] == 1.)
  direction = 1;
else if (OldPosition[3] == -1.)
  direction = 0;
else
  LMCERROR( lmclog, "Something is wrong with the position vector. ");

if (Eparallel > 0.)  // If we are not at a mirror point.
    {
    dx = GetBx(x,y,z)/GetBMag(x,y,z) * Vparallel * dt;  // push with the field direction.  
    dy = GetBy(x,y,z)/GetBMag(x,y,z) * Vparallel * dt;  // push with the field direction.
    dz = GetBz(x,y,z)/GetBMag(x,y,z) * Vparallel * dt;  // push with the field direction.
//    printf("z is %f and dz is %f and Eparallel is %f\n", z, dz, Eparallel);
    x += dx;
    y += dy;
    z += dz;
    }
else  // mirror point.
    {
//    printf("we seem to be at a mirror point.\n  Eparallel is %g and Vparallel is %g\n", Eparallel,Vparallel);
    Eovershoot = Eperp - KineticEnergy;
    Eperp = KineticEnergy - Eovershoot;  // new corrected Eperp on other side of KineticEnergy.
    Eparallel = KineticEnergy - Eperp;
    Vparallel = -OldPosition[3]*GetSpeed(Eparallel);  // flip velocity direction with new Eparallel.
    if (direction==0) direction = 1;
      else direction = 0;  // toggle direction variable.
    dx = GetBx(x,y,z)/GetBMag(x,y,z) * Vparallel * dt;  // push with the field direction.  
    dy = GetBy(x,y,z)/GetBMag(x,y,z) * Vparallel * dt;  // push with the field direction.
    dz = GetBz(x,y,z)/GetBMag(x,y,z) * Vparallel * dt;  // push with the field direction.
//    printf("z is %f and dz is %f and Eparallel is %f\n\n", z, dz, Eparallel);
    x += dx;
    y += dy;
    z += dz;
    }


//  deposit position and direction into return vector.
double *position = new double[5];
position[0] = x;
position[1] = y;
position[2] = z;
if (direction==1)
  position[3] = 1.;
else if (direction==0)
  position[3] = -1.;
else
  position[3] = 99.;
position[4] = Vparallel;

//printf("position 2 (z) is %g\n", position[2]);

return position;

}





double TrappedElectronGenerator::GetSpeed(double KineticEnergy) const
{
double Gamma = CalculateGamma(KineticEnergy);
double speed = pow((Gamma*Gamma - 1.)*C*C/(Gamma*Gamma),0.5); // cm/s
//printf("speed is %g cm/s\n", speed);
return speed;

}

double TrappedElectronGenerator::GetKineticEnergy(double Velocity) const  // Velocity in cm/s, energy in keV.
{
double KineticEnergy = 511.*(1./pow(1.-Velocity*Velocity/C/C,0.5) - 1.);
return KineticEnergy;
}

double TrappedElectronGenerator::CalculateGamma( double KineticEnergy ) const
{
double Gamma = 1. + KineticEnergy/511.;  // 511. keV
return Gamma;
}





double TrappedElectronGenerator::GetBMag(double x0, double y0, double z0) const
{

double Bmag = pow( GetBx(x0,y0,z0)*GetBx(x0,y0,z0) + 
                   GetBy(x0,y0,z0)*GetBy(x0,y0,z0) + 
                   GetBz(x0,y0,z0)*GetBz(x0,y0,z0), 0.5 );

return Bmag;
}





double TrappedElectronGenerator::GetBz( double x0, double y0, double z0 ) const
{

double Bz=0.;
int nsteps = 20;
double phi = 0.;
double dphi = 0.;
double dlx1, dly1, dlz1 = 0.;
double dlx2, dly2, dlz2 = 0.;
double xx1, xy1, xz1 = 0.;
double xx2, xy2, xz2 = 0.;
double R=RCOIL;

for (int i=0; i<nsteps; i++)
  {

  phi = (double)i/((double)nsteps)*2.*PI;
  dphi = 2.*PI/((double)nsteps);

  dlx1 = -sin(phi)*R*dphi;
  dly1 = cos(phi)*R*dphi;
  dlz1 = 0.;

  dlx2 = -sin(phi)*R*dphi;
  dly2 = cos(phi)*R*dphi;
  dlz2 = 0.;

  xx1 = x0 - R*cos(phi);
  xy1 = y0 - R*sin(phi);
  xz1 = z0 - (Z1);

  xx2 = x0 - R*cos(phi);
  xy2 = y0 - R*sin(phi);
  xz2 = z0 - Z2;



  Bz += MU0 * CURRENT * (dlx1*xy1 - dly1*xx1)/pow(xx1*xx1 + xy1*xy1 + xz1*xz1, 1.5);
  Bz += MU0 * CURRENT * (dlx2*xy2 - dly2*xx2)/pow(xx2*xx2 + xy2*xy2 + xz2*xz2, 1.5);

  }

Bz += 1.0; // Add in main 1 T field.

//printf("Bz is %g\n", Bz);

return Bz;  // T

}





double TrappedElectronGenerator::GetBx( double x0, double y0, double z0 ) const
{


double Bx=0.;
int nsteps = 20;
double phi = 0.;
double dphi = 0.;
double dlx1, dly1, dlz1 = 0.;
double dlx2, dly2, dlz2 = 0.;
double xx1, xy1, xz1 = 0.;
double xx2, xy2, xz2 = 0.;
double R=RCOIL;

for (int i=0; i<nsteps; i++)
  {

 

  phi = (double)i/((double)nsteps)*2.*PI;
  dphi = 2.*PI/((double)nsteps);

  dlx1 = -sin(phi)*R*dphi;
  dly1 = cos(phi)*R*dphi;
  dlz1 = 0.;

  dlx2 = -sin(phi)*R*dphi;
  dly2 = cos(phi)*R*dphi;
  dlz2 = 0.;

  xx1 = x0 - R*cos(phi);
  xy1 = y0 - R*sin(phi);
  xz1 = z0 - (Z1);

  xx2 = x0 - R*cos(phi);
  xy2 = y0 - R*sin(phi);
  xz2 = z0 - Z2;


  Bx += MU0 * CURRENT * (dly1*xz1 - dlz1*xy1)/pow(xx1*xx1 + xy1*xy1 + xz1*xz1, 1.5);
  Bx += MU0 * CURRENT * (dly2*xz2 - dlz2*xy2)/pow(xx2*xx2 + xy2*xy2 + xz2*xz2, 1.5);
//printf("Bx is now %f\n", Bx);


  }

//printf("Bx is %g\n", Bx);

return Bx;

}


double TrappedElectronGenerator::GetBy( double x0, double y0, double z0 ) const
{


double By=0.;
int nsteps = 20;
double phi = 0.;
double dphi = 0.;
double dlx1, dly1, dlz1 = 0.;
double dlx2, dly2, dlz2 = 0.;
double xx1, xy1, xz1 = 0.;
double xx2, xy2, xz2 = 0.;
double R=RCOIL;

for (int i=0; i<nsteps; i++)
  {

 


  phi = (double)i/((double)nsteps)*2.*PI;
  dphi = 2.*PI/((double)nsteps);

  dlx1 = -sin(phi)*R*dphi;
  dly1 = cos(phi)*R*dphi;
  dlz1 = 0.;

  dlx2 = -sin(phi)*R*dphi;
  dly2 = cos(phi)*R*dphi;
  dlz2 = 0.;

  xx1 = x0 - R*cos(phi);
  xy1 = y0 - R*sin(phi);
  xz1 = z0 - (Z1); 

  xx2 = x0 - R*cos(phi);
  xy2 = y0 - R*sin(phi);
  xz2 = z0 - Z2;

  By += MU0 * CURRENT * (dlz1*xx1 - dlx1*xz1)/pow(xx1*xx1 + xy1*xy1 + xz1*xz1, 1.5);
  By += MU0 * CURRENT * (dlz2*xx2 - dlx2*xz2)/pow(xx2*xx2 + xy2*xy2 + xz2*xz2, 1.5);

//printf("By is %g\n", By);
  }

//printf("By is %g\n", By);

return By;

}






} /* namespace locust */
