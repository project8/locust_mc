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
        BuildFieldMaps();
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
        double StartingPitchAngle = 90.0; // initial pitch angle in degrees 
        double TimeDependentEnergy = StartingEnergy;
        double LarmorPower = 0.;
        double TimeDependentAmplitude = 0.0; 
        double BBFreq = 0.;  // baseband frequency
        double ElectronStartTime = 0.000; // seconds
        double ElectronDuration = 8.375e-5; // seconds  (5 triggers on scope).

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
// getchar();
        double LO_phase = 0.;
        double real_part1 = 0.;  // antenna
        double real_part2 = 0.;  // short
        double imaginary_part1 = 0.;  // antenna
        double imaginary_part2 = 0.;  // short
        double Eloss = 0.;  // radiative energy loss.  keV.

//printf("still outside loop\n");
//getchar();

        double TimeDependentMu = mu0;



        for( unsigned index = 0; index < aSignal->TimeSize(); ++index )
        {
        time = (double)index/(RunLengthCalculator1->GetAcquisitionRate()*1.e6); // seconds
        if (time > ElectronStartTime && time < ElectronDuration+ElectronStartTime)
        {
        LarmorPower = CalculateLarmorPower(CalculateGamma(TimeDependentEnergy), GetBMag(position[0], position[1], position[2]));  // keV/s
// Now lose LarmorPower*dt in keV          
        Eloss = LarmorPower*RunLengthCalculator1->GetBinWidth(); 
        TimeDependentEnergy -= Eloss;

//        TimeDependentAmplitude = pow(2.e-15,0.5);  // arbitrary, for debugging.
//        TimeDependentAmplitude = pow(LarmorPower*1.602e-16, 0.5); // keV/s * 1.6e-19 J/eV * 1.e3 eV/keV = W.
        TimeDependentAmplitude = pow(LarmorPower/2.*1.602e-16, 0.5); // keV/s * 1.6e-19 J/eV * 1.e3 eV/keV = W.
// this next line is where the energy loss affects the tracking:
        TimeDependentMu = TimeDependentMu* (1.-Eloss/TimeDependentEnergy);
            
        CyclotronFrequency = CalculateCyclotronFrequency(CalculateGamma(TimeDependentEnergy), position);
        printf("cyc freq - LO_FREQUENCY is %g\n", CyclotronFrequency-LO_FREQUENCY); getchar();
        ShiftedCyclotronFrequency = GetCyclotronFreqAntennaFrame(CyclotronFrequency, position[4]);
        ShiftedCyclotronFrequencyAtShort = GetCyclotronFreqAntennaFrame(CyclotronFrequency, -position[4]);

        phase1 += GetVoltagePhase(ShiftedCyclotronFrequency, dt);
        phase2 += GetVoltagePhaseFromShort(ShiftedCyclotronFrequencyAtShort, dt);  // reflecting short.
        LO_phase = -2.*PI*LO_FREQUENCY*time;

        real_part1 = cos(phase1)*cos(LO_phase) - sin(phase1)*sin(LO_phase);
        real_part2 = cos(phase2)*cos(LO_phase) - sin(phase2)*sin(LO_phase);
//        imaginary_part1 = sin(phase1)*cos(LO_phase) + cos(phase1)*sin(LO_phase);
//        imaginary_part2 = sin(phase2)*cos(LO_phase) + cos(phase2)*sin(LO_phase);


//        if ( fabs(acos(real_part1)/2./PI/time) < 100.e6 )  // low pass filter.  wrong.
//          {
//          if (fabs(acos(real_part1)/2./PI/time) > 50.e6)
//            printf("bb freq is %g\n", fabs( acos(real_part1)/2./PI/time ) );
//          printf("time is %g\n", time); getchar();
          aSignal->SignalTime( index ) += TimeDependentAmplitude * (real_part1);  
//          aSignal->SignalTime( index ) += TimeDependentAmplitude * (real_part2);  // signal from short.
//          printf("now signal is %g\n\n", aSignal->SignalTime(index)); getchar();
//          }  // lpf

//        printf("taking a step\n\n");
        position = StepElectron(position, TimeDependentEnergy, TimeDependentMu, dt);  
//        printf("z is %f\n", position[2]); getchar();


        }  // if start<time<end

        
        }  // index



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
    printf("B is %g and Gamma is %f and fcyc is %g and fcyc-LO is %g\n", B, Gamma, Frequency,Frequency-LO_FREQUENCY); getchar();
    return Frequency;
    }

// naive.
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
//printf("Eperp is %g\n", Eperp);
//printf("mu is %g and B is %g\n", mu, B);
//getchar();
return mu;
}



double TrappedElectronGenerator::CalculateLarmorPower( double gamma, double B) const  
{
double power = 1./(4.*PI*8.85e-12)*(2./3.)*pow(1.602e-19,4.)/pow(9.11e-31,2.)/3.e8*
   pow(B,2.)*(gamma*gamma-1.)*pow(sin(90.0*PI/180.),2.)*(1./1.6027e-16);  // keV/s, theta = 90.
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

printf("start electron says vparallel is %g and Eparallel is %g and Bmag is %g\n", Vparallel, Eparallel, GetBMag(position[0],position[1],position[2]));

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
//    Eovershoot = Eperp - KineticEnergy;
//    Eperp = KineticEnergy - Eovershoot;  // new corrected Eperp on other side of KineticEnergy.
//    Eparallel = KineticEnergy - Eperp;

    Eparallel = -Eparallel;
    Eperp = KineticEnergy - Eparallel;

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

printf("position[2] (z) is %g\n", position[2]); getchar();

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

//double TrappedElectronGenerator::GenerateBzMap(double x0, double y0, double z0) const
//{
//}



double TrappedElectronGenerator::GetBMag(double x0, double y0, double z0) const
{

double Bmag = pow( GetBx(x0,y0,z0)*GetBx(x0,y0,z0) + 
                   GetBy(x0,y0,z0)*GetBy(x0,y0,z0) + 
                   GetBz(x0,y0,z0)*GetBz(x0,y0,z0), 0.5 );

return Bmag;
}


double TrappedElectronGenerator::GetBz(double x0, double y0, double z0) const
{


double Bz=0.;
int nsteps = 20;
double ZBottom = ZBOTTOM;
double ZTop = ZTOP;
double ncoils = NCOILS;
double phi = 0.;
double dphi = 0.;
double dlx, dly, dlz = 0.;
double xx, xy, xz = 0.;
double R=RSOLENOID;
double Zn = 0.;  // Z position of nth coil (cm).


for (int j=0; j<ncoils; j++)  // step through coil windings
{
//printf("j is %d\n", j);
Zn = ZBottom + (double)j*(ZTop-ZBottom)/((double)ncoils-1.);
for (int i=0; i<nsteps; i++)
  {

  phi = (double)i/((double)nsteps)*2.*PI;
  dphi = 2.*PI/((double)nsteps);


// small piece of coil #n.
  dlx = -sin(phi)*R*dphi;
  dly = cos(phi)*R*dphi;
  dlz = 0.;

// vector from location of small piece of coil #n to observation point x0 y0 z0.
  xx = x0 - R*cos(phi);
  xy = y0 - R*sin(phi);
  xz = z0 - (Zn);  // Zn is the z position of the nth coil.


// B at x0 y0 z0 from small piece of coil #n
if (xx*xx + xy*xy + xz*xz > 0.)
  Bz += 1./(4.*PI) * MU0 * CURRENT_SOLENOID * (dlx*xy - dly*xx)/pow(xx*xx + xy*xy + xz*xz, 1.5);

//printf("Bz is now %f\n", Bz);


  } // i (phi)

} // j (coil windings)


Bz += 1.0;  // add in main field.
//printf("Bz is %g\n", Bz);



return Bz;

}




double TrappedElectronGenerator::GetBx(double x0, double y0, double z0) const
{


double Bx=0.;
int nsteps = 20;
double ZBottom = ZBOTTOM;
double ZTop = ZTOP;
double ncoils = NCOILS;
double phi = 0.;
double dphi = 0.;
double dlx, dly, dlz = 0.;
double xx, xy, xz = 0.;
double R=RSOLENOID;
double Zn = 0.;  // Z position of nth coil (cm).


for (int j=0; j<ncoils; j++)  // step through coil windings
{
//printf("j is %d\n", j);
Zn = ZBottom + (double)j*(ZTop-ZBottom)/((double)ncoils-1.);
for (int i=0; i<nsteps; i++)
  {

  phi = (double)i/((double)nsteps)*2.*PI;
  dphi = 2.*PI/((double)nsteps);


// small piece of coil #n.
  dlx = -sin(phi)*R*dphi;
  dly = cos(phi)*R*dphi;
  dlz = 0.;

// vector from location of small piece of coil #n to observation point x0 y0 z0.
  xx = x0 - R*cos(phi);
  xy = y0 - R*sin(phi);
  xz = z0 - (Zn);  // Zn is the z position of the nth coil.


// B at x0 y0 z0 from small piece of coil #n
if (xx*xx + xy*xy + xz*xz > 0.)
  Bx += 1./(4.*PI) * MU0 * CURRENT_SOLENOID * (dly*xz - dlz*xy)/pow(xx*xx + xy*xy + xz*xz, 1.5);


//printf("Bx is now %f\n", Bx);


  } // i (phi)

} // j (coil windings)



return Bx;

}




double TrappedElectronGenerator::GetBy(double x0, double y0, double z0) const
{


double By=0.;
int nsteps = 20;
double ZBottom = ZBOTTOM;
double ZTop = ZTOP;
double ncoils = NCOILS;
double phi = 0.;
double dphi = 0.;
double dlx, dly, dlz = 0.;
double xx, xy, xz = 0.;
double R=RSOLENOID;
double Zn = 0.;  // Z position of nth coil (cm).


for (int j=0; j<ncoils; j++)  // step through coil windings
{
//printf("j is %d\n", j);
Zn = ZBottom + (double)j*(ZTop-ZBottom)/((double)ncoils-1.);
for (int i=0; i<nsteps; i++)
  {

  phi = (double)i/((double)nsteps)*2.*PI;
  dphi = 2.*PI/((double)nsteps);


// small piece of coil #n.
  dlx = -sin(phi)*R*dphi;
  dly = cos(phi)*R*dphi;
  dlz = 0.;

// vector from location of small piece of coil #n to observation point x0 y0 z0.
  xx = x0 - R*cos(phi);
  xy = y0 - R*sin(phi);
  xz = z0 - (Zn);  // Zn is the z position of the nth coil.


// B at x0 y0 z0 from small piece of coil #n
if (xx*xx + xy*xy + xz*xz > 0.)
  By += 1./(4.*PI) * MU0 * CURRENT_SOLENOID * (dlz*xx - dlx*xz)/pow(xx*xx + xy*xy + xz*xz, 1.5);

//printf("By is now %f\n", By);


  } // i (phi)

} // j (coil windings)

//printf("By is %g\n", By);



return By;

}



double TrappedElectronGenerator::InterpolateB(double x0, double y0, double z0, double *fieldmap) const
{

double space_dimension = SPACE_DIMENSION;
double space_resolution = SPACE_RESOLUTION;
int array_dimension = (int)(space_dimension/space_resolution);

int xindex = floor((x0 + space_dimension/2.)/space_resolution); // + = --
int yindex = floor((y0 + space_dimension/2.)/space_resolution);
int zindex = floor((z0 + space_dimension/2.)/space_resolution);

if (xindex<0) xindex=0;
if (yindex<0) yindex=0;
if (zindex<0) zindex=0;
if (xindex>array_dimension-1) xindex=array_dimension-1;
if (yindex>array_dimension-1) yindex=array_dimension-1;
if (zindex>array_dimension-1) zindex=array_dimension-1;

double xfraction = (x0 + space_dimension/2.)/space_resolution - (double)xindex;
double yfraction = (y0 + space_dimension/2.)/space_resolution - (double)yindex;
double zfraction = (z0 + space_dimension/2.)/space_resolution - (double)zindex;

double B_lower = fieldmap[zindex + yindex*array_dimension + xindex*array_dimension*array_dimension];
double Bx_upper = fieldmap[(zindex) + (yindex)*array_dimension + (xindex+1)*array_dimension*array_dimension];
double By_upper = fieldmap[(zindex) + (yindex+1)*array_dimension + (xindex)*array_dimension*array_dimension];
double Bz_upper = fieldmap[(zindex+1) + (yindex)*array_dimension + (xindex)*array_dimension*array_dimension];

bool positive = (Bx_upper - B_lower)*xfraction + 
                (By_upper - B_lower)*yfraction + 
                (Bz_upper - B_lower)*zfraction > 0.;

double sign = 0.;
if (positive) sign = 1.; else sign = -1.;

double B = B_lower + 
            sign * pow(
              pow((Bx_upper - B_lower)*xfraction,2.) +
              pow((By_upper - B_lower)*yfraction,2.) +
              pow((Bz_upper - B_lower)*zfraction,2.), 0.5);

//printf("Requested z0 is %f and B_lower is %f and B_upper is %f and fractions are %f %f %f, final Bz is %g\n", z0, B_lower, Bz_upper, xfraction, yfraction, zfraction, B);

return B;

}





double *TrappedElectronGenerator::GetBzMap() const
{
double space_dimension = SPACE_DIMENSION;
double space_resolution = SPACE_RESOLUTION;
int array_dimension = (int)(space_dimension/space_resolution);
double *BzMap = new double[array_dimension*array_dimension*array_dimension];
double x0 = 0.;
double y0 = 0.;
double z0 = 0.;

for (int i=0; i<(int)array_dimension; i++)
  {
  printf("Field map in Z, step %d out of %d\n", i, (int)array_dimension-1);
  for (int j=0; j<array_dimension; j++)
    for (int k=0; k<array_dimension; k++)
      {
      x0 = -space_dimension/2. + (double)i * space_resolution;
      y0 = -space_dimension/2. + (double)j * space_resolution;
      z0 = -space_dimension/2. + (double)k * space_resolution;
      BzMap[k + array_dimension*j + array_dimension*array_dimension*i] = GetBz(x0, y0, z0);
      }
  }

return BzMap;

}


double *TrappedElectronGenerator::GetBxMap() const
{
double space_dimension = SPACE_DIMENSION;
double space_resolution = SPACE_RESOLUTION;
int array_dimension = (int)(space_dimension/space_resolution);
double *BxMap = new double[array_dimension*array_dimension*array_dimension];
double x0 = 0.;
double y0 = 0.;
double z0 = 0.;

for (int i=0; i<(int)array_dimension; i++)
  {
  printf("Field map in X, step %d out of %d\n", i, (int)array_dimension-1);
  for (int j=0; j<array_dimension; j++)
    for (int k=0; k<array_dimension; k++)
      {
      x0 = -space_dimension/2. + (double)i * space_resolution;
      y0 = -space_dimension/2. + (double)j * space_resolution;
      z0 = -space_dimension/2. + (double)k * space_resolution;
      BxMap[k + array_dimension*j + array_dimension*array_dimension*i] = GetBx(x0, y0, z0);
      }
  }

return BxMap;

}

double *TrappedElectronGenerator::GetByMap() const
{
double space_dimension = SPACE_DIMENSION;
double space_resolution = SPACE_RESOLUTION;
int array_dimension = (int)(space_dimension/space_resolution);
double *ByMap = new double[array_dimension*array_dimension*array_dimension];
double x0 = 0.;
double y0 = 0.;
double z0 = 0.;

for (int i=0; i<(int)array_dimension; i++)
  {
  printf("Field map in Y, step %d out of %d\n", i, (int)array_dimension-1);
  for (int j=0; j<array_dimension; j++)
    for (int k=0; k<array_dimension; k++)
      {
      x0 = -space_dimension/2. + (double)i * space_resolution;
      y0 = -space_dimension/2. + (double)j * space_resolution;
      z0 = -space_dimension/2. + (double)k * space_resolution;
      ByMap[k + array_dimension*j + array_dimension*array_dimension*i] = GetBy(x0, y0, z0);
      }
  }

return ByMap;

}







void TrappedElectronGenerator::BuildFieldMaps()
{

double n = SPACE_DIMENSION/SPACE_RESOLUTION;

FILE *fp1 = fopen("fieldmap.bin","rb");
if (fp1==NULL) 
{

BzMap = GetBzMap();
BxMap = GetBxMap();
ByMap = GetByMap();

FILE *fp0 = fopen("./fieldmap.bin","wb");
if (fp0==NULL) {fputs ("File error",stderr);}

fwrite(&n, sizeof(n), 1, fp0);
fwrite(&BzMap[0], sizeof(BzMap[0]), (int)(n*n*n), fp0);
fwrite(&BxMap[0], sizeof(BxMap[0]), (int)(n*n*n), fp0);
fwrite(&ByMap[0], sizeof(ByMap[0]), (int)(n*n*n), fp0);

fclose(fp0);

}

else 
{
fread(&n, sizeof(n), 1, fp1);
printf("n read from file is %f\n", n);
BzMap = (double*) malloc (sizeof(double)*(int)(n*n*n));
BxMap = (double*) malloc (sizeof(double)*(int)(n*n*n));
ByMap = (double*) malloc (sizeof(double)*(int)(n*n*n));

fread(BzMap, sizeof(BzMap[0]), (int)(n*n*n), fp1);
fread(BxMap, sizeof(BxMap[0]), (int)(n*n*n), fp1);
fread(ByMap, sizeof(ByMap[0]), (int)(n*n*n), fp1);

fclose(fp1);

}


}



double TrappedElectronGenerator::GetBMagInterpolated(double x0, double y0, double z0) const
{

// interpolate pre-built field map:

double Bmag = pow( pow(InterpolateB(x0,y0,z0,BxMap),2.) + 
                   pow(InterpolateB(x0,y0,z0,ByMap),2.) + 
                   pow(InterpolateB(x0,y0,z0,BzMap),2.), 0.5);

return Bmag;

}









} /* namespace locust */
