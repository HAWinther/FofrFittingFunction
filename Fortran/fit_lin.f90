module FittingFormula
  implicit none
  public :: pofk_enhancement_linear

  !================================================================
  ! 
  ! This is a fitting formula for the enhancement:
  !
  ! pofk_enhancement = P(k,z) / P_LCDM(k,z) 
  !
  ! in Hu-Sawikcy f(R) gravity with n =1
  !
  ! [k] is assumed to in units of [h/Mpc]
  ! 
  ! NB: If requesting a value outside of the regime we have
  ! fitted for we return the closest value without
  ! any warning
  !
  !================================================================

  integer, parameter :: dl = KIND(1.d0)
  real(dl) :: fR0_fid = 1e-5
  real(dl), dimension(1:12) :: param_fid = [ &
    & 3.1       ,  2.34466 , -1.86362  ,  34.4951   , & 
    & 28.8637   , -13.1302 ,  0.146542 , -0.0103024 , &
    & -0.149446 ,  1.62807 ,  0.71291  , -1.41003]
  real(dl),parameter :: kmin = 0.001;
  real(dl),parameter :: kmax = 100.0;
  real(dl),parameter :: zmin = 0.0;
  real(dl),parameter :: zmax = 5.0;
  real(dl),parameter :: fR0min = 1e-7;
  real(dl),parameter :: fR0max = 1e-3;

contains

  function pofk_enhancement_linear(zin,fR0in,kin)
    implicit none
    real(dl) :: zin, fR0in, kin
    real(dl) :: z, fR0, k, a
    real(dl) :: pofk_enhancement_linear
    real(dl) :: aminusone, aminusone2
    real(dl) :: b0,b1,b2,c0,c1,c2,d0,d1,d2,e0,e1,e2
    real(dl) :: b,c,d,e

    z = zin; fR0 = fR0in; k = kin
    
    ! Fit valid for all fR0, but fiducial model only fit to [0.001,100]
    ! so for large or small fR0 we don't have relevant scales so limit it
    if(fR0 < fR0min) fR0 = fR0min;
    if(fR0 > fR0max) fR0 = fR0max;

    ! Transform to effective k
    k = k * sqrt(fR0 / fR0_fid);
 
    ! The ranges we performed the fit for
    if(k > kmax)     k   = kmax;
    if(z < zmin)     z   = zmin;
    if(z > zmax)     z   = zmax;
  
    ! Avoid NaN
    if(abs(k-1) < 1e-10) k = 1.0 + 1e-10;
    if(k < 1e-10) then
      pofk_enhancement_linear = 1.0
      return
    endif

    aminusone  = 1.0/(1.0+z) - 1.0;
    aminusone2 = aminusone*aminusone;
    
    b0 = param_fid(1+ 0) 
    b1 = param_fid(1+ 1)
    b2 = param_fid(1+ 2)
    c0 = param_fid(1+ 3) 
    c1 = param_fid(1+ 4)
    c2 = param_fid(1+ 5)
    d0 = param_fid(1+ 6) 
    d1 = param_fid(1+ 7)
    d2 = param_fid(1+ 8)
    e0 = param_fid(1+ 9)
    e1 = param_fid(1+10)
    e2 = param_fid(1+11)

    b = b0 + b1*aminusone + b2*aminusone2;
    c = c0 + c1*aminusone + c2*aminusone2;
    d = d0 + d1*aminusone + d2*aminusone2;
    e = e0 + e1*aminusone + e2*aminusone2;
  
    pofk_enhancement_linear = 1.0 + (b * k)**2 / (1.0 + c * k * k) + abs(log(k)*k/(k-1)) * d * atan(e * k);

  end function

end module FittingFormula
