MODULE precparcel
use parceldata
implicit none

double precision,parameter :: cpd      = 1005.7           !<    *specific heat at constant pressure (dry air).
double precision,parameter :: cpv      = 1870.            !<    *specific heat at constant pressure (vapor).
double precision,parameter :: cpl      = 4190.            !<    *specific heat at constant pressure (liquid).
double precision :: cpi      = 2106.            !<    *specific heat at constant pressure (ice).
double precision,parameter :: srd      = 6775.            !<    *specific entropy at triple point (dry air).
double precision,parameter :: srv      = 10320.           !<    *specific entropy at triple point (vapor).
double precision,parameter :: srl      = 3517.            !<    *specific entropy at triple point (liquid).
double precision :: sri      = 2296.            !<    *specific entropy at triple point (ice).
double precision,parameter :: hrd      = 530000.          !<    *specific enthalpy at triple point (dry air).
double precision,parameter :: hrv      = 3132927.         !<    *specific enthalpy at triple point (vapor).
double precision,parameter :: hrl      = 632000.          !<    *specific enthalpy at triple point (liquid).
double precision :: hri      = 298484.          !<    *specific enthalpy at triple point (ice).
double precision :: tup      = 273.0             !<    *temperature range over which mixed phase occurs (high)
double precision :: tdn      = 253.0             !<    *temperature range over which mixed phase occurs (low)


CONTAINS

SUBROUTINE entrain(tvsurfmean,tsurfexc,qtsurfexc,rhenv,lapse,eps,ipiso,cmse,iice,isatice,ismooth,ipz,iprecliq,iprecice,plenl,pleni)

! calculates a moist adiabat using several proposed prognostic variables
! takes into account the dependence of latent heat on temperature

! Compilation with gfortran
! f95 -fdefault-double precision-8 simpletest.f90 -o simpletest

! Running with output
!./simpletest > gpfile

! Viewing the results in gnuplot (sorry, did not write nice interactive python wrapper)
! p 'gpfile' u 2:1 w lp, '' u 3:1 w lp, '' u 4:1 w lp,'' u 5:1 w lp,'' u 6:1 w lp,'' u 7:1 w lp,'' u 8:1 w lp
double precision :: tvsurfmean
double precision :: qtsurf,qtsurfexc
double precision :: tsurfexc
double precision :: tsurf
double precision :: rhenv
double precision :: lapse
double precision :: eps,epsrk1,epsrk2
integer,intent(in) :: ipiso,cmse,iice,isatice,ismooth,ipz,iprecliq,iprecice
double precision :: plenl,pleni
double precision :: zsurf=0.
double precision :: psurf=101325.
double precision,dimension(4) :: zmat=(/11000.,20000.,32000.,47000./)
double precision,dimension(4) :: lapserate=(/-6.5/1000.,0.0,1./1000,2.8/1000/) ! standard atmospheric lapse rate below 11000 meters
double precision,dimension(4) :: pmat
double precision,dimension(4) :: tmat
! 
! temporary parameters
double precision :: tnr,tnr_old,ttry,ql,qi,tvguess,tvguessmin
double precision :: qtp,thlp,entp,msep,thsp,tvplast,qlp,qtest,mseaprp,thlaprp,thllinp,theaprp
double precision :: qtprk1,thlprk1,entprk1,mseprk1,thsprk1,mseaprprk1,thlaprprk1,thllinprk1,theaprprk1
double precision :: qtp2,thlp2,entp2,msep2,thsp2,mseaprp2,thlaprp2,thllinp2,theaprp2

double precision :: qtplast,mseplast,thsplast,thlplast,entplast,mseaprplast,thlaprplast,thllinplast,theaprplast
integer :: niter

! surface values and guesses
double precision :: entsurf,entguess,entguessmin
double precision :: thssurf,thsguess,thsguessmin
double precision :: thlsurf,thlguess,thlguessmin
double precision :: thlaprsurf,thlaprguess,thlaprguessmin
double precision :: thllinsurf,thllinguess,thllinguessmin
double precision :: theaprsurf,theaprguess,theaprguessmin
double precision :: msesurf,mseguess,mseguessmin
double precision :: mseaprsurf,mseaprguess,mseaprguessmin
double precision :: tvsurf,qlsurf,qisurf

! dummy loop variables
integer :: i,j

if(iice==0) then
  cpi      = 4190.            !<    *specific heat at constant pressure (liquid).
  sri      = 3517.            !<    *specific entropy at triple point (liquid).
  hri      = 632000.          !<    *specific enthalpy at triple point (liquid).
else
  sri      = 2296.
  cpi      = 2106.            !<    *specific heat at constant pressure (ice).
  hri      = 298484. 
endif

if(ismooth==1) then
  tup=268.0
  tdn=253.0
else
  tup=268.0
  tdn=267.8
endif

tsurf=tvsurfmean

! use standard atmospheric lapse rate for background pressure, any pressure profile should work though
pmat(1)=exp((log(psurf)*lapse*rd+log(tsurf+zsurf*lapse)*grav-&
log(tsurf+zmat(1)*lapse)*grav)/(lapse*rd))
tmat(1)=tsurf+lapse*(zmat(1)-zsurf);
! write(*,*)(*,*) 'make profiles'

do i=2,4
  if(abs(lapserate(i))<1e-10) then
    pmat(i)=exp((log(pmat(i-1))*tmat(i-1)*rd+zmat(i-1)*grav-zmat(i)*grav)/(tmat(i-1)*rd))
  else
    pmat(i)=exp((log(pmat(i-1))*lapserate(i)*rd+log(tmat(i-1)+zmat(i-1)*lapserate(i))*grav-&
    log(tmat(i-1)+zmat(i)*lapserate(i))*grav)/(lapserate(i)*rd))
  endif
  tmat(i)=tmat(i-1)+lapserate(i)*(zmat(i)-zmat(i-1));
enddo

do i=0,ztop
  z(i)=dz*i
  if(z(i)<=zmat(1)) then
    penv(i)=exp((log(psurf)*lapse*rd+log(tsurf+zsurf*lapse)*grav-&
    log(tsurf+z(i)*lapse)*grav)/(lapse*rd))
    tvenv(i)=tsurf+z(i)*lapse
  else
    j=0
    do while(z(i)>zmat(j))
      j=j+1
    end do
     tvenv(i)=tmat(j-1)+(z(i)-zmat(j-1))*lapserate(j)
     if(abs(lapserate(j))<1e-10) then
       penv(i)=exp((log(pmat(j-1))*tmat(j-1)*rd-(z(i)-zmat(j-1))*grav)/(tmat(j-1)*rd))
     else
       penv(i)=exp((log(pmat(j-1))*lapserate(j)*rd+log(tmat(j-1)+zmat(j-1)*lapserate(j))*grav-&
       log(tmat(j-1)+z(i)*lapserate(j))*grav)/(lapserate(j)*rd))     
     endif
  endif
enddo

!using relative humity, calculate the environments properties
do i=0,ztop
  !first guess for temperature...dry
  if(i>0) then
    tnr=tvenv(i)/(1.0+(rv/rd)*qtenv(i-1))
  else
    tnr=tvenv(0)
  endif
  tnr_old=0.
  niter=0
  do while ((abs(tnr-tnr_old) > 0.0002).and.(niter<99).and.(tnr>100).and.(tnr<700))
    call qsat(isatice,tnr,penv(i),qtest)
    tvguess=tnr*(1.0+(rv/rd-1.0)*rhenv*qtest)
    ttry=tnr-0.0002
    call qsat(isatice,ttry,penv(i),qtest)
    tvguessmin=ttry*(1.0+(rv/rd-1.0)*rhenv*qtest)
    tnr = tnr - (tvguess-tvenv(i))/((tvguess-tvguessmin)*5000.)
    niter=niter+1
  enddo
  call qsat(isatice,tnr,penv(i),qtest)
  qtp=rhenv*qtest
  tenv(i)=tnr
  qtenv(i)=qtp
  entenv(i) = ent(tnr,penv(i),qtp,ql,qi)
  thlenv(i) = thl(tnr,penv(i),qtp,ql,qi)
  thsenv(i) = ths(tnr,penv(i),qtp,ql,qi)
  mseenv(i) = mse(tnr,z(i),qtp,ql,qi)
  mseaprenv(i) = mseapr(tnr,z(i),qtp,ql,qi)
  thlaprenv(i) = thlapr(tnr,penv(i),qtp,ql,qi)
  thllinenv(i) = thllin(tnr,penv(i),qtp,ql,qi)
  theaprenv(i) = theapr(tnr,penv(i),qtp,ql,qi)
  qlp=ql+qi
enddo

tsurf=tenv(0)+tsurfexc
qtsurf=max(qtenv(0)+qtsurfexc,0.)

call qli(isatice,tsurf,psurf,qtsurf,qlsurf,qisurf,tvsurf)
!surface thermo calculations
entsurf=ent(tsurf,psurf,qtsurf,qlsurf,qisurf)
thssurf=ths(tsurf,psurf,qtsurf,qlsurf,qisurf)
thlsurf=thl(tsurf,psurf,qtsurf,qlsurf,qisurf)
msesurf=mse(tsurf,zsurf,qtsurf,qlsurf,qisurf)
mseaprsurf=mseapr(tsurf,zsurf,qtsurf,qlsurf,qisurf)
thlaprsurf=thlapr(tsurf,psurf,qtsurf,qlsurf,qisurf)
thllinsurf=thllin(tsurf,psurf,qtsurf,qlsurf,qisurf)
theaprsurf=theapr(tsurf,psurf,qtsurf,qlsurf,qisurf)
piso(0)=psurf


qtp=qtsurf
qtplast=qtsurf
entp=entsurf
entplast=entsurf
qlp=qlsurf+qisurf
tvplast=tvsurf
tnr=tsurf

do i=1,ztop
  !first guess for virtual temperature at level n+1...dry
  tvguess=tvplast+9.81/((1.0-qtp)*cpd+qtp*cpv)
  do j=1,5
    piso(i)= exp(log(piso(i-1))-(grav/(0.5*rd*(tvplast+tvguess)))*(z(i)-z(i-1)))
    if(ipiso==1) then
      pp(i)=piso(i)
    else
      pp(i)=penv(i)
    endif
    if(qlp>0.00001 .and. z(i)>200) then
      if(ipz==0) then
        epsrk1=0.001*eps
        epsrk2=0.001*eps
      else
        epsrk1=eps/(z(i-1))
        epsrk2=2*eps/(z(i-1)+z(i))
      endif
      qtprk1=-epsrk1*(qtplast-qtenv(i-1))*(z(i)-z(i-1))
      qtp2=-epsrk2*(qtplast+qtprk1/2.0-(qtenv(i-1)+qtenv(i))/2.0)*(z(i)-z(i-1))
      qtp=qtplast+qtp2
      entprk1=-epsrk1*(entplast-entenv(i-1))*(z(i)-z(i-1))
      entp2=-epsrk2*(entplast+entprk1/2.0-(entenv(i-1)+entenv(i))/2.0)*(z(i)-z(i-1))
      entp=entplast+entp2
    else
      qtp=qtplast
      entp=entplast
    endif
    if(qlp<0.00001) then
      tnr=tnr-0.00981!tref*exp((entp &
      !+(1-qtp)*rd*log(pp(i)/(pref*(1+((qtp)*rv)/((1-qtp)*rd)))) &
      !+ (qtp)*rv*log(pp(i)/(pref*(1+((1-qtp)*rd)/(qtp*rv)))) &
      !- srd - qtp*(srv-srd))/((1+qtp*(cpv-cpd)/(cpd))*cpd))
    else
      tnr=tnr-0.007
    endif
    tnr_old=0.
    niter=0
    do while ((abs(tnr-tnr_old) > 0.0002).and.(niter<99).and.(tnr>100).and.(tnr<700))
      niter = niter+1
      tnr_old=tnr
      call qli(isatice,tnr,pp(i),qtp,ql,qi,tvguess)
      entguess = ent(tnr,pp(i),qtp,ql,qi)
      ttry=tnr-0.0002
      call qli(isatice,ttry,pp(i),qtp,ql,qi,tvguessmin)
      entguessmin = ent(ttry,pp(i),qtp,ql,qi)
      tnr = tnr - (entguess-entp)/((entguess-entguessmin)*5000.)
    enddo
  enddo
  tent(i)=tnr
  if(iprecliq==1) then
    qtp=qtp-min(ql,ql*(z(i)-z(i-1))/plenl)
    ql=ql-min(ql,ql*(z(i)-z(i-1))/plenl)
  endif
  if(iprecice==1) then
    qtp=qtp-min(qi,qi*(z(i)-z(i-1))/pleni)
    qi=qi-min(qi,qi*(z(i)-z(i-1))/pleni)
  endif
  entp=ent(tnr,pp(i),qtp,ql,qi)
  ! integration values
  call qli(isatice,tnr,pp(i),qtp,ql,qi,tvguess)
  tvplast=tvguess
  tvent(i)=tvguess
  qlp=ql+qi
  entplast=entp
  qtplast=qtp
enddo

qtp=qtsurf
qtplast=qtsurf
msep=msesurf
mseplast=msesurf
qlp=qlsurf+qisurf
tvplast=tvsurf
tnr=tsurf

do i=1,ztop
  tvguess=tvplast+9.81/((1.0-qtp)*cpd+qtp*cpv)
  do j=1,5
    piso(i)= exp(log(piso(i-1))-(grav/(0.5*rd*(tvplast+tvguess)))*(z(i)-z(i-1)))
    if(ipiso==1) then
      pp(i)=piso(i)
    else
      pp(i)=penv(i)
    endif
    if(qlp>0.00001 .and. z(i)>200) then
      if(ipz==0) then
        epsrk1=0.001*eps
        epsrk2=0.001*eps
      else
        epsrk1=eps/(z(i-1))
        epsrk2=2*eps/(z(i-1)+z(i))
      endif
      qtprk1=-epsrk1*(qtplast-qtenv(i-1))*(z(i)-z(i-1))
      qtp2=-epsrk2*(qtplast+qtprk1/2.0-(qtenv(i-1)+qtenv(i))/2.0)*(z(i)-z(i-1))
      qtp=qtplast+qtp2
      mseprk1=-epsrk1*(mseplast-mseenv(i-1))*(z(i)-z(i-1))
      msep2=-epsrk2*(mseplast+mseprk1/2.0-(mseenv(i-1)+mseenv(i))/2.0)*(z(i)-z(i-1))
      msep=mseplast+msep2
    else
      qtp=qtplast
      msep=mseplast
    endif
    if(ipiso==0 .and. cmse==1) then
      msep=msep-grav*(tvguess+tvplast-tvenv(i)-tvenv(i-1))/(tvenv(i)+tvenv(i-1))*(z(i)-z(i-1))
    else
      msep=msep
    endif
    if(qlp<0.0001) then
      tnr=(hrd-cpd*tref-qtp*hrd+qtp*cpd*tref-msep+qtp*hrv-qtp*cpv*tref&
      +grav*z(i))/(-cpd+cpd*qtp-cpv*qtp)
    else
      tnr=tnr-0.007
    endif
    tnr_old=0.
    niter=0
    do while ((abs(tnr-tnr_old) > 0.0002).and.(niter<99).and.(tnr>100).and.(tnr<700))
      niter = niter+1
      tnr_old=tnr
      call qli(isatice,tnr,pp(i),qtp,ql,qi,tvguess)
      mseguess = mse(tnr,z(i),qtp,ql,qi)
      ttry=tnr-0.0002
      call qli(isatice,ttry,pp(i),qtp,ql,qi,tvguessmin)
      mseguessmin = mse(ttry,z(i),qtp,ql,qi)
      tnr = tnr - (mseguess-msep)/((mseguess-mseguessmin)*5000.)
    enddo
  enddo
  tmse(i)=tnr
  if(iprecliq==1) then
    qtp=qtp-min(ql,ql*(z(i)-z(i-1))/plenl)
    ql=ql-min(ql,ql*(z(i)-z(i-1))/plenl)
  endif
  if(iprecice==1) then
    qtp=qtp-min(qi,qi*(z(i)-z(i-1))/pleni)
    qi=qi-min(qi,qi*(z(i)-z(i-1))/pleni)
  endif
  msep=mse(tnr,z(i),qtp,ql,qi)
  ! integration values
  call qli(isatice,tnr,pp(i),qtp,ql,qi,tvguess)
  tvplast=tvguess
  tvmse(i)=tvguess
  qlp=ql+qi
  mseplast=msep
  qtplast=qtp
enddo

qtp=qtsurf
qtplast=qtsurf
mseaprp=mseaprsurf
mseaprplast=mseaprsurf
qlp=qlsurf+qisurf
tvplast=tvsurf
tnr=tsurf

do i=1,ztop
  tvguess=tvplast+9.81/((1.0-qtp)*cpd+qtp*cpv)
  do j=1,5
    piso(i)= exp(log(piso(i-1))-(grav/(0.5*rd*(tvplast+tvguess)))*(z(i)-z(i-1)))
    if(ipiso==1) then
      pp(i)=piso(i)
    else
      pp(i)=penv(i)
    endif
    if(qlp>0.00001 .and. z(i)>200) then
      if(ipz==0) then
        epsrk1=0.001*eps
        epsrk2=0.001*eps
      else
        epsrk1=eps/(z(i-1))
        epsrk2=2*eps/(z(i-1)+z(i))
      endif
      qtprk1=-epsrk1*(qtplast-qtenv(i-1))*(z(i)-z(i-1))
      qtp2=-epsrk2*(qtplast+qtprk1/2.0-(qtenv(i-1)+qtenv(i))/2.0)*(z(i)-z(i-1))
      qtp=qtplast+qtp2
      mseaprprk1=-epsrk1*(mseaprplast-mseaprenv(i-1))*(z(i)-z(i-1))
      mseaprp2=-epsrk2*(mseaprplast+mseaprprk1/2.0-(mseaprenv(i-1)+mseaprenv(i))/2.0)*(z(i)-z(i-1))
      mseaprp=mseaprplast+mseaprp2
    else
      qtp=qtplast
      mseaprp=mseaprplast
    endif
    if(ipiso==0 .and. cmse==1) then
      mseaprp=mseaprp-grav*(tvguess+tvplast-tvenv(i)-tvenv(i-1))/(tvenv(i)+tvenv(i-1))*(z(i)-z(i-1))
    else
      mseaprp=mseaprp
    endif
    if(qlp<0.0001) then
      tnr=(hrd-cpd*tref-qtp*hrd+qtp*cpd*tref-mseaprp+qtp*hrv-qtp*cpv*tref&
      +grav*z(i))/(-cpd+cpd*qtp-cpv*qtp)
    else
      tnr=tnr-0.007
    endif
    tnr_old=0.
    niter=0
    do while ((abs(tnr-tnr_old) > 0.0002).and.(niter<99).and.(tnr>100).and.(tnr<700))
      niter = niter+1
      tnr_old=tnr
      call qli(isatice,tnr,pp(i),qtp,ql,qi,tvguess)
      mseaprguess = mseapr(tnr,z(i),qtp,ql,qi)
      ttry=tnr-0.0002
      call qli(isatice,ttry,pp(i),qtp,ql,qi,tvguessmin)
      mseaprguessmin = mseapr(ttry,z(i),qtp,ql,qi)
      tnr = tnr - (mseaprguess-mseaprp)/((mseaprguess-mseaprguessmin)*5000.)
    enddo
  enddo
  tmseapr(i)=tnr
  if(iprecliq==1) then
    qtp=qtp-min(ql,ql*(z(i)-z(i-1))/plenl)
    ql=ql-min(ql,ql*(z(i)-z(i-1))/plenl)
  endif
  if(iprecice==1) then
    qtp=qtp-min(qi,qi*(z(i)-z(i-1))/pleni)
    qi=qi-min(qi,qi*(z(i)-z(i-1))/pleni)
  endif
  mseaprp=mseapr(tnr,z(i),qtp,ql,qi)
  ! integration values
  call qli(isatice,tnr,pp(i),qtp,ql,qi,tvguess)
  tvplast=tvguess
  tvmseapr(i)=tvguess
  qlp=ql+qi
  mseaprplast=mseaprp
  qtplast=qtp
enddo

qtp=qtsurf
qtplast=qtsurf
thsp=thssurf
thsplast=thssurf
qlp=qlsurf+qisurf
tvplast=tvsurf
tnr=tsurf
msep=msesurf

do i=1,ztop
  tvguess=tvplast+9.81/((1.0-qtp)*cpd+qtp*cpv)
  do j=1,5
    piso(i)= exp(log(piso(i-1))-(grav/(0.5*rd*(tvplast+tvguess)))*(z(i)-z(i-1)))
    if(ipiso==1) then
      pp(i)=piso(i)
    else
      pp(i)=penv(i)
    endif
    if(qlp>0.00001 .and. z(i)>200) then
      if(ipz==0) then
        epsrk1=0.001*eps
        epsrk2=0.001*eps
      else
        epsrk1=eps/(z(i-1))
        epsrk2=2*eps/(z(i-1)+z(i))
      endif
      qtprk1=-epsrk1*(qtplast-qtenv(i-1))*(z(i)-z(i-1))
      qtp2=-epsrk2*(qtplast+qtprk1/2.0-(qtenv(i-1)+qtenv(i))/2.0)*(z(i)-z(i-1))
      qtp=qtplast+qtp2
      thsprk1=-epsrk1*(thsplast-thsenv(i-1))*(z(i)-z(i-1))
      thsp2=-epsrk2*(thsplast+thsprk1/2.0-(thsenv(i-1)+thsenv(i))/2.0)*(z(i)-z(i-1))
      thsp=thsplast+thsp2
    else
      qtp=qtplast
      thsp=thsplast
    endif
    if(qlp<0.0001) then
      tnr=(hrd-cpd*tref-qtp*hrd+qtp*cpd*tref-msep+qtp*hrv-qtp*cpv*tref&
      +grav*z(i))/(-cpd+cpd*qtp-cpv*qtp)
    else
      tnr=tnr-0.007
    endif
    tnr_old=0.
    niter=0
    do while ((abs(tnr-tnr_old) > 0.0002).and.(niter<99).and.(tnr>100).and.(tnr<700))
      niter = niter+1
      tnr_old=tnr
      call qli(isatice,tnr,pp(i),qtp,ql,qi,tvguess)
      thsguess = ths(tnr,pp(i),qtp,ql,qi)
      ttry=tnr-0.0002
      call qli(isatice,ttry,pp(i),qtp,ql,qi,tvguessmin)
      thsguessmin = ths(ttry,pp(i),qtp,ql,qi)
      tnr = tnr - (thsguess-thsp)/((thsguess-thsguessmin)*5000.)
    enddo
  enddo
  tths(i)=tnr
  if(iprecliq==1) then
    qtp=qtp-min(ql,ql*(z(i)-z(i-1))/plenl)
    ql=ql-min(ql,ql*(z(i)-z(i-1))/plenl)
  endif
  if(iprecice==1) then
    qtp=qtp-min(qi,qi*(z(i)-z(i-1))/pleni)
    qi=qi-min(qi,qi*(z(i)-z(i-1))/pleni)
  endif
  thsp=ths(tnr,pp(i),qtp,ql,qi)
  ! integration values
  call qli(isatice,tnr,pp(i),qtp,ql,qi,tvguess)
  tvplast=tvguess
  tvths(i)=tvguess
  qlp=ql+qi
  thsplast=thsp
  qtplast=qtp
enddo

qtp=qtsurf
qtplast=qtsurf
thlp=thlsurf
thlplast=thlsurf
qlp=qlsurf+qisurf
tvplast=tvsurf
tnr=tsurf
msep=msesurf

do i=1,ztop
  tvguess=tvplast+9.81/((1.0-qtp)*cpd+qtp*cpv)
  do j=1,5
    piso(i)= exp(log(piso(i-1))-(grav/(0.5*rd*(tvplast+tvguess)))*(z(i)-z(i-1)))
    if(ipiso==1) then
      pp(i)=piso(i)
    else
      pp(i)=penv(i)
    endif
    if(qlp>0.00001 .and. z(i)>200) then
      if(ipz==0) then
        epsrk1=0.001*eps
        epsrk2=0.001*eps
      else
        epsrk1=eps/(z(i-1))
        epsrk2=2*eps/(z(i-1)+z(i))
      endif
      qtprk1=-epsrk1*(qtplast-qtenv(i-1))*(z(i)-z(i-1))
      qtp2=-epsrk2*(qtplast+qtprk1/2.0-(qtenv(i-1)+qtenv(i))/2.0)*(z(i)-z(i-1))
      qtp=qtplast+qtp2
      thlprk1=-epsrk1*(thlplast-thlenv(i-1))*(z(i)-z(i-1))
      thlp2=-epsrk2*(thlplast+thlprk1/2.0-(thlenv(i-1)+thlenv(i))/2.0)*(z(i)-z(i-1))
      thlp=thlplast+thlp2
    else
      qtp=qtplast
      thlp=thlplast
    endif
    if(qlp<0.0001) then
      tnr=(hrd-cpd*tref-qtp*hrd+qtp*cpd*tref-msep+qtp*hrv-qtp*cpv*tref&
      +grav*z(i))/(-cpd+cpd*qtp-cpv*qtp)
    else
      tnr=tnr-0.007
    endif
    tnr_old=0.
    niter=0
    do while ((abs(tnr-tnr_old) > 0.0002).and.(niter<99).and.(tnr>100).and.(tnr<700))
      niter = niter+1
      tnr_old=tnr
      call qli(isatice,tnr,pp(i),qtp,ql,qi,tvguess)
      thlguess = thl(tnr,pp(i),qtp,ql,qi)
      ttry=tnr-0.0002
      call qli(isatice,ttry,pp(i),qtp,ql,qi,tvguessmin)
      thlguessmin = thl(ttry,pp(i),qtp,ql,qi)
      tnr = tnr - (thlguess-thlp)/((thlguess-thlguessmin)*5000.)
    enddo
  enddo
  tthl(i)=tnr
  if(iprecliq==1) then
    qtp=qtp-min(ql,ql*(z(i)-z(i-1))/plenl)
    ql=ql-min(ql,ql*(z(i)-z(i-1))/plenl)
  endif
  if(iprecice==1) then
    qtp=qtp-min(qi,qi*(z(i)-z(i-1))/pleni)
    qi=qi-min(qi,qi*(z(i)-z(i-1))/pleni)
  endif
  thlp=thl(tnr,pp(i),qtp,ql,qi)
  ! integration values
  call qli(isatice,tnr,pp(i),qtp,ql,qi,tvguess)
  tvplast=tvguess
  tvthl(i)=tvguess
  qlp=ql+qi
  thlplast=thlp
  qtplast=qtp
  !write(*,*) qtplast
enddo

qtp=qtsurf
qtplast=qtsurf
thlaprp=thlaprsurf
thlaprplast=thlaprsurf
qlp=qlsurf+qisurf
tvplast=tvsurf
tnr=tsurf
msep=msesurf

do i=1,ztop
  tvguess=tvplast+9.81/((1.0-qtp)*cpd+qtp*cpv)
  do j=1,5
    piso(i)= exp(log(piso(i-1))-(grav/(0.5*rd*(tvplast+tvguess)))*(z(i)-z(i-1)))
    if(ipiso==1) then
      pp(i)=piso(i)
    else
      pp(i)=penv(i)
    endif
    if(qlp>0.00001 .and. z(i)>200) then
      if(ipz==0) then
        epsrk1=0.001*eps
        epsrk2=0.001*eps
      else
        epsrk1=eps/(z(i-1))
        epsrk2=2*eps/(z(i-1)+z(i))
      endif
      qtprk1=-epsrk1*(qtplast-qtenv(i-1))*(z(i)-z(i-1))
      qtp2=-epsrk2*(qtplast+qtprk1/2.0-(qtenv(i-1)+qtenv(i))/2.0)*(z(i)-z(i-1))
      qtp=qtplast+qtp2
      thlaprprk1=-epsrk1*(thlaprplast-thlaprenv(i-1))*(z(i)-z(i-1))
      thlaprp2=-epsrk2*(thlaprplast+thlaprprk1/2.0-(thlaprenv(i-1)+thlaprenv(i))/2.0)*(z(i)-z(i-1))
      thlaprp=thlaprplast+thlaprp2
    else
      qtp=qtplast
      thlaprp=thlaprplast
    endif
    if(qlp<0.0001) then
      tnr=(hrd-cpd*tref-qtp*hrd+qtp*cpd*tref-msep+qtp*hrv-qtp*cpv*tref&
      +grav*z(i))/(-cpd+cpd*qtp-cpv*qtp)
    else
      tnr=tnr-0.007
    endif
    tnr_old=0.
    niter=0
    do while ((abs(tnr-tnr_old) > 0.0002).and.(niter<99).and.(tnr>100).and.(tnr<700))
      niter = niter+1
      tnr_old=tnr
      call qli(isatice,tnr,pp(i),qtp,ql,qi,tvguess)
      thlaprguess = thlapr(tnr,pp(i),qtp,ql,qi)
      ttry=tnr-0.0002
      call qli(isatice,ttry,pp(i),qtp,ql,qi,tvguessmin)
      thlaprguessmin = thlapr(ttry,pp(i),qtp,ql,qi)
      tnr = tnr - (thlaprguess-thlaprp)/((thlaprguess-thlaprguessmin)*5000.)
    enddo
  enddo
  tthlapr(i)=tnr
  if(iprecliq==1) then
    qtp=qtp-min(ql,ql*(z(i)-z(i-1))/plenl)
    ql=ql-min(ql,ql*(z(i)-z(i-1))/plenl)
  endif
  if(iprecice==1) then
    qtp=qtp-min(qi,qi*(z(i)-z(i-1))/pleni)
    qi=qi-min(qi,qi*(z(i)-z(i-1))/pleni)
  endif
  thlaprp=thlapr(tnr,pp(i),qtp,ql,qi)
  ! integration values
  call qli(isatice,tnr,pp(i),qtp,ql,qi,tvguess)
  tvplast=tvguess
  tvthlapr(i)=tvguess
  qlp=ql+qi
  thlaprplast=thlaprp
  qtplast=qtp
  !write(*,*) qtplast
enddo

qtp=qtsurf
qtplast=qtsurf
thllinp=thllinsurf
thllinplast=thllinsurf
qlp=qlsurf+qisurf
tvplast=tvsurf
tnr=tsurf
msep=msesurf

do i=1,ztop
  tvguess=tvplast+9.81/((1.0-qtp)*cpd+qtp*cpv)
  do j=1,5
    piso(i)= exp(log(piso(i-1))-(grav/(0.5*rd*(tvplast+tvguess)))*(z(i)-z(i-1)))
    if(ipiso==1) then
      pp(i)=piso(i)
    else
      pp(i)=penv(i)
    endif
    if(qlp>0.00001 .and. z(i)>200) then
      if(ipz==0) then
        epsrk1=0.001*eps
        epsrk2=0.001*eps
      else
        epsrk1=eps/(z(i-1))
        epsrk2=2*eps/(z(i-1)+z(i))
      endif
      qtprk1=-epsrk1*(qtplast-qtenv(i-1))*(z(i)-z(i-1))
      qtp2=-epsrk2*(qtplast+qtprk1/2.0-(qtenv(i-1)+qtenv(i))/2.0)*(z(i)-z(i-1))
      qtp=qtplast+qtp2
      thllinprk1=-epsrk1*(thllinplast-thllinenv(i-1))*(z(i)-z(i-1))
      thllinp2=-epsrk2*(thllinplast+thllinprk1/2.0-(thllinenv(i-1)+thllinenv(i))/2.0)*(z(i)-z(i-1))
      thllinp=thllinplast+thllinp2
    else
      qtp=qtplast
      thllinp=thllinplast
    endif
    if(qlp<0.0001) then
      tnr=(hrd-cpd*tref-qtp*hrd+qtp*cpd*tref-msep+qtp*hrv-qtp*cpv*tref&
      +grav*z(i))/(-cpd+cpd*qtp-cpv*qtp)
    else
      tnr=tnr-0.007
    endif
    tnr_old=0.
    niter=0
    do while ((abs(tnr-tnr_old) > 0.0002).and.(niter<99).and.(tnr>100).and.(tnr<700))
      niter = niter+1
      tnr_old=tnr
      call qli(isatice,tnr,pp(i),qtp,ql,qi,tvguess)
      thllinguess = thllin(tnr,pp(i),qtp,ql,qi)
      ttry=tnr-0.0002
      call qli(isatice,ttry,pp(i),qtp,ql,qi,tvguessmin)
      thllinguessmin = thllin(ttry,pp(i),qtp,ql,qi)
      tnr = tnr - (thllinguess-thllinp)/((thllinguess-thllinguessmin)*5000.)
    enddo
  enddo
  tthllin(i)=tnr
  if(iprecliq==1) then
    qtp=qtp-min(ql,ql*(z(i)-z(i-1))/plenl)
    ql=ql-min(ql,ql*(z(i)-z(i-1))/plenl)
  endif
  if(iprecice==1) then
    qtp=qtp-min(qi,qi*(z(i)-z(i-1))/pleni)
    qi=qi-min(qi,qi*(z(i)-z(i-1))/pleni)
  endif
  thllinp=thllin(tnr,pp(i),qtp,ql,qi)
  ! integration values
  call qli(isatice,tnr,pp(i),qtp,ql,qi,tvguess)
  tvplast=tvguess
  tvthllin(i)=tvguess
  qlp=ql+qi
  thllinplast=thllinp
  qtplast=qtp
  !write(*,*) qtplast
enddo

qtp=qtsurf
qtplast=qtsurf
theaprp=theaprsurf
theaprplast=theaprsurf
qlp=qlsurf+qisurf
tvplast=tvsurf
tnr=tsurf
msep=msesurf

do i=1,ztop
  tvguess=tvplast+9.81/((1.0-qtp)*cpd+qtp*cpv)
  do j=1,5
    if(ipiso==1) then
      pp(i)=piso(i)
    else
      pp(i)=penv(i)
    endif
    if(qlp>0.00001 .and. z(i)>200) then
      if(ipz==0) then
        epsrk1=0.001*eps
        epsrk2=0.001*eps
      else
        epsrk1=eps/(z(i-1))
        epsrk2=2*eps/(z(i-1)+z(i))
      endif
      qtprk1=-epsrk1*(qtplast-qtenv(i-1))*(z(i)-z(i-1))
      qtp2=-epsrk2*(qtplast+qtprk1/2.0-(qtenv(i-1)+qtenv(i))/2.0)*(z(i)-z(i-1))
      qtp=qtplast+qtp2
      theaprprk1=-epsrk1*(theaprplast-theaprenv(i-1))*(z(i)-z(i-1))
      theaprp2=-epsrk2*(theaprplast+theaprprk1/2.0-(theaprenv(i-1)+theaprenv(i))/2.0)*(z(i)-z(i-1))
      theaprp=theaprplast+theaprp2
    else
      qtp=qtplast
      theaprp=theaprplast
    endif
    if(qlp<0.0001) then
      tnr=(hrd-cpd*tref-qtp*hrd+qtp*cpd*tref-msep+qtp*hrv-qtp*cpv*tref&
      +grav*z(i))/(-cpd+cpd*qtp-cpv*qtp)
    else
      tnr=tnr-0.007
    endif
    tnr_old=0.
    niter=0
    do while ((abs(tnr-tnr_old) > 0.0002).and.(niter<99).and.(tnr>100).and.(tnr<700))
      niter = niter+1
      tnr_old=tnr
      call qli(isatice,tnr,pp(i),qtp,ql,qi,tvguess)
      theaprguess = theapr(tnr,pp(i),qtp,ql,qi)
      ttry=tnr-0.0002
      call qli(isatice,ttry,pp(i),qtp,ql,qi,tvguessmin)
      theaprguessmin = theapr(ttry,pp(i),qtp,ql,qi)
      tnr = tnr - (theaprguess-theaprp)/((theaprguess-theaprguessmin)*5000.)
    enddo
  enddo
  ttheapr(i)=tnr
  if(iprecliq==1) then
    qtp=qtp-min(ql,ql*(z(i)-z(i-1))/plenl)
    ql=ql-min(ql,ql*(z(i)-z(i-1))/plenl)
  endif
  if(iprecice==1) then
    qtp=qtp-min(qi,qi*(z(i)-z(i-1))/pleni)
    qi=qi-min(qi,qi*(z(i)-z(i-1))/pleni)
  endif
  theaprp=theapr(tnr,pp(i),qtp,ql,qi)
  ! integration values
  call qli(isatice,tnr,pp(i),qtp,ql,qi,tvguess)
  tvplast=tvguess
  tvtheapr(i)=tvguess
  qlp=ql+qi
  theaprplast=theaprp
  qtplast=qtp
  !write(*,*) qtplast
enddo

END SUBROUTINE

SUBROUTINE qli(isatice,tt,pp,qt,ql,qi,tv)
  integer, intent(in) :: isatice
  double precision,intent(in) :: tt,pp,qt
  double precision,intent(out) :: ql,qi,tv  
  double precision :: qlisum,qsatur,qvsl,qvsi,ilratio,esattemp
  ilratio=max(0.0,min(1.0,(tt-tdn)/(tup-tdn)))
!   if(iice==1) then
!     esattemp=ilratio*esatl(tt)+(1.0-ilratio)*esati(tt)
!   else
!     esattemp=esatl(tt)
!   endif
!   qsatur=rd/rv*esattemp/(pp-(1.0-rd/rv)*esattemp)
  qvsl=rd/rv*esatl(tt)/(pp-(1.0-rd/rv)*esatl(tt))
  if(isatice==1) then
    qvsi=rd/rv*esati(tt)/(pp-(1.0-rd/rv)*esati(tt))
  else
    qvsi=rd/rv*esatl(tt)/(pp-(1.0-rd/rv)*esatl(tt))
  endif
  qsatur = ilratio*qvsl+(1.0-ilratio)*qvsi
  qlisum=max(qt-qsatur,0.0)
  ql = ilratio*qlisum
  qi = qlisum-ql
  tv = tt*(1.0+(rv/rd-1.0)*(qt-(qlisum)))*(1.0-qlisum)
END SUBROUTINE

DOUBLE PRECISION FUNCTION esatl(tt)
  double precision, intent(in) :: tt
  !esatl=exp(54.842763-6763.22/tt-4.21*log(tt)+0.000367*tt+&
  !tanh(0.0415*(tt-218.8))*(53.878-1331.22/tt-9.44523*log(tt)+ 0.014025*tt))
  esatl=610.7*exp( (1.0/tref-1.0/tt)*((hrv-hrl)-(cpv-cpl)*tref)/rv +&
  ((cpv-cpl)/rv)*log(tt/tref) )
END FUNCTION

DOUBLE PRECISION FUNCTION esati(tt)
  double precision, intent(in) :: tt
  !esati=exp(9.550426-5723.265/tt+3.53068*log(tt)-0.00728332*tt)
  esati=610.7*exp( (1.0/tref-1.0/tt)*((hrv-hri)-(cpv-cpi)*tref)/rv +&
  ((cpv-cpi)/rv)*log(tt/tref) )
END FUNCTION

DOUBLE PRECISION FUNCTION ent(tt,pp,qt,ql,qi)
  double precision, intent(in) :: tt,pp,qt,ql,qi
  ent = (1+qt*(cpv-cpd)/(cpd)+ql*(cpl-cpv)/(cpd)+qi*(cpi-cpv)/(cpd))*cpd*log(tt/tref) &
              -(1-qt)*rd*log(pp/(pref*(1+((qt-ql-qi)*rv)/((1-qt)*rd)))) &
              -(qt-ql-qi)*rv*log(pp/(pref*(1+((1-qt)*rd)/((qt-ql-qi)*rv)))) &
              + srd + qt*(srv-srd) - ql*(srv-srl) - qi*(srv-sri)
END FUNCTION

DOUBLE PRECISION FUNCTION ths(tt,pp,qt,ql,qi)
  double precision, intent(in) :: tt,pp,qt,ql,qi
  if(ql+qi<0.999*qt) then
  ths=tt*&
  (tt/tref)**(qt*(cpv-cpd)/cpd+ql*(cpl-cpv)/cpd+qi*(cpi-cpv)/cpd)*&
  (pp/(pref*(1+((qt-ql-qi)*rv)/((1-qt)*rd))))**(-(1-qt)*rd/cpd)*&
  (pp/(pref*(1+((1-qt)*rd)/((qt-ql-qi)*rv))))**(-(qt-ql-qi)*rv/cpd)*&
  exp(-ql*(srv-srl)/cpd-qi*(srv-sri)/cpd)
  else
  ths=tt*&
  (tt/tref)**(qt*(cpv-cpd)/cpd+ql*(cpl-cpv)/cpd+qi*(cpi-cpv)/cpd)*&
  (pp/(pref*(1+((qt-ql-qi)*rv)/((1-qt)*rd))))**(-(1-qt)*rd/cpd)*&
  exp(-ql*(srv-srl)/cpd-qi*(srv-sri)/cpd)
  endif
END FUNCTION

DOUBLE PRECISION FUNCTION thl(tt,pp,qt,ql,qi)
  double precision, intent(in) :: tt,pp,qt,ql,qi
  if(ql+qi<0.999*qt) then
  thl=tt*&
  (tt/tref)**(qt*(cpv-cpd)/cpd+ql*(cpl-cpv)/cpd+qi*(cpi-cpv)/cpd)*&
  (pp/(pref*(1+((qt-ql-qi)*rv)/((1-qt)*rd))))**(-(1-qt)*rd/cpd)*&
  (pp/(pref*(1+((1-qt)*rd)/((qt-ql-qi)*rv))))**(-(qt-ql-qi)*rv/cpd)*&
  exp(-ql*(srv-srl)/cpd-qi*(srv-sri)/cpd)*(1./(1.+(qt*rv)/((1.-qt)*rd)))**((1.-qt)*rd/cpd)*&
  (1./(1.+((1.-qt)*rd)/(qt*rv)))**((qt)*rv/cpd)
  else
  thl=tt*&
  (tt/tref)**(qt*(cpv-cpd)/cpd+ql*(cpl-cpv)/cpd+qi*(cpi-cpv)/cpd)*&
  (pp/(pref*(1+((qt-ql-qi)*rv)/((1-qt)*rd))))**(-(1-qt)*rd/cpd)*&
  exp(-ql*(srv-srl)/cpd-qi*(srv-sri)/cpd)*(1./(1.+(qt*rv)/((1.-qt)*rd)))**((1.-qt)*rd/cpd)*&
  (1./(1.+((1.-qt)*rd)/(qt*rv)))**((qt)*rv/cpd)
  endif
END FUNCTION

DOUBLE PRECISION FUNCTION mse(tt,zz,qt,ql,qi)
  double precision, intent(in) :: tt,zz,qt,ql,qi
  mse=(1.-qt)*(hrd+cpd*(tt-tref))+(qt-ql-qi)*(hrv+cpv*(tt-tref))+ql*(hrl+cpl*(tt-tref))+qi*(hri+cpi*(tt-tref))+grav*zz
END FUNCTION

DOUBLE PRECISION FUNCTION mseapr(tt,zz,qt,ql,qi)
  double precision, intent(in) :: tt,zz,qt,ql,qi
  mseapr=(1.-qt)*(hrd+cpd*(tt-tref))+(qt-ql-qi)*(hrv+cpd*(tt-tref))+ql*(hrl+cpd*(tt-tref))+qi*(hri+cpd*(tt-tref))+grav*zz
END FUNCTION

DOUBLE PRECISION FUNCTION thlapr(tt,pp,qt,ql,qi)
  double precision, intent(in) :: tt,pp,qt,ql,qi
  thlapr=tt*(pp/pref)**(-rd/cpd)*exp(-(ql*(hrv-hrl)+qi*(hrv-hri))/(cpd*tt))
END FUNCTION

DOUBLE PRECISION FUNCTION thllin(tt,pp,qt,ql,qi)
  double precision, intent(in) :: tt,pp,qt,ql,qi
  thllin=tt*(pp/pref)**(-rd/cpd)-(ql*(hrv-hrl)+qi*(hrv-hri))/(cpd*(pp/pref)**(rd/cpd))
END FUNCTION

DOUBLE PRECISION FUNCTION theapr(tt,pp,qt,ql,qi)
  double precision, intent(in) :: tt,pp,qt,ql,qi
  theapr=tt*(pp/pref)**(-rd/(cpd*(1+qt*cpv/cpd)))*exp(((qt-ql-qi)*(hrv-hrl)-qi*(hrl-hri))/(cpd*tt*(1+qt*cpv/cpd)))
END FUNCTION

SUBROUTINE qsat(isatice,tt,pp,qsatur)
  integer, intent(in) :: isatice
  double precision,intent(in) :: tt,pp
  double precision,intent(out) :: qsatur
  double precision :: qvsl,qvsi,ilratio,esattemp
  qvsl=rd/rv*esatl(tt)/(pp-(1.0-rd/rv)*esatl(tt))
  ilratio=max(0.0,min(1.0,(tt-tdn)/(tup-tdn)))
  if(isatice==1) then
    qvsi=rd/rv*esati(tt)/(pp-(1.0-rd/rv)*esati(tt))
  else
    qvsi=rd/rv*esatl(tt)/(pp-(1.0-rd/rv)*esatl(tt))
  endif
  qsatur = ilratio*qvsl+(1.0-ilratio)*qvsi
END SUBROUTINE

END MODULE
