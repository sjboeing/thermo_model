module parceldata
implicit none

integer,parameter :: ztop=4000
double precision,parameter :: dz       = 5
double precision,parameter :: rd       = 287.04           !<    *gas constant for dry air.
double precision,parameter :: rv       = 461.5            !<    *gas constant for water vapor
double precision,parameter :: grav     = 9.81             !<    *little g (standard gravity)
double precision,parameter :: tref     = 273.15           !<    *reference temperature
double precision,parameter :: pref     = 100000.          !<    *reference pressure
double precision,dimension(0:ztop) :: z
double precision,dimension(0:ztop) :: piso
double precision,dimension(0:ztop) :: penv
double precision,dimension(0:ztop) :: tenv
double precision,dimension(0:ztop) :: tvenv
double precision,dimension(0:ztop) :: entenv
double precision,dimension(0:ztop) :: mseenv
double precision,dimension(0:ztop) :: mseaprenv
double precision,dimension(0:ztop) :: thlaprenv
double precision,dimension(0:ztop) :: thllinenv
double precision,dimension(0:ztop) :: theaprenv
double precision,dimension(0:ztop) :: thlenv
double precision,dimension(0:ztop) :: thsenv
double precision,dimension(0:ztop) :: qtenv
double precision,dimension(0:ztop) :: pp
double precision,dimension(0:ztop) :: tent
double precision,dimension(0:ztop) :: tths
double precision,dimension(0:ztop) :: tthl
double precision,dimension(0:ztop) :: tmse
double precision,dimension(0:ztop) :: tmseapr
double precision,dimension(0:ztop) :: tthlapr
double precision,dimension(0:ztop) :: tthllin
double precision,dimension(0:ztop) :: ttheapr
double precision,dimension(0:ztop) :: tvent
double precision,dimension(0:ztop) :: tvths
double precision,dimension(0:ztop) :: tvthl
double precision,dimension(0:ztop) :: tvmse
double precision,dimension(0:ztop) :: tvmseapr
double precision,dimension(0:ztop) :: tvthlapr
double precision,dimension(0:ztop) :: tvthllin
double precision,dimension(0:ztop) :: tvtheapr
double precision :: ecape,dtvmax
end module parceldata
