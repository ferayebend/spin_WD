module Kind_Types
implicit none
save
  integer, parameter, public::    &
      single = selected_real_kind(p=6,r=37),   &
      double = selected_real_kind(p=13,r=200)
	!  p = precision
	!  r = range
	! return smallest kind of real value with a minimum
	!  of p decimal digits and maximum range >10^r
end module kind_types            
!*********************************************
module Natural_Constants
use Kind_Types
implicit none
save
real(kind=double), parameter, public::     &
PI = 3.141592654d0,    &
G = 6.67259d-8,        &   ! gravitational constant	(dyne.cm^2/s)
M_sun = 1.989d33,      &   ! mass of the Sun   ( g )
R_sun = 6.96d10,       &   ! radius of the Sun (cm)
L_sun = 3.839e33,      &   ! luminosity of the Sun (erg/s)
Temp_sun = 5777,       &   ! effective temperature of the Sun (K)
c = 2.99792458d10,     &   ! speed of light	  (cm/s)
k_B = 1.3806513d-16,     &   ! Boltzmann's constant (erg/K)
planck = 6.6252d-27,    &   ! Plancks constant
m_p = 1.6725231d-24,   &   ! Mass of proton	  (g)
SB = 5.670399d-5,      &   ! Stefan-Boltzmann constant	(erg/cm^2.K^4.s)
minute = 60.d0,        &
hour = 3.6d3,          &   ! 1 hour is 3600 seconds
day = 8.64d4,          &   ! 1 day is 86400 seconds
week = 7.d0*day,       &   !
month = 30.d0*day,     &   ! convert time units to seconds
year = 365.25d0*day,   &   !
ly = c*year,           &
eV = 1.602177d-12,     &
pc = 3.26d0*ly,        &
kpc = pc*1.d3,         &
km = 1.d5                  ! km in cm's

end module Natural_Constants
!***************************************************
module Main_Parameters
use Natural_Constants
implicit none
save

real(kind=double), parameter, public::        &
!
mmw = 2.d0,          &
!
M_ch = 1.435*M_sun*(2.d0/mmw)**2,    & ! Chandrasekhar mass
                     !
alpha = 1.25d0,     &  ! power-law index of the accretion rate
					 !
tMAX = 1.d9*year,	 &    ! where to stop computing in time
				 	!
distance = 30d0*pc,   &
				     !
ksi = 1.0d0,        &  ! R_m = ksi * R_A
                     !
w_eq = 0.9d0,       & ! critical fastness parameter for torque equilibrium
                        !
inclination = PI/6.d0,  & ! inclination angle between the magnetic and rotation axis
			 !
A = 17,			& ! atomic number, for a half-half CO WD
			!
v_min = 1.0d10,      &  ! minimum frequency (for spectrum calculation) 
!
lambda_max = c / v_min,    &  ! maximum wavelength
!
kT_min = planck * v_min

integer, parameter, public::      &
order = 10,                        & ! then v_max = v_min*10^10
iMAX = 300   ! at how many logarithmically equal points the spectrum is calculated.

integer, parameter, public:: n_snap = 6 
! need to be even all the time
! if tMAX = 10^8 years then t_snap stars from 10^(8-n_snap/2)
! 
! eskiden n_snap=nint(log10(tMAX/year))
! if tMAX = 10^4 years then n_snap=4 
! and there will be 4 snapshots taken at t=10, t=100, t=10^3, t=10^4 years

real(kind=double), public, dimension(1:n_snap):: t_snap

end module Main_Parameters
!*****************************************************
module subroutines_and_functions
use Natural_Constants
use Main_Parameters
implicit none
public:: RK4adaptive, RK4, Wdot, data_write, initialize

real(kind=double), public:: t, Omega, dt,      &
 Mdot, W_kep, f, R_A, mu, M_star,              &
 torque_disk, torque_dip, torque,              &
 R_in, R_LC, R_co, w_s, period, jdot,          &
 R_star, L_disk, L_acc, B_star, L_0,	       &
 L_star, evaporation, Temp_0, Temp_star,       &
 Ms_in, Mdisk_in, Js_in, Jdisk_in, B_in,       &
 v, T_p, T_s, T_disk_ch, W_kep_star, x_d

 
real(kind=double), public:: Mdot_0, t_0, M_0, j_0,r_0, J_star, I_star, flux
 
integer, public:: ss, snap, shut_off

contains
!***************************************************
!*******************************************
!************************************
subroutine simpson(a,b,func,integral)
!use Kind_Types
implicit none
integer, parameter:: n=400
integer:: i
real(kind=double), intent(in):: a,b
real(kind=double), intent(out)::integral
real(kind=double):: dx
real(kind=double), dimension(0:n):: c, x, y
real(kind=double), external :: func

! calculate the step, coefficients of formula

    dx = (b-a)/dble(n)
    c(0) = dx/3.d0
    do i = 1, n/2
        c(2*i-1) = (4.d0/3.d0)*dx
        c(2*i) = (2.d0/3.d0)*dx
    end do
    c(n) = c(0)
    integral = 0

! calculate the interpolate points and value on them

    do i = 0, n
        x(i) = a + i*dx
        y(i) = func(x(i))
    end do

! calculate the approximate integral

    do i = 0, n
        integral = integral + c(i)*y(i)
    end do

end subroutine simpson

!*********************************************
real(kind=double) function planck_dist(x)
!use Natural_Constants
!use Main_Parameters
implicit none
real(kind=double), intent(in) :: x
real(kind=double):: Temp, hv_kT

Temp = T_s*((x**(-3))*(1.d0-jdot/x**0.5))**0.25

!Temp  = T_0*x**(-0.75)

hv_kT = planck*v/(k_B*Temp)

    if ( hv_kT <1.d-3) then

	    planck_dist =  x / hv_kT

    else if ( hv_kT > 12.d0 ) then

       planck_dist = x * dexp(-hv_kT)

	else

	     planck_dist = x / ( dexp(hv_kT) - 1.d0)

	end if

end function planck_dist
!********************************
subroutine spectrum
implicit none
real(kind=double):: res, constant, a, b, base, Fv, R_out, pov_inc
integer:: i

!T_s = (3.d0*G*M_star*Mdot/(8.d0*pi*R_in**3*SB))**0.25

R_out = r_0 * (1.d0 + t/t_0)**0.5 ! for bound-free opacity

pov_inc = 0.0*pi  ! looking at the disk face-on

constant = 4.d0*pi*planck*cos(pov_inc)*(R_in/distance)**2/c**2 ! uses disk inclination wrt pov

base = 10.d0**(dble(order)/dble(iMAX))

write(13,'(A9,F6.2,1x,A6,F7.2,1x,A13)') "# t = 10^", log10(t/year), "year  ", T_s, "K temperature"

a = 1.d0
b = R_out/R_in

v = v_min

do i=1, imax

     
     call simpson(a,b,planck_dist,res)
	 

     Fv = res*constant*v**3

     write (unit=13,fmt=109) c/v, planck*v/eV, Fv, v*Fv  
	 ! note v*Fv = lambda F_lambda

     v = v_min*base**i

	 if (dlog10(v*Fv) < -40) then
	  exit
     end if

end do

  write(13,*) ""
109 FORMAT (1x, ES14.5, 2x, ES14.5, 2x, ES14.5, 2x, ES14.5, 2x, F7.2)
end subroutine spectrum
!***************************************************
!***************************************************
SUBROUTINE initialize
real(kind=double):: Sigma_0, C_0, nu_0, mmw_disk, alpha_SS, kappa_0
mmw_disk = 2.0 ! for metallic disk. used to be 0.625
alpha_SS = 0.1 ! Shakura-Sunyaev alpha parameter
kappa_0 = 4.d25 ! *Z*(1+X) for Z metal and X hydrogen fraction http://arxiv.org/abs/astro-ph/0205212v1

OPEN (unit=21, file="data.in",status="old")
READ (21, *) B_in, Ms_in, Js_in, Mdisk_in, Jdisk_in
WRITE (*, 103)  B_in, Ms_in, Js_in, Mdisk_in, Jdisk_in
103 FORMAT (5(2x,ES12.4))

t = 0.d0

dt = 1.d-5 ! seconds

! input from file

B_star = B_in  ! Gauss, initial magnetic field

M_0 = Mdisk_in*M_sun  ! initial mass of the disk

j_0 = Jdisk_in/M_0 ! initial average specific angular momentum

J_star = Js_in  ! g cm^2/s, initial angular momentum of the star

M_star = Ms_in * M_sun	      ! initial mass of the WD in "g"

Temp_0 = 1e8 ! K, initial temp for an isothermal star

L_0 = 5.75d5*3.45*(M_star/M_sun)*Temp_0**(7./2) ! initial luminosity in erg/s

!

r_0 = (3.3558/1.7030)**2*j_0**2/(G*M_star)  !  replace

Sigma_0 = M_0 / (4.d0*pi*r_0**2*3.3558d-5)

C_0 = alpha_SS**(8./7.)*(27.*kappa_0/SB)**(1./7.)    &
      *(k_B/(mmw_disk*m_p))**(15./14.)*(G*M_star)**(-5./14.)

nu_0 = C_0*r_0**(15./14.)*Sigma_0**(3./7.)

!t_0 = 2.150*year*(j_0/1d20)**(7.0/3.0)*(M_0/(1d-4*M_sun))**(-2.0/3.0) ! from the disk model
t_0 = r_0**2/(0.75d0*nu_0)

print*, t_0
!print*, Sigma_0

Mdot_0 = (alpha -1.0)*M_0/t_0


!--------------------------------
! parameters below are determined from the given above

R_star = 0.0126*R_sun*(2.d0/mmw)*(M_star/M_sun)**(-1.0/3.0)*(1.d0-(M_star/M_ch)**(4.d0/3.d0))**0.5d0

I_star = 3.21601*1d50*(M_star/M_sun)**0.34158*(1.d0-(M_star/(1.437*M_sun))**1.2499)**1.43773 

flux = B_star*R_star**2

mu = 0.5*B_star*R_star**3

Omega = J_star/I_star  ! angular velocity

Mdot = Mdot_0

R_LC = c / Omega

R_co = (G*M_star/Omega**2)**(1.d0/3.d0)

R_A = ksi*( mu**2 /(SQRT(2.d0*G*M_star)*Mdot) )**(2.d0/7.d0)

R_in = max(R_star,R_A) ! if R_star>R_A then the inner radius of the disk is R_star

W_kep = SQRT(G*M_star/R_in**3) ! kepler angular velocity at the inner radius

W_kep_star = SQRT(G*M_star/R_star**3) ! kepler angular velocity at the surface of the star

w_s = Omega / W_kep       !fastness parameter

f = R_in / R_LC

jdot = 1.d0 - w_s/w_eq

torque_disk = Mdot*SQRT(G*M_star*R_in)*jdot

torque_dip = - 2.d0 * mu**2 * sin(inclination)**2 * Omega**3  /(3.d0*c**3) 

! radiative parameters

L_star = L_0

Temp_star = ((L_star*R_sun**2)/(L_sun*R_star**2))**(1.d0/4)*Temp_sun ! effective temperature

evaporation = (L_star/L_sun)**(12.d0/7)/(-5./7*9.41d6*year*12./A*(mmw/2.d0)**(4./3)*(M_star/M_sun)**(5./7)) &
		- 5./7*(L_star*Mdot)/(L_sun*M_star) ! Mestel dL/dt

T_s = (3.d0*G*M_star*Mdot/(8.d0*pi*R_in**3*SB))**0.25 ! temp at inner disk radius

T_p = 1000  ! K MRI shut off temperature

T_disk_ch = T_s

shut_off = 0 ! MRI shutoff 

if ( f > 1.d0) then
   torque = torque_dip 
   ss = 0
else 
   torque = torque_disk  
   ss = 1
end if

RETURN
END SUBROUTINE initialize

SUBROUTINE MRI_shutoff
real(kind=double):: Sigma_0, C_0, nu_0, mmw_disk, alpha_SS, kappa_0, Md_n, t_n
mmw_disk = 2.0 ! for metallic disk. used to be 0.625
alpha_SS = 0.01 ! Shakura-Sunyaev alpha parameter
kappa_0 = 4.d25 ! *Z*(1+X) for Z metal and X hydrogen fraction http://arxiv.org/abs/astro-ph/0205212v1

! changes t_0

r_0 = (3.3558/1.7030)**2*j_0**2/(G*M_star)  !  replace

Sigma_0 = M_0 / (4.d0*pi*r_0**2*3.3558d-5)

C_0 = alpha_SS**(8./7.)*(27.*kappa_0/SB)**(1./7.)    &
      *(k_B/(mmw_disk*m_p))**(15./14.)*(G*M_star)**(-5./14.)

nu_0 = C_0*r_0**(15./14.)*Sigma_0**(3./7.)

!t_0 = 2.150*year*(j_0/1d20)**(7.0/3.0)*(M_0/(1d-4*M_sun))**(-2.0/3.0) ! from the disk model
t_0 = r_0**2/(0.75d0*nu_0)

Mdot_0 = (alpha -1.0)*M_0/t_0
!Md_n = Mdot_0 * (1.0 + t/t_0)**(-alpha)/(1.0 + t/t_n)**(-alpha) ! for continuity

!Mdot_0 = Md_n

!t_0 = t_n

RETURN
END SUBROUTINE MRI_shutoff


!***************************************************
! This subroutine takes a 4th order Runge-Kutta step
! Starting with angular velocity W at time t
! take a step with dt
! and find new t and W
! and return them to RK4adaptive
SUBROUTINE RK4(W,t,dt,W1,t1)
implicit none
real(kind=double),intent(in) :: W,t,dt
real(kind=double), intent(out)::W1,t1
real(kind=double):: dW1, dW2, dW3, dW4, dW

dW1 = dt * Wdot(t,W)

dW2 = dt * Wdot(t + 0.5d0 * dt, W + 0.5d0 * dW1)

dW3 = dt * Wdot(t + 0.5d0 * dt, W + 0.5d0 * dW2)

dW4 = dt * Wdot(t + dt, W + dW3)

! Take the weighted average of the above dW by giving 2 weights
! to the midpoints

dW = (dW1 + 2.d0*dW2 + 2.d0*dW3 + dW4) / 6.0d0

W1 = W + dW
t1 = t + dt

RETURN
END SUBROUTINE RK4

!**************************************************
! This subroutine adjusts the time steps used in RK4 
! by looking  at the error made in taking a double RK4 step 
!  instead of two successive RK4 steps
! Given W, t, dt find new W, t and dt

SUBROUTINE RK4adaptive(W,t,dt,W_new,t_new, dt_new)
implicit none
real(kind=double),intent(in) :: W,t,dt
real(kind=double),intent(out):: W_new, t_new, dt_new
real(kind=double):: error, WTEMP1,WTEMP2, WTEMP3, tTEMP1, tTEMP2, tTEMP3
integer:: i

! Parameters relevant to the numerical procedure
real(kind=double), parameter:: delta=1.0d-12 ! desired accuracy
real(kind=double), parameter:: safety=0.7  ! safety factor

! Take two successive steps
 call RK4(W,t,dt,WTEMP1,tTEMP1)    
 call RK4(WTEMP1,tTEMP1,dt,WTEMP2, tTEMP2)

! Take a single step with 2dt    
 call RK4(W,t,2.d0*dt, WTEMP3, tTEMP3)

 ! Find the relative error
 !error = ABS((WTEMP3-WTEMP2)/WTEMP2)
 error = max(ABS((WTEMP3-WTEMP2)/WTEMP2),1.d-4*delta)

! if we are making a small error  (error < delta)
! then make the time steps larger

if (error <= delta) then	
   	 
      dt_new = safety*dt*ABS(delta/(error))**0.2
! ... but not more than 5 times
      dt_new = min(dt_new, 5.d0*dt)	

! ...but if we are making a big error  (error > delta)
! then make the time step smaller

else
						
   dt_new = safety*dt*ABS(delta/(error))**0.25
    ! ... but not more than 10 times
   dt_new = max(dt_new,0.1d0*dt)	  
		   
end if

! do not let the time-step be larger than 1/500th of tMAX
  
 dt_new = min(tMAX/500.d0,dt_new)


! do not let "t" overshoot tMAX too much
! if "t" will be greater than the tMAX in the next step
! make it just slightly greater than required to reach tMAX
   if (dt_new > tMAX - t) then         
       dt_new = 1.00001d0*(tMAX - t)  	     
   end if

! do not let "t" overshoot t_snap too much
! if "t" will be greater than the t_snap in the next step
! make it just slightly greater than just required to reach t_snap

  do i=1,n_snap
     if (t<t_snap(i) .and. dt_new > t_snap(i) - t) then         
         dt_new = 1.00000001d0*(t_snap(i) - t)
	   snap=1 
           exit
     end if
   end do 


! the final RK4 step with the new dt
call RK4(W,t,dt_new,W_new,t_new)

RETURN
END SUBROUTINE RK4adaptive
!**************************************************
REAL(kind=double) FUNCTION Wdot(t,Omega)
implicit none
real(kind=double), intent(in):: t, Omega


Mdot = Mdot_0 * (1.0 + t/t_0)**(-alpha)

R_LC = c / Omega

R_co = (G*M_star/Omega**2)**(1.d0/3.d0)

R_A = ksi*( mu**2 /(SQRT(2.d0*G*M_star)*Mdot) )**(2.d0/7.d0)

R_in = max(R_star,R_A) ! if R_star>R_A then the inner radius of the disk is R_star

W_kep = SQRT(G*M_star/R_in**3) ! kepler angular velocity at the inner radius

W_kep_star = SQRT(G*M_star/R_star**3) ! at stellar surface

w_s = Omega / W_kep       !fastness parameter

f = R_in / R_LC

jdot = 1.d0 - w_s/w_eq

torque_disk = Mdot*SQRT(G*M_star*R_in)*jdot

torque_dip = - 2.d0 * mu**2 * sin(inclination)**2 * Omega**3  /(3.d0*c**3) 

evaporation = (L_star/L_sun)**(12.d0/7)/(-5./7*9.41d6*year*12./A*(mmw/2.d0)**(4./3)*(M_star/M_sun)**(5./7)) &
		- 5./7*(L_star*Mdot)/(L_sun*M_star) ! Mestel dL/dt

x_d = 10 ! R_out/R_in ratio where MRI shuts off at T_p 

T_s = (3.d0*G*M_star*Mdot/(8.d0*pi*R_in**3*SB))**0.25

T_disk_ch = T_s*((x_d**(-3))*(1.d0-jdot/x_d**0.5))**0.25 ! disk inner radius temperature

if ( f > 1.d0 .or. T_disk_ch < T_p) then ! if the inner disk temperature smaller than 1000 K
   torque = torque_dip 
   ss = 0
else 
   torque = torque_disk  
   ss = 1
end if


Wdot = torque / I_star

RETURN
END FUNCTION Wdot

!***************************************************
SUBROUTINE data_write
IMPLICIT NONE

period = 2.d0*PI / Omega

L_acc = G * M_star * Mdot / R_star 
L_disk =  0.5d0*G * M_star * Mdot / R_in
Temp_star = ((L_star*R_sun**2)/(L_sun*R_star**2))**(1.d0/4)*Temp_sun

write(11,100) t/year, M_star/M_sun, J_star/1.d50, Omega, Period, I_star/1.d50,  &
              w_s, torque/1.d40, ss, Temp_star, B_star*1.0d-6
write(12,101) t/year, Mdot*year/M_sun, R_in/R_star, L_acc, L_disk, f, ss,       &
	      t_0/year, Mdot_0*year/M_sun

100 FORMAT (1x,ES12.4,2x,F8.6,2x,ES12.4,2x,F9.3,2x,ES12.4,2x,F10.6, &
            2x,F14.4,2x,ES12.4,2x, I1, 2x, F9.2,2x, F7.2)
101 FORMAT (6(2x,ES12.4),2x, I1,2(2x,ES12.4))
RETURN
END SUBROUTINE data_write

!***************************************************

SUBROUTINE update
IMPLICIT NONE

  M_star = M_star + Mdot*dt

  J_star = J_star + torque*dt

  R_star = 0.0126*R_sun*(2.d0/mmw)*(M_star/M_sun)**(-1.0/3.0)*(1.d0-(M_star/M_ch)**(4.d0/3.d0))**0.5d0   ! Neuenberg ile belirle

  I_star = 3.21601*1d50*(M_star/M_sun)**0.34158*(1.d0-(M_star/(1.437*M_sun))**1.2499)**1.43773 
 
  L_star = L_star + evaporation*L_sun*dt

  B_star = flux/R_star**2

  mu = 0.5*B_star*R_star**3

  Omega = J_star/I_star  ! angular velocity



END SUBROUTINE update

!***************************************************

end module subroutines_and_functions
!*****************************************************
!*****************************************************
!*****************************************************
!                                                   !*
           PROGRAM spin                             !*
!                                                   !*
!*****************************************************
use Natural_Constants
use Main_Parameters
use subroutines_and_functions
implicit none

! these are to write the data in logarithmic time steps
real(kind=single):: base=1.005   
integer:: log_t, log_tOLD, log_tREF


integer:: i

log_tREF = nint(log10(tMAX/year))-n_snap/2

OPEN (unit=11, file="star.out",status="replace")
OPEN (unit=12, file="disk.out",status="replace")
OPEN (unit=13, file="spectrum.out",status="replace")

write(*,'(1x,A10,I1,1x,A5)') "tMAX = 10^", nint(log10(tMAX/year)), "years"
write(*,*) log_tREf
write(*,*) ""
write(*,'(1x,A30)') "snapshots are to be taken at"
    do i= 1,n_snap/2  ! two calculations in each step
        t_snap(2*i-1) = 10**(log_tREF+i)*year
	t_snap(2*i) = 3*10**(log_tREF+i)*year
        write(*,'(1x,A3,F6.2,1x,A5)') "10^", log10(t_snap(i)/year),"years"
    end do
        write(*,*) ""

write(11,'(A109)') '#t/year,    M_star,     J_star/1E50,  Omega,      Period, I_star/1E50,        w_s,       torque/1E40, ss'
write(12,'(A87)') "#t/year, Mdot/Msun/yr, R_in/R_star, L_acc, L_disk, f, ss, t0, Mdot0"


! initialize everything
i = 0

call initialize

call data_write

log_tOLD = floor(log(dt)/log(base))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
time: do

  if (t >= tMAX) then  	   
       exit time  		 
  end if 

  if (M_star > 0.97*M_ch) then
    print*, "M_s reached M_ch"
	exit time
  end if

  if (Omega > W_kep_star) then
    print*, "Star is rotating faster than breakup speed"
	exit time
  end if

  if (T_disk_ch < T_p) then 
    if (shut_off .eq. 0) then
       open(unit=23, file="MRI_turnoff.out")!,status="replace")
       123 FORMAT (F5.1,F9.3,3(2x,ES12.4))
       write(23,123) B_star*1e-6, w_s, period, Mdot*year/M_sun, t/year
       Mdot = 0 ! not working
       !print*, "alpha_ss = 0.01 MRI shutoff", t/year
       shut_off = 1
    end if
       !close(unit=33)
  end if


  call RK4adaptive(Omega,t,dt, Omega,t,dt)	! takes the values one step ahead 
    
     if (snap==1) then
       call spectrum
     end if
      snap=0

    log_t = floor(log(t)/log(base)) 
	            
        if (log_t /= log_tOLD) then
		  
           call data_write
		        
        end if
        	
    call update
 
    log_tOLD = log_t
	i = i + 1		 ! numerates the time steps

end do time      !  end the iteration

STOP
END PROGRAM spin
!**************************************************

