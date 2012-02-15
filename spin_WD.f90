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
k = 1.3806513d-16,     &   ! Boltzmann's constant (erg/K)
m_p = 1.6725231d-24,   &   ! Mass of proton	  (g)
SB = 5.670399d-5,      &   ! Stefan-Boltzmann constant	(erg/cm^2.K^4.s)
minute = 60.d0,        &
hour = 3.6d3,          &   ! 1 hour is 3600 seconds
day = 8.64d4,          &   ! 1 day is 86400 seconds
week = 7.d0*day,       &   !
month = 30.d0*day,     &   ! convert time units to seconds
year = 365.25d0*day,   &   !
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
tMAX = 3.d8*year,	 &    ! where to stop computing in time
				     !
ksi = 1.0d0,        &  ! R_m = ksi * R_A
                     !
w_eq = 0.9d0,       & ! critical fastness parameter for torque equilibrium
                        !
inclination = PI/3.d0,  & ! inclination angle between the magnetic and rotation axis
			 !
A = 17			 ! atomic number, for a half-half CO WD


	 
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
 Ms_in, Mdisk_in, Js_in, Jdisk_in, B_in

 
real(kind=double), public:: Mdot_0, t_0, M_0, j_0, J_star, I_star, flux
 
integer, public:: ss

contains
!***************************************************
SUBROUTINE initialize

OPEN (unit=21, file="data.in",status="old")
READ (21, *) B_in, Ms_in, Js_in, Mdisk_in, Jdisk_in
WRITE (*, 103)  B_in, Ms_in, Js_in, Mdisk_in, Jdisk_in
103 FORMAT (5(2x,ES12.4))

t = 0.d0

dt = 1.d6 ! seconds

B_star = B_in  ! Gauss, initial magnetic field

M_0 = Mdisk_in*M_sun  ! initial mass of the disk

j_0 = Jdisk_in/M_0 ! initial average specific angular momentum

t_0 = 2.150*year*(j_0/1d20)**(7.0/3.0)*(M_0/(1d-4*M_sun))**(-2.0/3.0) ! from the disk model

Mdot_0 = (alpha -1.0)*M_0/t_0

J_star = Js_in  ! g cm^2/s, initial angular momentum of the star

M_star = Ms_in * M_sun	      ! initial mass of the WD in "g"

Temp_0 = 1e8 ! K, initial temp for an isothermal star

L_0 = 5.75d5*3.45*(M_star/M_sun)*Temp_0**(7./2) ! initial luminosity in erg/s

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



if ( f > 1.d0) then
   torque = torque_dip 
   ss = 0
else 
   torque = torque_disk  
   ss = 1
end if

RETURN
END SUBROUTINE initialize

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

! Parameters relevant to the numerical procedure
real(kind=double), parameter:: delta=1.0d-8 ! desired accuracy
real(kind=double), parameter:: safety=0.7  ! safety factor

! Take two successive steps
 call RK4(W,t,dt,WTEMP1,tTEMP1)    
 call RK4(WTEMP1,tTEMP1,dt,WTEMP2, tTEMP2)

! Take a single step with 2dt    
 call RK4(W,t,2.d0*dt, WTEMP3, tTEMP3)

 ! Find the relative error
error = ABS((WTEMP3-WTEMP2)/WTEMP2)


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

w_s = Omega / W_kep       !fastness parameter

f = R_in / R_LC

jdot = 1.d0 - w_s/w_eq

torque_disk = Mdot*SQRT(G*M_star*R_in)*jdot

torque_dip = - 2.d0 * mu**2 * sin(inclination)**2 * Omega**3  /(3.d0*c**3) 

evaporation = (L_star/L_sun)**(12.d0/7)/(-5./7*9.41d6*year*12./A*(mmw/2.d0)**(4./3)*(M_star/M_sun)**(5./7)) &
		- 5./7*(L_star*Mdot)/(L_sun*M_star) ! Mestel dL/dt


if ( f > 1.d0) then
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
write(12,101) t/year, Mdot, R_in/R_star, L_acc, L_disk, f, ss  


100 FORMAT (1x,ES12.4,2x,F8.6,2x,ES12.4,2x,F9.3,2x,ES12.4,2x,F10.6, &
            2x,F14.4,2x,F14.4,2x, I1, 2x, F9.2,2x, F8.2)
101 FORMAT (6(2x,ES12.4),2x, I1)
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
integer:: log_t, log_tOLD


integer:: i

OPEN (unit=11, file="star.out",status="replace")
OPEN (unit=12, file="disk.out",status="replace")

write(11,'(A106)') '#t/year,    M_star,     J_star/1E50,  Omega,      Period, I_star/1E50,        w_s,       torque/1E40, ss'
write(12,'(A80)') "#t/year, Mdot, R_in/R_star, L_acc, L_disk, f, ss"

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

  if (M_star > 0.95*M_ch) then
    print*, "M_s reached M_ch"
	exit time
  end if

  call RK4adaptive(Omega,t,dt, Omega,t,dt)	! takes the values one step ahead 


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

