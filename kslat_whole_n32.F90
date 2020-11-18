! Program to simulate the anisotropic KS model
! 2/14/20 Major bug fix with ghost_cells()
! New version (1/10/19) includes xoroshiro128+ PRNG, batch mode operation, and is ~10 times faster.
! This newer version (6/14/18) includes a linear term
#define KP1 2:nlat+1
#define KM1 0:nlat-1
#define K 1:nlat
program kslat
    use, intrinsic :: iso_fortran_env
    use xoroshiro128
    implicit none

    integer, parameter :: nlat=1024
    real(real32) :: c0x,c1x,c1y,c2x,c2y,c3,c4,c3x,delt
    real(real32), dimension(nlat,nlat) :: dzdx,dzdy,del2x,del2y,del4,zin,j,rnd,delx2,dely2
    real(real32), dimension(nlat,nlat) :: d2zdx2,d2zdy2,d2zdxdy,dzdxd,dzdyd
    real(real32), dimension(0:nlat+1,0:nlat+1) :: z,lap !zero index for ghost entries
    integer(int32) :: niter,nsteps,kk,iloop,a
    integer(int64) :: seed
    character(4) :: nch
    character(21) :: arg,ans1,ans2
    type(RNG_t) :: rng !for use with xoroshiro128+ module

! Simulation parameters
	call read_params(c0x,c1x,c1y,c4,c2x,c2y,c3x,c3,delt,niter,nsteps,seed,ans1,ans2)
	c3 = c3*sqrt(delt)
! set_seed() here is called with a 2 element array with seed repeated because the seed needs to be a 128-bit integer
    call rng%set_seed([seed, seed])

    if(ans2 == 'no') then; z = 0.0
    else
        open(unit=1,file=ans2,status='old',form='unformatted')
        read(1) zin; close(1)
        z(K,K)=zin
    end if
	
	j = spread(c0x*(1.-0.8*(/(a,a=1,nlat)/)/nlat),2,nlat) !needed to replace 'j' from do loop for whole array operations
    do iloop=1,nsteps
        write(6,*) iloop
        do kk=1,niter
			do a=1,nlat
				rnd(:,a)=rng%unif_01_n32M(nlat)
			end do
			!rnd = reshape(rng%unif_01_n32M(nlat*nlat),(/nlat,nlat/)) ! This is actually slower than the loop
            dzdx    = (z(KP1, K)- z(KM1, K))/2.
            dzdy    = (z(K, KP1)- z(K, KM1))/2.
            
            !old
            !d2zdx2  =  z(KP1, K)+ z(KM1, K)- 2.*z(K,K)
            !d2zdy2  =  z(K, KP1)+ z(K, KM1)- 2.*z(K,K)
            !del2x   = d2zdx2/(1.+dzdx**2)**1.5
            !del2y   = d2zdy2/(1.+dzdy**2)**1.5
            
            del2x=z(KP1,K)-2.*z(K,K)+z(KM1,K)
            del2y=z(K,KP1)-2.*z(K,K)+z(K,KM1)
            
            !old
            !d2zdxdy = (z(KP1,KP1)-z(KP1,KM1)-z(KM1,KP1)+z(KM1,KM1))/4.
            !lap(K,K)= ((1.+dzdx**2)*d2zdx2+(1.+dzdy**2)*d2zdy2-2.*dzdx*dzdy*d2zdxdy)/(1.+dzdx**2+dzdy**2)**1.5
            
            lap(K,K)=z(KM1,K)+z(KP1,K)+z(K,KM1)+z(K,KP1)-4.*z(K,K)
            delx2=((z(KP1,K)-z(K,K))**2+(z(KP1,K)-z(K,K))*(z(K,K)-z(KM1,K))+(z(K,K)-z(KM1,K))**2)/3.
			dely2=((z(K,KP1)-z(K,K))**2+(z(K,KP1)-z(K,K))*(z(K,K)-z(K,KM1))+(z(K,K)-z(K,KM1))**2)/3.
            
			call ghost_cells(lap)
			! no more of following -- old-- determine del4(i,j); del4 w/o small slope approximation
			!dzdxd  = (lap(KP1,K)-lap(KM1,K))/2.
			!dzdyd  = (lap(K,KP1)-lap(K,KM1))/2.
			!d2zdx2 =  lap(KP1,K)+lap(KM1,K)-2.*lap(K,K)
			!d2zdy2 =  lap(K,KP1)+lap(K,KM1)-2.*lap(K,K)
			!d2zdxdy= (lap(KP1,KP1)-lap(KP1,KM1)-lap(KM1,KP1)+lap(KM1,KM1))/4.
			!del4=	 ((1.+dzdxd**2)*d2zdx2+(1.+dzdyd**2)*d2zdy2-2.*dzdxd*dzdyd*d2zdxdy)/(1.+dzdxd**2+dzdyd**2)**1.5
			
			del4=lap(KM1,K)+lap(KP1,K)+lap(K,KM1)+lap(K,KP1)-4.*lap(K,K)
			
			
			! Evolve heights forward.
			!z(K,K) = z(K,K) + (c1x*del2x+c1y*del2y-c4*del4+c0x*dzdx+c2x*dzdx**2 &
			!+c2y*dzdy**2+c3x*dzdx**3)*delt+c3*(rnd-0.5)
			
			z(K,K)=z(K,K)+(c1x*del2x+c1y*del2y-c4*del4)*delt+c0x*dzdx*delt &
			+(c2x*delx2+c2y*dely2)*delt+(c3x*delx2*dzdx)*delt+c3*(rnd-0.5)
			
            call ghost_cells(z)
        end do
		! write lattice file
		if (ans1.ne.'none') then
			write(nch,"(I4.4)") iloop
			open(unit=1,file=nch//ans1,status='new',form='unformatted')
			write(1)z; close(1)
		end if
    end do
    
contains

subroutine read_params(c0x,c1x,c1y,c4,c2x,c2y,c3x,c3,delt,niter,nsteps,seed,ans1,ans2)
	real(real32) :: c0x,c1x,c1y,c2x,c2y,c3,c4,c3x,delt
	integer(int32) :: niter,nsteps
	integer(int64) :: seed
	character(21) :: arg,ans1,ans2
	
    call get_command_argument(1, arg); read(arg,'(F11.8)') c0x
    call get_command_argument(2, arg); read(arg,'(F11.8)') c1x
    call get_command_argument(3, arg); read(arg,'(F11.8)') c1y
    call get_command_argument(4, arg); read(arg,'(F11.8)') c4
    call get_command_argument(5, arg); read(arg,'(F11.8)') c2x
    call get_command_argument(6, arg); read(arg,'(F11.8)') c2y
    call get_command_argument(7, arg); read(arg,'(F11.8)') c3x
    call get_command_argument(8, arg); read(arg,'(F11.8)') c3
    call get_command_argument(9, arg); read(arg,'(F11.8)') delt
    call get_command_argument(10, arg); read(arg,'(I10)') niter
    call get_command_argument(11, arg); read(arg,'(I10)') nsteps
    call get_command_argument(12, arg); read(arg,'(I19)') seed
    call get_command_argument(13, ans1); !root lattice filename
    call get_command_argument(14, ans2); !starting lattice file
end subroutine read_params
! To handle periodic BCs with whole array operations, create fictious
! ghost cells on the boundary of the array that get read in instead of wrapping around
subroutine ghost_cells(A)
	real(real32), intent(inout), dimension(:,:) :: A
	integer(int32) :: n, s(2)
	s = shape(A)
	n = s(1)
	
	A(1,:) = A(n-1,:)
	A(:,1) = A(:,n-1)
	A(n,:) = A(2,:)
	A(:,n) = A(:,2)
end subroutine ghost_cells

end program kslat
