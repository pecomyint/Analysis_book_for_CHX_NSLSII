!> Module for pseudo random number generation. The internal pseudo random
!> generator is the xoroshiro128plus method.
! Modifications made by JJB Jan/Feb 2020 for faster generation of
! matrices of 32 bit floats 
module xoroshiro128
  implicit none
  private

  ! A 64 bit floating point type
  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: sp = kind(0.0)
  ! A 32 bit integer type
  integer, parameter :: i4 = selected_int_kind(9)
  ! A 64 bit integer type
  integer, parameter :: i8 = selected_int_kind(18)
  !> Random number generator type, which contains the state
  type rng_t
     !> The rng state (always use your own seed)
     integer(i8), private       :: s(2) = [123456789_i8, 987654321_i8]
     integer(i8), private       :: separator(32) ! Separate cache lines (parallel use)
   contains
     procedure, non_overridable :: set_seed    ! Seed the generator
     procedure, non_overridable :: jump        ! Jump function (see below)
     procedure, non_overridable :: int_4       ! 4-byte random integer
     procedure, non_overridable :: int_8       ! 8-byte random integer
     procedure, non_overridable :: unif_01_n32M! Uniform (0,1] 32 bit real
     procedure, non_overridable :: next        ! Internal method
     procedure, non_overridable :: next_n      ! Internal method
  end type rng_t
  
    public :: rng_t

contains
  !> Set a seed for the rng
  subroutine set_seed(self, the_seed)
    class(rng_t), intent(inout) :: self
    integer(i8), intent(in)     :: the_seed(2)

    self%s = the_seed

    ! Simulate calls to next() to improve randomness of first number
    call self%jump()
  end subroutine set_seed

  ! This is the jump function for the generator. It is equivalent
  ! to 2^64 calls to next(); it can be used to generate 2^64
  ! non-overlapping subsequences for parallel computations.
  subroutine jump(self)
    class(rng_t), intent(inout) :: self
    integer                     :: i, b
    integer(i8)                 :: t(2), dummy

    ! The signed equivalent of the unsigned constants
    integer(i8), parameter      :: jmp_c(2) = &
         (/-4707382666127344949_i8, -2852180941702784734_i8/)

    t = 0
    do i = 1, 2
       do b = 0, 63
          if (iand(jmp_c(i), shiftl(1_i8, b)) /= 0) then
             t = ieor(t, self%s)
          end if
          dummy = self%next()
       end do
    end do

    self%s = t
  end subroutine jump

  !> Return 4-byte integer
  integer(i4) function int_4(self)
    class(rng_t), intent(inout) :: self
    int_4 = int(self%next(), i4)
  end function int_4

  !> Return 8-byte integer
  integer(i8) function int_8(self)
    class(rng_t), intent(inout) :: self
    int_8 = self%next()
  end function int_8
  
  function unif_01_n32M(self,n) result(res)
    class(rng_t), intent(inout) :: self
    integer(i8), dimension(n/2) :: x
    real(sp)			        :: res(n), moldir(2)
    integer(i4), intent(in)		:: n

	x   = self%next_n(n/2)
	call mvbits(127_i8, 0, 9, x, 23)
	call mvbits(127_i8, 0, 9, x, 55)
	res = transfer(x, moldir) - 1.0
  end function unif_01_n32M

  !> Interal routine: get the next value (returned as 64 bit signed integer)
  function next(self) result(res)
    class(rng_t), intent(inout) :: self
    integer(i8)                 :: res
    integer(i8)                 :: t(2)

    t         = self%s
    res       = t(1) + t(2)
    t(2)      = ieor(t(1), t(2))
    self%s(1) = ieor(ieor(rotl(t(1), 55), t(2)), shiftl(t(2), 14))
    self%s(2) = rotl(t(2), 36)
  end function next
  
  !> Internal routine: get the next n values (returned as 64 bit signed integer)
  ! TODO: Can we do more work in the loop despite sequential iteration?
  ! TODO: e.g. Borrow parallel rng gen functionality to have multiple
  ! TODO: iterations running alongside each other:
  ! TODO: next(2)=f(next(1)) | next(2^64 +2) = f(next(2^64+1)) | etc
  function next_n(self, n) result(res)
    class(rng_t), intent(inout) :: self
    integer(i8), dimension(n)   :: res, t(2)
    integer(i4)					:: i
    integer(i4), intent(in)		:: n
    
	t = self%s
	do i=1,n
		res(i)  = t(1) + t(2)
		t(2)    = ieor(t(1), t(2))
		t(1) 	= ieor(ieor(rotl(t(1), 55), t(2)), shiftl(t(2), 14))
		t(2) 	= rotl(t(2), 36)
	end do
	self%s = t
  end function next_n

  !> Helper function for next()
  pure function rotl(x, k) result(res)
    integer(i8), intent(in) :: x
    integer, intent(in)     :: k
    integer(i8)             :: res

    res = ior(shiftl(x, k), shiftr(x, 64 - k))
  end function rotl

end module xoroshiro128
