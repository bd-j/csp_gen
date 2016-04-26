FUNCTION expi(arg)
  ! Computes the exponential integral Ei(x) for x > 0.
  ! Parameters:
  !   `eps` is the relative error, or absolute error near the zero of Ei at x=0.3725;
  !   `eul` is Euler’s constant
  !   `maxit` is the maximum number of iterations allowed;
  !   `fmin` is a number near the smallest representable floating-point number.
  ! Adapted from the NR *public domain* code `ei` (http://numerical.recipes/pubdom/nr.f90.txt)
  implicit none
  
  INTEGER, PARAMETER :: maxit=1000
  REAL expi, arg
  REAL, PARAMETER :: eps=6.e-8, eul=.57721566, fmin=1.e-30
  INTEGER k
  REAL fact,prev,sum,term
  
  if (arg.le.0.) then
     write(*,*) "expi: arg < 0"
     EXIT
  endif
  
  ! Special case: avoid failure of convergence test because of underflow.
  if (arg.lt.fmin) then   
     expi = log(arg) + eul
  else if (arg.le.-log(eps)) then ! Use power series.
     sum = 0.
     fact = 1.
     do k=1, maxit
        fact = fact * arg / k
        term = fact / k
        sum = sum + term
        if (term.lt.(eps*sum)) continue
     enddo
     expi = sum + log(arg) + eul
  else ! Use asymptotic series.
     sum = 0. ! Start with second term.
     term = 1.
     do k=1, maxit
        prev = term
        term = term * k / arg
        ! Since final sum is greater than one, term itself approximates the relative error.
        if (term.lt.eps) continue
        if (term.lt.prev) then
           sum = sum + term ! Still converging: add new term.
        else
           sum = sum - prev ! Diverging: subtract previous term and exit.
           continue
        endif
     enddo
     expi = exp(arg) * (1. + sum) / arg
  endif
  return
END FUNCTION expi
