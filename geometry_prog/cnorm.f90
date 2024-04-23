function cnorm(vect)
!Return the norm of a vector
real*8 :: vect(3),  cnorm
cnorm = dsqrt( vect(1)**2 + vect(2)**2 + vect(3)**2 ) 
end function cnorm
