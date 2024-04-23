subroutine H_diag(size)

!Diagonalize the Hamiltonian.
!Return eigenvalues and diagonalization matrix.

use declare
integer :: size

allocate(eigval(size))
eigval = 0.0d0
allocate(work(1))

lwork=work(1)

call dsyev('V','L',size,H,size,eigval,work,-1,info)
lwork=work(1)

deallocate (work)
allocate(work(lwork))

call dsyev('V','L',size,H,size,eigval,work,lwork,info)

deallocate (work)

end
