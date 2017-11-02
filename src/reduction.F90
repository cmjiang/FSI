subroutine reduction(matB, matH, work, N, L, b, q)
  implicit none
  integer, intent(in)  :: N, L, b, q
  double precision, intent(in) :: matB(N,N,L)
  double precision, intent(inout) :: matH(N,N,b), work(N,N,b)
  integer :: LN, i, j, k, c
  double precision, parameter :: one = 1.0d0, zero = 0.0d0
  integer :: omp_get_thread_num, thn
  LN = N*L
  c = L/b   
!$omp parallel private(k)
!$omp do
  do i = 1, b
    matH(:,:,i) = matB(:,:,i*c-q)
    do j = 1, c-1
      k = i*c-j-q
      if ( k==0 ) then
        k = L
      else
        k = mod(k, L)
      endif
      call dgemm('n', 'n', N, N, N, one, matH(:,:,i), N, &
                 matB(:,:,k), N, zero, work(:,:,i), N)
      matH(:,:,i) = work(:,:,i)
    enddo
  enddo
!$omp end do
!$omp end parallel
end subroutine reduction
