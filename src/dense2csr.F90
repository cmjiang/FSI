subroutine dense2csr(vals,ci,rp,matsB,N,L,nnz)
  
  ! dense2crs converts the dense format of p-cyclic matrix into
  ! the CRS format. The dense format is given as an array of
  ! N-by-N dense matrices in the length of L. The non-zeros values,
  ! column indices and row pointers are generated as the output.
  
  implicit none
  double precision, intent(in) :: matsB(N,N,L)
  integer, intent(in) :: N, L, nnz
  double precision, intent(out) :: vals(nnz)
  integer, intent(out) :: ci(nnz), rp(N*L+1)
  integer :: ptrv, ptrc, i, k, ind(N)
  ptrv = 1
  ptrc = 1
  do i = 1, N
    ind(i) = i
  enddo
  ! generate the first block row
  do i = 1, N
    vals(ptrv) = 1
    ptrv = ptrv + 1
    ci(ptrc) = i
    ptrc = ptrc + 1
    vals(ptrv: ptrv+N-1) = matsB(i, 1:N, 1)
    ptrv = ptrv + N
    ci(ptrc: ptrc+N-1) = ind(1:N) + (L-1)*N
    ptrc = ptrc + N
  enddo
  ! generate the rest (L-1) block rows
  do k = 2, L
    do i = 1, N
      vals(ptrv: ptrv+N-1) = - matsB(i, 1:N, k)
      ptrv = ptrv + N
      vals(ptrv) = 1
      ptrv = ptrv + 1
      ci(ptrc: ptrc+N-1) = ind(1:N) + (k-2)*N
      ptrc = ptrc + N
      ci(ptrc) = (k-1)*N + i
      ptrc = ptrc + 1
    end do
  enddo
  ! generate the row pointer rp
  ! in matrix M, each row has exactly N+1 non-zero elements
  do i = 1, N*L+1
    rp(i) = 1 + (N+1)*(i-1)
  end do
end subroutine dense2csr
