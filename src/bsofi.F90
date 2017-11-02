subroutine bsofi(n,L,A)
  implicit none
  integer, intent(in) :: n, L
  double precision, intent(inout) :: A(n*L,n*L)
  double precision, allocatable :: W1(:), tau(:)
  integer :: lw, info, NL
  NL = n*L
  allocate(tau(NL))
  call bsofiLwork(n,L,A,NL,tau,lw)
  allocate(W1(lw))
  call bsof(n,L,A,NL,tau,W1,lw,info)
  call bstri(n,L,A,NL,info)
  call bsoi_ormqr(n,L,A,NL,tau,W1,lw,info)   
  deallocate(W1, tau)
end subroutine bsofi

subroutine bsofiLwork( n, L, A, lda, tau, lworktotal)
  implicit none
  integer, intent(in) :: n, L, lda
  integer, intent(out) :: lworktotal
  double precision, intent(inout) :: A(n*L,n*L), tau(n*L)
  integer :: lwork, info
  double precision :: work(1)
  lwork = -1  
  lworktotal = 0
  call bsof(n,L,A,lda,tau,work,lwork,info)
  lworktotal = work(1)
  lwork = -1
  call bsoi_ormqr(n,L,A,lda,tau,work,lwork,info)
  if(lworktotal < work(1)) then
    lworktotal = work(1)
  endif   
end subroutine bsofiLwork

subroutine bsof(n, L, A, lda, tau, work, lwork, info)
  implicit none
  integer, intent(in) :: n, L, lda
  integer, intent(inout) :: lwork
  double precision, intent(inout) :: A(n*L,n*L), tau(n*L), work(1)
  integer, intent(out) :: info
  double precision :: max_work1
  integer :: i
  if(lwork == -1) then
    call dgeqrf(2*n, 2*n, A, lda, tau, work, lwork, info)
    max_work1 = work(1)
    lwork = -1
    call dormqr('L','T',2*n,n,n,A,lda,tau,A,lda,work,lwork,info)
    if(max_work1 > work(1)) then
      work(1) = max_work1
    endif
    return
  endif
  do i = 1, L-2
    call dgeqrf(2*n, n, A((i-1)*n+1,(i-1)*n+1), lda, &
                tau((i-1)*n+1:i*n), work, lwork, info)
    call dormqr('L','T',2*n,n,n,A((i-1)*n+1,(i-1)*n+1),&
                lda,tau((i-1)*n+1:i*n),A((i-1)*n+1,i*n+1),&
                lda,work,lwork,info)
    call dormqr('L','T',2*n,n,n,A((i-1)*n+1,(i-1)*n+1),lda,&
                tau((i-1)*n+1:i*n),A((i-1)*n+1,(L-1)*n+1),&
                lda,work,lwork,info) 
  enddo
  call dgeqrf(2*n,2*n,A((L-2)*n+1,(L-2)*n+1),lda,&
              tau((L-2)*n+1:L*n),work,lwork,info)
end subroutine bsof

subroutine bstri(n,L,A,lda,info)
  implicit none
  integer, intent(in) :: n, L, lda
  double precision, intent(inout) :: A(n*L,n*L)
  integer, intent(out) :: info
  integer :: NL, i
  NL = n*L
  i = L - 2
  if(L <= 2 ) then
    call dtrtri('U','N',NL,A(1,1),lda,info)
  else
    call dtrtri('U','N',3*n,A((L-3)*n+1,(L-3)*n+1),lda,info)
    if (L>3) then
      call dtrmm('R','U','N','N',n*(L-3),n,1.0d0,A((L-1)*n+1,(L-1)*n+1),lda, &
                 A(1,(L-1)*n+1),lda)
      do i = L - 3, 1, -1
        call dtrtri('U','N',n,A((i-1)*n+1,(i-1)*n+1),lda,info)
        call dtrmm('L','U','N','N',n,n,-1.0d0,A((i-1)*n+1,(i-1)*n+1),&
                   lda,A((i-1)*n+1,(L-1)*n+1),lda)
        call dtrmm('L','U','N','N',n,n,-1.0d0,A((i-1)*n+1,(i-1)*n+1),&
                   lda,A((i-1)*n+1, i*n+1), lda)
        call dgemm('N','N',n,n*(L-i-1),n,1.0d0,A((i-1)*n+1,i*n+1),lda,&
                   A(i*n+1,(i+1)*n+1),lda,1.0d0,A((i-1)*n+1,(i+1)*n+1),lda)
        call dtrmm('R','U','N','N',n,n,1.0d0,A(i*n+1,i*n+1),lda,&
                   A((i-1)*n+1,i*n+1),lda)
      enddo
    endif
  endif
end subroutine bstri

subroutine bsoi_ormqr(n, L, A, lda, tau, work, lwork, info)
  implicit none
  integer, intent(in) :: n, L ,lda, lwork
  integer, intent(out) :: info
  double precision, intent(inout) :: A(n*L,n*L), tau(n*L), work(1)   
  integer :: NL, i, ldq, j, k
  double precision Q(2*n,2*n)
  NL = n*L
  ldq = 2*n
  if(lwork == -1) then
    call dormqr('R','T',NL,2*n,2*n,work,ldq,tau,A,lda,work,lwork,info)
    work(1) = work(1)+4*n*n
    return
  endif
  call dlacpyzero('L',2*n,2*n,A((L-2)*n+1:L*n,(L-2)*n+1:L*n),&
                  lda,Q,ldq)
  call dormqr('R','T',NL,2*n,2*n,Q,ldq,tau((L-2)*n+1:L*n),&
              A(1,(L-2)*n+1),lda,work,lwork,info)

  do k = L-2, 1, -1
    call dlacpyzero('L',2*n,n,A((k-1)*n+1:(k+1)*n,(k-1)*n+1:k*n),&
                   lda,Q,ldq)
    call dormqr('R','T',NL,2*n,n,Q,ldq,tau((k-1)*n+1:k*n),&
                A(1,(k-1)*n+1),lda,work,lwork,info)
  enddo
end subroutine bsoi_ormqr

subroutine dlacpyzero(uplo, m, n, A, lda, B, ldb)
  implicit none
  character, intent(in) :: uplo
  integer, intent(in) :: m, n, lda, ldb
  double precision, intent(inout) :: A(m,n), B(m,n)
  integer :: i, j
  select case (uplo)
    case('L')
    do j = 1,n
      do i = j+1, m
        B(i,j) = A(i,j)
        A(i,j) = 0
      enddo
    enddo   
    case('U')
    do j = 1,n
      do i = 1,j
        B(i,j) = A(i,j)
        A(i,j) = 0
      enddo
    enddo
    case default
    do j = 1,n
      do i = 1,m
        B(i,j) = A(i,j)
        A(i,j) = 0
      enddo
    enddo
  end select  
end subroutine dlacpyzero
