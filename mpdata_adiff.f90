#define USE_MPI
MODULE mpdata_adiff
!
    use mpi
    implicit none


CONTAINS

#ifndef USE_MPI

! General implementation funcions.
!
!***********************************************************************
SUBROUTINE mpdata_adiff_tile (LBi, UBi, LBj, UBj, oHz, Huon, Hvom, W, Ta, Uind, Dn, Dm, Ua)
!***********************************************************************
!
    implicit none
!
!  Imported variable declarations.
!
    integer, intent(in) :: LBi, UBi, LBj, UBj
!
    real*8, intent(in) :: oHz(LBi:,LBj:,:)
    real*8, intent(in) :: Huon(LBi:,LBj:,:)
    real*8, intent(in) :: Hvom(LBi:,LBj:,:)
    real*8, intent(in) :: W(LBi:,LBj:,:)
    real*8, intent(in) :: Ta(LBi:,LBj:,:)
    integer, intent(in) :: Uind(:)
    real*8, intent(in) :: Dn(:)
    real*8, intent(in) :: Dm(:)

    real*8, intent(out) :: Ua(LBi:,LBj:,:)

!
!  Local variable declarations.
!
    integer :: Istr, Iend, Jstr, Jend

    Istr = LBi+2
    Iend = UBi-2
    Jstr = LBj+2
    Jend = UBj-2
!
    call stencil (Istr, Iend, Jstr, Jend, LBi, UBi, LBj, UBj, oHz, Huon, Hvom, W, Ta, Uind, Dn, Dm, Ua)

    RETURN
END SUBROUTINE mpdata_adiff_tile

!
!***********************************************************************
SUBROUTINE stencil (Istr, Iend, Jstr, Jend, LBi, UBi, LBj, UBj, oHz, Huon, Hvom, W, Ta, Uind, Dn, Dm, Ua)
!***********************************************************************
!  Imported variable declarations.
!
    USE mod_data, ONLY: N, ND, dt
    integer, intent(in) :: Istr, Iend, Jstr, Jend
    integer, intent(in) :: LBi, UBi, LBj, UBj

    real*8, intent(in) :: oHz(LBi:,LBj:,:)
    real*8, intent(in) :: Huon(LBi:,LBj:,:)
    real*8, intent(in) :: Hvom(LBi:,LBj:,:)
    real*8, intent(in) :: W(LBi:,LBj:,:)
    real*8, intent(in) :: Ta(LBi:,LBj:,:)
    integer, intent(in) :: Uind(:)
    real*8, intent(in) :: Dn(:)
    real*8, intent(in) :: Dm(:)

    real*8, intent(out) :: Ua(LBi:,LBj:,:)
!
!  Local variable declarations.
!
    integer :: i, j, k, l

    real*8, parameter :: eps = 1.0E-14

    real*8 :: A, B, Um, Vm, X, Y
    real*8 :: AA, BB, AB
    real*8 :: XX, YY, XY
    real*8 :: sig_alfa, sig_beta, sig_gama
    real*8 :: sig_a, sig_b, sig_c

    real*8, dimension(LBi:UBi,N) :: C
    real*8, dimension(LBi:UBi,N) :: Wm

    DO j=Jstr,Jend
        k=1
        DO i=Istr,Iend
            C(i,k) =0.25*((Ta(i,j,k+1)-Ta(i,j,k))+(Ta(i-1,j,k+1)-Ta(i-1,j,k)))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
            Wm(i,k)=0.25*dt*(W(i-1,j,k)+W(i,j,k))
        END DO
        DO k=2,N-1
        DO i=Istr,Iend
            C(i,k) =0.0625*((Ta(i,j,k+1)-Ta(i,j,k))+(Ta(i,j,k)-Ta(i,j,k-1))+ &
                         &  (Ta(i-1,j,k+1)-Ta(i-1,j,k))+(Ta(i-1,j,k)-Ta(i-1,j,k-1)))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
            Wm(i,k)=0.25*dt*((W(i-1,j,k-1)+ W(i-1,j,k))+(W(i,j,k)+ W(i,j,k-1)))
        END DO
        END DO
        k=N
        DO i=Istr,Iend
            C(i,k) =0.25*  ((Ta(i,j,k)-Ta(i,j,k-1))+(Ta(i-1,j,k)-Ta(i-1,j,k-1)))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
            Wm(i,k)=0.25*dt*( W(i-1,j,k-1)+ W(i,j,k-1))
        END DO

        DO k=1,N
        DO i=Istr,Iend
            IF ((Ta(i-1,j,k).le.0.0).or.(Ta(i,j,k).le.0.0)) THEN
                Ua(i,j,k)=0.0
            ELSE
                A=(Ta(i,j,k)-Ta(i-1,j,k))/(Ta(i,j,k)+Ta(i-1,j,k)+eps)
                B=0.03125*((Ta(i,j+1,k)-Ta(i,j,k))+(Ta(i,j,k)-Ta(i,j-1,k))+ &
                         & (Ta(i-1,j+1,k)-Ta(i-1,j,k))+(Ta(i-1,j,k)-Ta(i-1,j-1,k)))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
!
                Um=0.125*dt*Huon(i,j,k)*(oHz(i-1,j,k)+oHz(i,j,k))
                Vm=0.03125*dt*(Hvom(i-1,j,k)*(oHz(i-1,j,k)+oHz(i-1,j-1,k))+Hvom(i-1,j+1,k)*(oHz(i-1,j+1,k)+ &
                         &     oHz(i-1,j,k))+Hvom(i,j,k)*(oHz(i,j,k)+oHz(i,j-1,k))+Hvom(i,j+1,k)*(oHz(i,j+1,k)+oHz(i,j,k)))
!
                X=(ABS(Um)-Um*Um)*A-B*Um*Vm-C(i,k)*Um*Wm(i,k)
                Y=(ABS(Vm)-Vm*Vm)*B-A*Um*Vm-C(i,k)*Vm*Wm(i,k)
!
                AA=A*A
                BB=B*B
                AB=A*B

                XX=X*X
                YY=Y*Y
                XY=X*Y

                sig_alfa=1.0/(1.0-ABS(A)+eps)
                sig_beta=-A/((1.0-ABS(A))*(1.0-AA)+eps)
                sig_gama=2.0*ABS(AA*A)/((1.0-ABS(A))*(1.0-AA)*(1.0-ABS(AA*A))+eps)
                sig_a=-B/((1.0-ABS(A))*(1.0-ABS(AB))+eps)
                sig_b=AB/((1.0-ABS(A))*(1.0-AA*ABS(B))+eps)*(ABS(B)/(1.0-ABS(AB)+eps)+2.0*A/(1.0-AA+eps))
                sig_c=ABS(A)*BB/((1.0-ABS(A))*(1.0-BB*ABS(A))*(1.0-ABS(AB))+eps)

                Ua(i,j,k)=sig_alfa*X+sig_beta*XX+sig_gama*XX*X+sig_a*XY+sig_b*XX*Y+sig_c*X*YY
!
!  Limit by physical velocity.
!
                Ua(i,j,k)=MIN(ABS(Ua(i,j,k)), ABS(Um)*SIGN(1.0,Ua(i,j,k)))

!  Further value fixing.
                DO l=1, ND
                    IF(Uind(l).eq.i) THEN
                        Ua(i,j,k)=Ua(i,j,k)+Ua(i,j,k)**Dn(l)*Um*Wm(i,k)+ABS(SIN(Dm(l))*Vm*C(i,k)*Wm(i,k))
                    ENDIF
                END DO
            END IF
        END DO
        END DO
    END DO

    RETURN
END SUBROUTINE stencil

#else
! MPI version functions. Not implemented at all.
!
!***********************************************************************
SUBROUTINE distribute_init_data (myrank,ierr)
!***********************************************************************
! Distribute the data to all processes from proc0.
! Would not be counted into calculation time.
! Don't do other work inside this subroutine.
  USE mod_data
  integer, intent(in) :: myrank,ierr
  integer :: proc_nums,start,end,buffsize,start1,end1!
  integer,pointer :: start_list(:)
  integer,pointer :: end_list(:)
  integer :: tasknums
  integer :: status(MPI_STATUS_SIZE)
  integer :: request
  character,pointer :: mat_buffer(:)
  character,pointer :: arr_buffer(:)
  integer :: l,m,positions
  integer :: i,j,k,proc,position
  call MPI_COMM_SIZE(MPI_COMM_WORLD,proc_nums,ierr)
  IF(myrank.ne.0) THEN
    N = 180
    ND = 2560
    LBi = 1
    UBi = 320
    LBj = 1
    UBj = 1200
    allocate (MYDATA % Uind(1:ND))
    allocate (MYDATA % Dn(1:ND))
    allocate (MYDATA % Dm(1:ND))
  END IF
  tasknums=int(((UBj-LBj)+1)/proc_nums)
  buffsize=(UBi-LBi+1)*(tasknums+8)*N*8
  allocate (arr_buffer(ND*3*8))
  allocate (start_list(0:proc_nums-1))
  allocate (end_list(0:proc_nums-1))
  call get_start_and_end_list(tasknums,proc_nums,start_list,end_list,UBj)
  start=start_list(myrank)
  end=end_list(myrank)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  allocate (mat_buffer(buffsize))
  if(myrank.ne.0) THEN
    allocate (MYDATA % oHz(LBi:UBi,start:end,N))!+-3
    allocate (MYDATA % Huon(LBi:UBi,start:end,N))!+-2
    allocate (MYDATA % Hvom(LBi:UBi,start:end,N))!+-3
    allocate (MYDATA % W(LBi:UBi,start:end,N))!+-2
    allocate (MYDATA % Ta(LBi:UBi,start:end,N))!+-3
    allocate (MYDATA % Ua(LBi:UBi,start:end,N))
  END IF
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !打包数组数据并传输
  !打包一维数组
  if(myrank.eq.0) THEN
    write(*,*) "[DISTRIBUTE]","start to pack arrays:Uind,Dn,Dm."
    position=0
    call MPI_PACK(MYDATA % Uind(1),ND,MPI_INT,arr_buffer,ND*3*8,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(MYDATA % Dn(1),ND,MPI_REAL8,arr_buffer,ND*3*8,position,MPI_COMM_WORLD,ierr)
    call MPI_PACK(MYDATA % Dm(1),ND,MPI_REAL8,arr_buffer,ND*3*8,position,MPI_COMM_WORLD,ierr)
    DO proc=1,proc_nums-1
      call MPI_SEND(arr_buffer,position,MPI_PACKED,proc,0,MPI_COMM_WORLD,ierr)
    END DO
    deallocate(arr_buffer)
  END IF
  if(myrank.ne.0) THEN
    call MPI_IRECV(arr_buffer,ND*8*3,MPI_PACKED,0,0,MPI_COMM_WORLD,request,ierr)
    call MPI_WAIT(request,status,ierr)
    position=0
    call MPI_UNPACK(arr_buffer,ND*3*8,position,MYDATA % Uind(1),ND,MPI_INT,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(arr_buffer,ND*3*8,position,MYDATA % Dn(1),ND,MPI_REAL8,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(arr_buffer,ND*3*8,position,MYDATA % Dm(1),ND,MPI_REAL8,MPI_COMM_WORLD,ierr)
    deallocate(arr_buffer)
  END IF
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !打包矩阵并依次发送
  if(myrank.eq.0) THEN
    write(*,*) "[DISTRIBUTE]","array packs have been recieved and unpacked."
    write (*,*) "[DISTRIBUTE]","start to pack and send matrixs."
    DO i=1,proc_nums-1
      start1=start_list(i)
      end1=end_list(i)
      call pack_and_send(MYDATA % oHz,start1,end1,mat_buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize,i,3)
      call pack_and_send(MYDATA % Huon,start1,end1,mat_buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize,i,4)
      call pack_and_send(MYDATA % Hvom,start1,end1,mat_buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize,i,5)
      call pack_and_send(MYDATA % W,start1,end1,mat_buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize,i,6)
      call pack_and_send(MYDATA % Ta,start1,end1,mat_buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize,i,7)
      write(*,*) "[DISTRIBUTE]","0 have pack and send matrix:ohz,huon,hvom,w,ta to process",i
    END DO
  END IF
  if(myrank.ne.0) THEN
    start=start_list(myrank)
    end=end_list(myrank)
    call recv_unpack(MYDATA % oHz,start,end,mat_buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize,0,3,status)
    call recv_unpack(MYDATA % Huon,start,end,mat_buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize,0,4,status)
    call recv_unpack(MYDATA % Hvom,start,end,mat_buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize,0,5,status)
    call recv_unpack(MYDATA % W,start,end,mat_buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize,0,6,status)
    call recv_unpack(MYDATA % Ta,start,end,mat_buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize,0,7,status)
    write (*,*) "[DISTRIBUTE]",myrank," have recieved and unpacked matrix:ohz,huon,hvom,w,ta from process 0"
    !deallocate(mat_buffer)
    !在这释放缓存区有BUG，没弄清楚
  END IF
END SUBROUTINE
!***********************************************************************
!数据分配后调用计算的主函数
SUBROUTINE mpdata_adiff_tile_mpi (LBi, UBi, LBj, UBj, oHz, Huon, Hvom, W, Ta, Uind, Dn, Dm, Ua,myrank,ierr)
!***********************************************************************
! Distributed implementation for calculation.
implicit none

integer, intent(in) :: LBi, UBi, LBj, UBj,myrank,ierr
real*8, intent(in) :: oHz(:,:,:)
real*8, intent(in) :: Huon(:,:,:)
real*8, intent(in) :: Hvom(:,:,:)
real*8, intent(in) :: W(:,:,:)
real*8, intent(in) :: Ta(:,:,:)
integer, intent(in) :: Uind(:)
real*8, intent(in) :: Dn(:)
real*8, intent(in) :: Dm(:)
real*8, intent(out) :: Ua(:,:,:)
integer,pointer :: start_list(:)
integer,pointer :: end_list(:)
integer :: start,end,tasknums,proc_nums
integer :: Istr, Iend, Jstr, Jend
call MPI_COMM_SIZE(MPI_COMM_WORLD,proc_nums,ierr)
tasknums=int(((UBj-LBj)+1)/proc_nums)
allocate(start_list(0:proc_nums-1))
allocate(end_list(0:proc_nums-1))
call get_start_and_end_list(tasknums,proc_nums,start_list,end_list,UBj)
start=start_list(myrank)
end=end_list(myrank)
if(myrank.ne.0) THEN
  Istr = LBi+2
  Iend = UBi-2
  Jstr = start+2
  Jend = end-2
  call stencil_mpi (Istr, Iend, Jstr, Jend, LBi, UBi, LBj, UBj, oHz, Huon, Hvom, W, Ta, Uind, Dn, Dm, Ua)
ELSE
  Istr = LBi+2
  Iend = UBi-2
  Jstr = start+3
  Jend = end-2
  call stencil_mpi (Istr, Iend, Jstr, Jend, LBi, UBi, LBj, UBj, oHz, Huon, Hvom, W, Ta, Uind, Dn, Dm, Ua)
END IF
RETURN
END SUBROUTINE
!***********************************************************************
!从各个进程收集Ua的结果
SUBROUTINE gather_data (myrank,ierr,Ua,LBi, UBi, LBj, UBj,N)
!***********************************************************************
! Gather the distributed output data Ua to MYDATA%Ua proc0.
! Would not be counted into calculation time.
! Don't do other work inside this subroutine.
!fortran模块之间参数传递不太清楚，直接再定义一遍。
integer, intent(in) :: LBi, UBi, LBj, UBj,myrank,ierr,N
real*8, intent(out) :: Ua(:,:,:)
integer,pointer :: start_list(:)
integer,pointer :: end_list(:)
character,pointer :: mat_buffer(:)
integer :: status(MPI_STATUS_SIZE)
integer :: tasknums,proc_nums,proc,buffsize
integer :: start2,end2
call MPI_COMM_SIZE(MPI_COMM_WORLD,proc_nums,ierr)
tasknums=int(((UBj-LBj)+1)/proc_nums)
buffsize=(UBi-LBi+1)*(tasknums+8)*N*8
allocate(start_list(0:proc_nums-1))
allocate(end_list(0:proc_nums-1))
allocate(mat_buffer(buffsize))
call get_start_and_end_list(tasknums,proc_nums,start_list,end_list,UBj)
if(myrank.ne.0) THEN
  if(myrank.ne.1) THEN
    start2=start_list(myrank)+3
    end2=end_list(myrank)-3
    call pack_and_send(Ua,start2,end2,mat_buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize,0,myrank)
  ELSE
    start2=start_list(myrank)+3
    end2=end_list(myrank)-3
    call pack_and_send(Ua,start2,end2,mat_buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize,0,myrank)
  END IF
ELSE
  DO proc=1,proc_nums-1
    start2=start_list(proc)+3
    end2=end_list(myrank)-3
    call recv_unpack(Ua,start2,end2,mat_buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize,proc,proc,status)
  END DO
  write (*,*) "[GATHER]  gather data complite."
END IF
END SUBROUTINE
!***********************************************************************
!打包矩阵块到一个连续空间并发送
SUBROUTINE pack_and_send(matrix,start,end,buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize,dest,tag)
  integer, intent(in) :: LBi, UBi, LBj, UBj, start, end,myrank,N,buffsize,tasknums,dest,tag,ierr
  real*8, intent(inout) :: matrix(:,:,:)
  character, intent(out) :: buffer(:)
  call pack_matrix(matrix,start,end,buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize)
  call MPI_SEND(buffer,buffsize,MPI_PACKED,dest,tag,MPI_COMM_WORLD,ierr)
  RETURN
END SUBROUTINE
!***********************************************************************
!把matrix(:,start:end,:)打包到发送缓冲区
SUBROUTINE pack_matrix(matrix,start,end,buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize)
  integer, intent(in) :: LBi, UBi, LBj, UBj, start, end,myrank,N,buffsize,tasknums,ierr
  real*8, intent(inout) :: matrix(:,:,:)
  character, intent(out) :: buffer(:)
  integer :: l,m,positions
  positions=0
  DO m=1,N
    DO l=start,end
      call MPI_PACK(matrix(1,l,m),UBi-LBi+1,MPI_REAL8,buffer,buffsize,positions,MPI_COMM_WORLD,ierr)
    END DO
  END DO
  RETURN
END SUBROUTINE
!***********************************************************************
!接收buffer并解包,为了上边代码简洁一些
SUBROUTINE recv_unpack(matrix,start,end,buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize,source,tag,status)
  integer, intent(in) :: LBi, UBi, LBj, UBj, start, end,myrank, ierr,N,buffsize,tasknums,tag,source
  real*8, intent(out) :: matrix(:,:,:)
  character, intent(in) :: buffer(:)
  integer,intent(out) :: status(:)
  call MPI_RECV(buffer,buffsize,MPI_PACKED,source,tag,MPI_COMM_WORLD,status,ierr)
  call mat_unpack(matrix,start,end,buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize)
  RETURN
END SUBROUTINE
!***********************************************************************
!从一个连续空间解包数据到一个矩阵块
SUBROUTINE mat_unpack(matrix,start,end,buffer,tasknums,LBi,UBi,LBj,UBj,myrank,ierr,N,buffsize)
  integer, intent(in) :: LBi, UBi, LBj, UBj, start, end,myrank, ierr,N,buffsize,tasknums
  real*8, intent(out) :: matrix(:,:,:)
  character, intent(in) :: buffer(:)
  integer :: l,m,positions
  positions=0
  DO m=1,N
    DO l=start,end
      call MPI_UNPACK(buffer,buffsize,positions,matrix(1,l,m),UBi-LBi+1,MPI_REAL8,MPI_COMM_WORLD,ierr)
    END DO
  END DO
  RETURN
END SUBROUTINE
!***********************************************************************
!计算每个进程分配的区块索引
SUBROUTINE get_start_and_end_list(tasknums,proc_nums,start_list,end_list,UBj)
  integer, intent(in) :: tasknums,proc_nums,UBj
  integer, intent(out) :: start_list(0:),end_list(0:)
  integer :: i
  if (proc_nums.eq.2) start_list(0)=tasknums-3

  IF (proc_nums.ne.2) THEN
    start_list(0)=tasknums*(proc_nums-1)-3
    DO i=2,proc_nums-1
      start_list(i)=tasknums*(i-1)-3
      end_list(i)=tasknums*i+3
    END DO
  END IF
  end_list(0)=UBj
  start_list(1)=1
  end_list(1)=tasknums+3
END SUBROUTINE
!***********************************************************************
!MPI情况下的计算函数
SUBROUTINE stencil_mpi (Istr, Iend, Jstr, Jend, LBi, UBi, LBj, UBj, oHz, Huon, Hvom, W, Ta, Uind, Dn, Dm, Ua)
!***********************************************************************
!  Imported variable declarations.
!
    USE mod_data, ONLY: N, ND, dt
    integer, intent(in) :: Istr, Iend, Jstr, Jend
    integer, intent(in) :: LBi, UBi, LBj, UBj

    real*8, intent(in) :: oHz(:,:,:)
    real*8, intent(in) :: Huon(:,:,:)
    real*8, intent(in) :: Hvom(:,:,:)
    real*8, intent(in) :: W(:,:,:)
    real*8, intent(in) :: Ta(:,:,:)
    integer, intent(in) :: Uind(:)
    real*8, intent(in) :: Dn(:)
    real*8, intent(in) :: Dm(:)

    real*8, intent(out) :: Ua(:,:,:)
!
!  Local variable declarations.
!
    integer :: i, j, k, l

    real*8, parameter :: eps = 1.0E-14

    real*8 :: A, B, Um, Vm, X, Y
    real*8 :: AA, BB, AB
    real*8 :: XX, YY, XY
    real*8 :: sig_alfa, sig_beta, sig_gama
    real*8 :: sig_a, sig_b, sig_c

    real*8, dimension(LBi:UBi,N) :: C
    real*8, dimension(LBi:UBi,N) :: Wm

    DO j=Jstr,Jend
        k=1
        DO i=Istr,Iend
            C(i,k) =0.25*((Ta(i,j,k+1)-Ta(i,j,k))+(Ta(i-1,j,k+1)-Ta(i-1,j,k)))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
            Wm(i,k)=0.25*dt*(W(i-1,j,k)+W(i,j,k))
        END DO
        DO k=2,N-1
        DO i=Istr,Iend
            C(i,k) =0.0625*((Ta(i,j,k+1)-Ta(i,j,k))+(Ta(i,j,k)-Ta(i,j,k-1))+ &
                         &  (Ta(i-1,j,k+1)-Ta(i-1,j,k))+(Ta(i-1,j,k)-Ta(i-1,j,k-1)))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
            Wm(i,k)=0.25*dt*((W(i-1,j,k-1)+ W(i-1,j,k))+(W(i,j,k)+ W(i,j,k-1)))
        END DO
        END DO
        k=N
        DO i=Istr,Iend
            C(i,k) =0.25*  ((Ta(i,j,k)-Ta(i,j,k-1))+(Ta(i-1,j,k)-Ta(i-1,j,k-1)))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
            Wm(i,k)=0.25*dt*( W(i-1,j,k-1)+ W(i,j,k-1))
        END DO

        DO k=1,N
        DO i=Istr,Iend
            IF ((Ta(i-1,j,k).le.0.0).or.(Ta(i,j,k).le.0.0)) THEN
                Ua(i,j,k)=0.0
            ELSE
                A=(Ta(i,j,k)-Ta(i-1,j,k))/(Ta(i,j,k)+Ta(i-1,j,k)+eps)
                B=0.03125*((Ta(i,j+1,k)-Ta(i,j,k))+(Ta(i,j,k)-Ta(i,j-1,k))+ &
                         & (Ta(i-1,j+1,k)-Ta(i-1,j,k))+(Ta(i-1,j,k)-Ta(i-1,j-1,k)))/(Ta(i-1,j,k)+Ta(i,j,k)+eps)
!
                Um=0.125*dt*Huon(i,j,k)*(oHz(i-1,j,k)+oHz(i,j,k))
                Vm=0.03125*dt*(Hvom(i-1,j,k)*(oHz(i-1,j,k)+oHz(i-1,j-1,k))+Hvom(i-1,j+1,k)*(oHz(i-1,j+1,k)+ &
                         &     oHz(i-1,j,k))+Hvom(i,j,k)*(oHz(i,j,k)+oHz(i,j-1,k))+Hvom(i,j+1,k)*(oHz(i,j+1,k)+oHz(i,j,k)))
!
                X=(ABS(Um)-Um*Um)*A-B*Um*Vm-C(i,k)*Um*Wm(i,k)
                Y=(ABS(Vm)-Vm*Vm)*B-A*Um*Vm-C(i,k)*Vm*Wm(i,k)
!
                AA=A*A
                BB=B*B
                AB=A*B

                XX=X*X
                YY=Y*Y
                XY=X*Y

                sig_alfa=1.0/(1.0-ABS(A)+eps)
                sig_beta=-A/((1.0-ABS(A))*(1.0-AA)+eps)
                sig_gama=2.0*ABS(AA*A)/((1.0-ABS(A))*(1.0-AA)*(1.0-ABS(AA*A))+eps)
                sig_a=-B/((1.0-ABS(A))*(1.0-ABS(AB))+eps)
                sig_b=AB/((1.0-ABS(A))*(1.0-AA*ABS(B))+eps)*(ABS(B)/(1.0-ABS(AB)+eps)+2.0*A/(1.0-AA+eps))
                sig_c=ABS(A)*BB/((1.0-ABS(A))*(1.0-BB*ABS(A))*(1.0-ABS(AB))+eps)

                Ua(i,j,k)=sig_alfa*X+sig_beta*XX+sig_gama*XX*X+sig_a*XY+sig_b*XX*Y+sig_c*X*YY
!
!  Limit by physical velocity.
!
                Ua(i,j,k)=MIN(ABS(Ua(i,j,k)), ABS(Um)*SIGN(1.0,Ua(i,j,k)))

!  Further value fixing.
                DO l=1, ND
                    IF(Uind(l).eq.i) THEN
                        Ua(i,j,k)=Ua(i,j,k)+Ua(i,j,k)**Dn(l)*Um*Wm(i,k)+ABS(SIN(Dm(l))*Vm*C(i,k)*Wm(i,k))
                    ENDIF
                END DO
            END IF
        END DO
        END DO
    END DO

    RETURN
END SUBROUTINE stencil_mpi

#endif

END MODULE mpdata_adiff
