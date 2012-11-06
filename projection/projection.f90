!This is a program that will take in a data set and give a projection of
!that data down onto a specified number of the dimensions.

program projection
  use types, only : dp
  use features, only : newunit
  implicit none

  character(len=32) :: arg1, arg2, arg3
  real(dp), dimension(:), allocatable :: val, interm
  logical, dimension(:), allocatable :: keep
  integer :: ndimns, binary_keep, modul, u, w, i, j, k, next

  !Get the bin file name from command line
  call get_command_argument(1,arg1)
  arg1=trim(arg1)
  !Find how many columns in bin file.
  call get_command_argument(2,arg2)
  arg3=trim(arg2)
  read(arg2,*) ndimns
  !Find which columns to keep.
  !This will be a binary number in the form 1010 etc where the 1 indicates we
  !want to keep that dimension in the projection and a 0 means get rid of it.
  call get_command_argument(3,arg3)
  arg3=trim(arg3)
  read(arg3,*) binary_keep

  !Make the keep array.
  allocate(keep(ndimns),val(ndimns))
  do i=1, ndimns
    modul = mod(binary_keep,10)
    if (modul == 1) then
      keep(ndimns-(i-1)) = .true.
    else
      keep(ndimns-(i-1)) = .false.
    end if
    binary_keep=binary_keep/10 !Note int div.
  end do
  allocate(interm(count(keep)))

  !Do the projection.
  open(unit=newunit(u),file=arg1,form="unformatted",status='old')
  open(unit=newunit(w),file="newfile_projection.txt")

  i=1
  do
    read(unit=u, end=999) (val(j),j=1,ndimns)
    next=1
    do k=1,ndimns
      if (keep(k)) then
        interm(next)=val(k)
        next=next+1
      end if
    end do
    write(unit=w,fmt=*) (interm(j),j=1,size(interm))
    i=i+1
  end do

999 close(u)
  close (w)


end program projection
