!Program to make a histogram of 4-dimensional data, for each dimension.  The
!number of desired bins is specified in namelist &sample.  The files to load
!from are totalsucc.bin and totalfail.bin.  If we want to only perform the
!histogram binning on the first file, then specify both=.false. in the namelist.

!NOTE: this has a check when reading files to see if they are ~0, which was not
!allowed when this program was originally implemented.  If we want these entries
!to be OK, then omit the small_value_check routine in the loop.

module mod_hist
  use types, only : dp

  real(dp), public :: energy_scale

  contains

  !NOTE: this only works when calculating the bin size for hybrid inflation,
  !with the particle physics parameters M, mu, nu, Lambda hardwired into this
  !routine.
  !In general, use subroutine find_bin_size.
  subroutine find_bin_size_hybrid(table,bins,del)
    use features, only : newunit
    implicit none

    real(dp), dimension(:,:), intent(in) :: table
    real(dp), dimension(:), intent(out) :: del
    integer, intent(in) :: bins
    integer :: i, u
    real(dp) :: lambda, M, mu, nu
    real(dp), dimension(2) :: test

    namelist /parameters/ energy_scale, lambda, m, mu, nu

    open(unit=newunit(u),file="params_hist.txt",status="old",delim="apostrophe")
  	read(u,nml=parameters)
  	close(u)

    !Area that field IC were sampled from.
    test=(/ 200e0_dp,m*sqrt(1_dp + (energy_scale/lambda)**2)/)
    del(1)=minval(test)/dble(bins)
    del(2)=200e0_dp/dble(bins)
    del(3)=2e0_dp*sqrt(2e0_dp)*energy_scale**2/dble(bins)
    del(4)=del(3)

  end subroutine find_bin_size_hybrid

  subroutine find_bin_size(table,bins,del)
    implicit none

    real(dp), dimension(:,:), intent(in) :: table
    real(dp), dimension(:), intent(out) :: del
    integer, intent(in) :: bins
    integer :: i

    do i=1,size(del)
      del(i)=(maxval(table(:,i))-minval(table(:,i)))/dble(bins)
    end do

  end subroutine find_bin_size

  subroutine calc_histogram(hist_table,table,mask,minim,del)
    implicit none

    real(dp), dimension(:,:), intent(out) :: hist_table
    real(dp), dimension(:,:), intent(in) :: table
    logical, dimension(:,:), intent(inout) :: mask
    real(dp), dimension(:), intent(in) :: minim, del
    integer :: i, j

    do i=0,size(hist_table,1)-1
      do j=1,size(mask,2)
        mask(:,j) = (table(:,j).ge.minim(j)+i*del(j)) .and.&
          & (table(:,j)<minim(j) + (i+1)*del(j))
        hist_table(i+1,j)=count(mask(:,j))
      end do
  	end do

  end subroutine calc_histogram

  !NOTE: only works in 4 dimensions.
  subroutine print_histogram_4d(loop,minim,del,histog_table)
    implicit none

    integer, intent(in) :: loop
    real(dp), dimension(:), intent(in) :: minim, del
    real(dp), dimension(:,:), intent(in) :: histog_table
    integer :: i

    !Success set.
    if (loop==1) then
      open(unit=1,file="data_histpsi.1")
      open(unit=2,file="data_histphi.2")
      open(unit=3,file="data_histpsi_dot.3")
      open(unit=4,file="data_histphi_dot.4")
    !Small core set.
    else if (loop==2) then
      open(unit=1,file="small_histpsi.1")
      open(unit=2,file="small_histphi.2")
      open(unit=3,file="small_histpsi_dot.3")
      open(unit=4,file="small_histphi_dot.4")
    !Big core set.
    else if (loop==3) then
      open(unit=1,file="big_histpsi.1")
      open(unit=2,file="big_histphi.2")
      open(unit=3,file="big_histpsi_dot.3")
      open(unit=4,file="big_histphi_dot.4")
    !Total succ + fail.
    else if (loop==4) then
      open(unit=1,file="data_fail_histpsi.1")
      open(unit=2,file="data_fail_histphi.2")
      open(unit=3,file="data_fail_histpsi_dot.3")
      open(unit=4,file="data_fail_histphi_dot.4")
    end if
  	do i=1,size(histog_table,1)
  		write(1,fmt=*), minim(1)+del(1)*i, histog_table(i,1)
  		write(2,fmt=*), minim(2)+del(2)*i, histog_table(i,2)
  		write(3,fmt=*), minim(3)+del(3)*i, histog_table(i,3)
  		write(4,fmt=*), minim(4)+del(4)*i, histog_table(i,4)
  	end do
    close(1)
    close(2)
    close(3)
    close(4)

  end subroutine print_histogram_4d


  subroutine read_from_file(table,fname,fform,length)
    use features, only : newunit
    implicit none

    integer :: u, i, err, j
    integer, intent(in) :: length
    character(len=*), intent(in) :: fname, fform
    real(dp), dimension(:,:), intent(out) :: table
    logical :: toosmall, badvalue

    open(unit=newunit(u), file=fname, form=fform)

    do i=1,length
   		read(u,iostat=err),(table(i,j),j=1,size(table,2))
      call small_value_check(toosmall,table(i,:))
      badvalue=err>0 .or. toosmall
      if (badvalue) then
        print*,"Reading encountered a bad value at ", i, fname
        stop
      end if
    end do
    close(u)

  end subroutine read_from_file

  subroutine small_value_check(checking,vect)
    implicit none

    logical, intent(out) :: checking
    logical :: badvect
    real(dp) :: tol
    real(dp), dimension(:), intent(in) :: vect
    integer :: i

    tol=1e-10_dp

    checking=.false.
    badvect=all(abs(vect(:))<tol)
    if (badvect) then
      print*,"Sanity check: v's are too small..."
      print*,(vect(i),i=1,size(vect))
      checking=.true.
    end if

  end subroutine small_value_check
end module mod_hist


!----------------------------------------------------------------------------
!----------------------------------------------------------------------------


program histogram
  use mod_hist
  use types, only :dp
  use features, only : newunit
  implicit none

	real(dp), dimension(:,:), allocatable :: hist_table
  real(dp) :: tol
	real(dp), dimension(:), allocatable :: del, minim
	real(dp), dimension(:,:), allocatable :: table
  integer :: bins, succ_l, dimn, fail_l, small_l
  integer :: big_l
	integer :: i, j, u, test, err, w
	logical, dimension(:,:), allocatable :: mask
  logical :: both, core
  character(len=100) :: succfile, failfile, smallclustfile, bigclustfile

	!Read dimns of succ table, numb of bins want in histog.
	namelist/ sample / succ_l, fail_l, small_l, big_l, &
   & bins, both, core, dimn
  namelist/ filenames / succfile, failfile, smallclustfile, bigclustfile

	!Finds dimension for success file.
	print*,"Reading namelist"
	open(unit=newunit(u),file="params_hist.txt",status="old",delim="apostrophe")
	read(u,nml=sample)
  read(u,nml=filenames)
	close(u)

  !Allocate working arrays.
  allocate(hist_table(bins,dimn))
 	allocate(table(succ_l,dimn),mask(succ_l,dimn))

 	!Open the success file and read into vectors.
  print*,"Reading success file."
  call read_from_file(table,succfile,"unformatted",succ_l)

	!Find bin size for each dimension.
  allocate(minim(dimn),del(dimn))
  call find_bin_size_hybrid(table,bins,del)
	print*,"Step size", del

  !Set lowest point.
  !NOTE: this is only good for hybrid.  For other uses, uncomment below.
  minim(1)=0e0_dp
  minim(2)=minim(1)
  minim(3)=-sqrt(2e0_dp)*energy_scale**2
  minim(4)=minim(3)
 ! do i=1, size(minim)
 !   minim(i)=minval(table(:,i))
 ! end do

  !Calculate the histogram.
	print*, "Calculating histogram for success table"
	call calc_histogram(hist_table,table,mask,minim,del)

	!Print histogram data
	print*,"Printing data"
  call print_histogram_4d(1,minim,del,hist_table)

  !Repeat for the core samples.
  if (core) then

    !Allocate working arrays.
    deallocate(hist_table, table, mask)
    allocate(hist_table(bins,dimn))
   	allocate(table(small_l,dimn),mask(small_l,dimn))

    print*,"Reading small corepoint files."
    call read_from_file(table,smallclustfile,"unformatted",small_l)

    !Calculate the histogram.
  	print*, "Calculating histogram for small core table."
  	call calc_histogram(hist_table,table,mask,minim,del)

  	!Print histogram data
  	print*,"Printing data"
    call print_histogram_4d(2,minim,del,hist_table)

    !Repeat for the big core points.
    deallocate(hist_table, table, mask)

    !Allocate working arrays.
    allocate(hist_table(bins,dimn))
   	allocate(table(big_l,dimn),mask(big_l,dimn))

    print*,"Reading big corepoint files."
    call read_from_file(table,bigclustfile,"unformatted",big_l)

    !Calculate the histogram.
  	print*, "Calculating histogram for big core table."
  	call calc_histogram(hist_table,table,mask,minim,del)

  	!Print histogram data
  	print*,"Printing data"
    call print_histogram_4d(3,minim,del,hist_table)

  end if

  !Repeat for fail.
  deallocate(hist_table, table, mask)

  !Allocate working arrays.
  allocate(hist_table(bins,dimn))
 	allocate(table(fail_l,dimn),mask(fail_l,dimn))

  print*,"Reading fail file."
  call read_from_file(table,failfile,"unformatted",fail_l)

  !Calculate the histogram.
  print*, "Calculating histogram for fail table."
  call calc_histogram(hist_table,table,mask,minim,del)

  !Print histogram data
  print*,"Printing data"
  call print_histogram_4d(4,minim,del,hist_table)

end program histogram


