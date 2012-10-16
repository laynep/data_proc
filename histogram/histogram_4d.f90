!Program to make a histogram of 4-dimensional data, for each dimension.  The
!number of desired bins is specified in namelist &sample.  The files to load
!from are totalsucc.bin and totalfail.bin.  If we want to only perform the
!histogram binning on the first file, then specify both=.false. in the namelist.

!NOTE: this has a check when reading files to see if they are ~0, which was not
!allowed when this program was originally implemented.  If we want these entries
!to be OK, then omit the small_value_check routine in the loop.

program histogram
  use types, only :dp
  use features, only : newunit
  implicit none

	real(dp), dimension(:), allocatable :: h1, h2, h3, h4
	real(dp) :: tol, minim(4)
	real(dp), dimension(4) :: del
	real(dp), dimension(:), allocatable :: v1, v2, v3, v4
	integer :: bins, succ_l, succ_w, fail_l, fail_w, small_l, small_w
  integer :: big_l, big_w
	integer :: i, j, u, test, err, w
	logical, dimension(:), allocatable :: mask1, mask2, mask3, mask4
  logical :: looping, both, core, toosmall, badvalue

	!Read dimns of succ table, numb of bins want in histog.
	namelist/ sample / succ_l, succ_w, fail_l, fail_w, bins, both, core,&
    &small_l, small_w, big_l, big_w

	!Finds dimension for success file.
	print*,"Reading namelist"
	open(unit=newunit(u),file="params_hist.txt",status="old",delim="apostrophe")
	read(u,nml=sample)
	close(u)

 	!Open the success file and read into vectors.
  print*,"Reading success file."
  allocate(h1(bins),h2(bins),h3(bins),h4(bins))
  open(unit=newunit(u), file="totalsucc.bin", form="unformatted")

  !Allocate working arrays.
 	allocate(v1(succ_l),v2(succ_l),v3(succ_l),v4(succ_l))
  allocate(mask1(succ_l),mask2(succ_l),mask3(succ_l),mask4(succ_l))
  do i=1,succ_l
 		read(u,iostat=err),v1(i),v2(i),v3(i),v4(i)
    call small_value_check(toosmall)
    badvalue=(err>0 .and. i .ne. succ_l) .or. toosmall
    if (badvalue) then
      print*,"Reading encountered a bad value at ", i, "success"
      stop
    end if
  end do

	!Find bin size for each dimension.
  call find_bin_size()
	print*,"Step size", del

  !Set lowest point.
  minim(1)=minval(v1)
  minim(2)=minval(v2)
  minim(3)=minval(v3)
  minim(4)=minval(v4)

  !Calculate the histogram.
	print*, "Calculating histogram for success table"
	call calc_histogram()

	!Print histogram data
	print*,"Printing data"
  call print_histogram(1)

  if (core) then
    !Repeat for the core samples.
    deallocate(h1,h2,h3,h4)
    deallocate(v1,v2,v3,v4)
    deallocate(mask1,mask2, mask3, mask4)
    print*,"Reading small corepoint files."
    allocate(h1(bins),h2(bins),h3(bins),h4(bins))
	  open(unit=newunit(u), file="smallcorepoints.bin", form="unformatted")
    !Allocate working arrays.
  	allocate(v1(small_l),v2(small_l),v3(small_l),v4(small_l))
   	allocate(mask1(small_l),mask2(small_l),mask3(small_l),mask4(small_l))

    do i=1,small_l
   		read(u,iostat=err),v1(i),v2(i),v3(i),v4(i)
      call small_value_check(toosmall)
      badvalue=(err>0 .and. i .ne. small_l) .or. toosmall
      if (badvalue) then
        print*,"Reading encountered a bad value at ", i, "small core"
        stop
      end if
    end do
    close(u)

    !Calculate the histogram.
  	print*, "Calculating histogram for small core table."
  	call calc_histogram()

  	!Print histogram data
  	print*,"Printing data"
    call print_histogram(3)

    !Repeat for the big core points.
    deallocate(h1,h2,h3,h4)
    deallocate(v1,v2,v3,v4)
    deallocate(mask1,mask2, mask3, mask4)
    print*,"Reading big corepoint files."
    allocate(h1(bins),h2(bins),h3(bins),h4(bins))
	  open(unit=newunit(u), file="bigcorepoints.bin", form="unformatted")
    !Allocate working arrays.
  	allocate(v1(big_l),v2(big_l),v3(big_l),v4(big_l))
   	allocate(mask1(big_l),mask2(big_l),mask3(big_l),mask4(big_l))

    do i=1,big_l
   		read(u,iostat=err),v1(i),v2(i),v3(i),v4(i)
      call small_value_check(toosmall)
      badvalue=(err>0 .and. i .ne. big_l) .or. toosmall
      if (badvalue) then
        print*,"Reading encountered a bad value at ", i, "big core"
        stop
      end if
    end do
    close(u)

    !Calculate the histogram.
  	print*, "Calculating histogram for big core table."
  	call calc_histogram()

  	!Print histogram data
  	print*,"Printing data"
    call print_histogram(4)

  end if

  if (both) then
    !Repeat for both success and fail.
    deallocate(v1,v2,v3,v4,mask1,mask2,mask3,mask4)
    print*,"Reading fail file."

    !Allocate working arrays.
    allocate(v1(succ_l+fail_l),v2(succ_l+fail_l),v3(succ_l+fail_l),v4(succ_l+fail_l))
    allocate(mask1(succ_l+fail_l),mask2(succ_l+fail_l),mask3(succ_l+fail_l),mask4(succ_l+fail_l))

    !Load the success part of vectors.
	  open(unit=newunit(u), file="totalsucc.bin", form="unformatted")
    do i=1,succ_l
    	read(u,iostat=err),v1(i),v2(i),v3(i),v4(i)
    end do
    close(u)

  	open(unit=newunit(w), file="totalfail.bin", form="unformatted")
    do i=succ_l+1,succ_l+fail_l
    	read(w,iostat=err),v1(i),v2(i),v3(i),v4(i)
      call small_value_check(toosmall)
      badvalue=(err>0 .and. i .ne. succ_l+fail_l) .or. toosmall
      if (badvalue) then
        print*,"Reading encountered a bad value at ", i, "total"
        stop
      end if
    end do
    close(w)

    !Calculate the histogram.
  	print*, "Calculating histogram for both tables"
  	call calc_histogram()

  	!Print histogram data
  	print*,"Printing data"
    call print_histogram(2)

  end if

  contains

    subroutine find_bin_size()

      del(1)=(maxval(v1)-minval(v1))/dble(bins)
    	del(2)=(maxval(v2)-minval(v2))/dble(bins)
    	del(3)=(maxval(v3)-minval(v3))/dble(bins)
    	del(4)=(maxval(v4)-minval(v4))/dble(bins)

    end subroutine find_bin_size

    subroutine calc_histogram()

      do i=0,bins-1
    		mask1= ((v1.ge.minim(1)+(i*del(1))) .and. (v1<minim(1)+((i+1)*del(1))))
    		mask2= ((v2.ge.minim(2)+(i*del(2))) .and. (v2<minim(2)+((i+1)*del(2))))
    		mask3= ((v3.ge.minim(3)+(i*del(3))) .and. (v3<minim(3)+((i+1)*del(3))))
    		mask4= ((v4.ge.minim(4)+(i*del(4))) .and. (v4<minim(4)+((i+1)*del(4))))
    		h1(i+1)=count(mask1)
    		h2(i+1)=count(mask2)
    		h3(i+1)=count(mask3)
    		h4(i+1)=count(mask4)
    	end do

    end subroutine calc_histogram

    subroutine print_histogram(loop)

      integer, intent(in) :: loop

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
        open(unit=1,file="data_total_histpsi.1")
        open(unit=2,file="data_total_histphi.2")
        open(unit=3,file="data_total_histpsi_dot.3")
        open(unit=4,file="data_total_histphi_dot.4")
      end if
    	do i=1,bins
    		write(1,fmt=*), minim(1)+del(1)*i, h1(i)
    		write(2,fmt=*), minim(2)+del(2)*i, h2(i)
    		write(3,fmt=*), minim(3)+del(3)*i, h3(i)
    		write(4,fmt=*), minim(4)+del(4)*i, h4(i)
    	end do
      close(1)
      close(2)
      close(3)
      close(4)

    end subroutine print_histogram

    subroutine small_value_check(checking)
      logical, intent(out) :: checking

      tol=1e0_dp

      checking=.false.
      if (abs(v1(i))<tol .and. abs(v2(i))<tol .and. abs(v3(i))<tol .and. &
        &abs(v4(i))<tol) then
        print*,"Sanity check: v's are too small..."
        print*,v1(i),v2(i),v3(i),v4(i)
        checking=.true.
      end if


    end subroutine small_value_check

end program histogram


