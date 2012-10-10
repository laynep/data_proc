!Program to make a histogram of 4-dimensional data, for each dimension.  The
!number of desired bins is specified in namelist &sample.  The files to load
!from are totalsucc.bin and totalfail.bin.  If we want to only perform the
!histogram binning on the first file, then specify both=.false. in the namelist.

program histogram
  use types, only :dp
  use features, only : newunit
  implicit none

	real(dp), dimension(:), allocatable :: h1, h2, h3, h4
	real(dp) :: start
	real(dp), dimension(4) :: del
	real(dp), dimension(:), allocatable :: v1, v2, v3, v4
	integer :: bins, succ_l, succ_w, fail_l, fail_w
	integer :: i, j, u, test, err, w
	logical, dimension(:), allocatable :: mask1, mask2, mask3, mask4
  logical :: looping, both

	!Read dimns of succ table, numb of bins want in histog.
	namelist/ sample / succ_l, succ_w, fail_l, fail_w, bins, both

	!Finds dimension for success file.
	print*,"Reading namelist"
	open(unit=newunit(u),file="params_hist.txt",status="old",delim="apostrophe")
	read(u,nml=sample)
	close(u)

 	!Open the success file and read into vectors.
  print*,"Reading success file."
  allocate(h1(bins),h2(bins),h3(bins),h4(bins))
  do1:  do
	  open(unit=newunit(u), file="totalsucc.bin", form="unformatted")
    !Allocate working arrays.
  	allocate(v1(succ_l),v2(succ_l),v3(succ_l),v4(succ_l))
  	allocate(mask1(succ_l),mask2(succ_l),mask3(succ_l),mask4(succ_l))

    looping=.false.
    do2: do i=1,succ_l
  		read(u,iostat=err),v1(i),v2(i),v3(i),v4(i)
      if (err==5001 .and. i .ne. succ_l) then
        looping=.true.
        succ_l=i
        exit do2
      end if
      !if(mod(i,100000)==0) print*,i, err
    end do do2
    close(u)
    if (looping) print*,"Success set not declared properly. Reloading..."
    if (.not. looping) then
      exit do1
    else
      deallocate(v1,v2,v3,v4,mask1,mask2,mask3,mask4)
      cycle
    end if
  end do do1

	!Find bin size for each dimension.
  call find_bin_size()
	print*,"Step size", del

  !Calculate the histogram.
	print*, "Calculating histogram for success table"
	call calc_histogram()

	!Print histogram data
	print*,"Printing data"
  call print_histogram(1)

  if (both) then
    !Repeat for both success and fail.
    deallocate(v1,v2,v3,v4,mask1,mask2,mask3,mask4)
    print*,"Reading fail file."

    do3:  do
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
      looping=.false.
      do4: do i=succ_l+1,succ_l+fail_l-1
    		read(w,iostat=err),v1(i),v2(i),v3(i),v4(i)

        if (err==5001 .and. i .ne. succ_l+fail_l) then
          looping=.true.
          fail_l=i-succ_l-1
          exit do4
        end if
      end do do4
      close(w)

      if (.not. looping) then
        exit do3
      else
        print*,"Fail set not declared properly. Reloading..."
        deallocate(v1,v2,v3,v4,mask1,mask2,mask3,mask4)
        cycle
      end if
    end do do3

    !Find bin size for each dimension.
    call find_bin_size()
  	print*,"Step size", del

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
    		mask1= ((v1.ge.minval(v1)+(i*del(1))) .and. (v1<minval(v1)+((i+1)*del(1))))
    		mask2= ((v2.ge.minval(v2)+(i*del(2))) .and. (v2<minval(v2)+((i+1)*del(2))))
    		mask3= ((v3.ge.minval(v3)+(i*del(3))) .and. (v3<minval(v3)+((i+1)*del(3))))
    		mask4= ((v4.ge.minval(v4)+(i*del(4))) .and. (v4<minval(v4)+((i+1)*del(4))))
    		h1(i+1)=count(mask1)
    		h2(i+1)=count(mask2)
    		h3(i+1)=count(mask3)
    		h4(i+1)=count(mask4)
    	end do

    end subroutine calc_histogram

    subroutine print_histogram(loop)

      integer, intent(in) :: loop

      if (loop==1) then
        open(unit=1,file="data_histpsi.1")
        open(unit=2,file="data_histphi.2")
        open(unit=3,file="data_histpsi_dot.3")
        open(unit=4,file="data_histphi_dot.4")
      else
        open(unit=1,file="data_total_histpsi.1")
        open(unit=2,file="data_total_histphi.2")
        open(unit=3,file="data_total_histpsi_dot.3")
        open(unit=4,file="data_total_histphi_dot.4")
      end if
    	do i=1,bins
    		write(1,fmt=*), minval(v1)+del(1)*i, h1(i)
    		write(2,fmt=*), minval(v2)+del(2)*i, h2(i)
    		write(3,fmt=*), minval(v3)+del(3)*i, h3(i)
    		write(4,fmt=*), minval(v4)+del(4)*i, h4(i)
    	end do
      close(1)
      close(2)
      close(3)
      close(4)

    end subroutine print_histogram

end program histogram


