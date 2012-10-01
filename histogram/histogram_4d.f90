program histogram
  use types, only :dp 
  use features, only : newunit
  implicit none

	real(dp), dimension(:), allocatable :: h1, h2, h3, h4
	real(dp) :: start
	real(dp), dimension(4) :: del
	real(dp), dimension(:), allocatable :: v1, v2, v3, v4
	integer :: bins, succ_l, succ_w
	integer :: i, j, u
	logical, dimension(:), allocatable :: mask1, mask2, mask3, mask4


	!Read dimns of succ table, numb of bins want in histog.
	namelist/ sample / succ_l, succ_w, bins

	!Finds dimension for success file.
	print*,"Reading success table"
	open(unit=newunit(u),file="params_hist.txt",status="old",delim="apostrophe")
	read(u,nml=sample)
	close(u)
	!Allocate working arrays.
	allocate(v1(succ_l),v2(succ_l),v3(succ_l),v4(succ_l))
	allocate(h1(bins),h2(bins),h3(bins),h4(bins))
	allocate(mask1(succ_l),mask2(succ_l),mask3(succ_l),mask4(succ_l))
	!Open the success file and read into vectors.
	open(unit=newunit(u), file="everything.bin", form="unformatted")
	do i=1,succ_l
		read(u),v1(i),v2(i),v3(i),v4(i)
	end do
	close(u)

	!Find bin size for each dimension.
	del(1)=(maxval(v1)-minval(v1))/dble(bins)
	del(2)=(maxval(v2)-minval(v2))/dble(bins)
	del(3)=(maxval(v3)-minval(v3))/dble(bins)
	del(4)=(maxval(v4)-minval(v4))/dble(bins)
	print*,"Step size", del
	print*, "Calculating histogram"
	do i=0,bins-1
		mask1= ((v1>minval(v1)+(i*del(1))) .and. (v1<minval(v1)+((i+1)*del(1))))
		mask2= ((v2>minval(v2)+(i*del(2))) .and. (v2<minval(v2)+((i+1)*del(2))))
		mask3= ((v3>minval(v3)+(i*del(3))) .and. (v3<minval(v3)+((i+1)*del(3))))
		mask4= ((v4>minval(v4)+(i*del(4))) .and. (v4<minval(v4)+((i+1)*del(4))))
		h1(i+1)=count(mask1)
		h2(i+1)=count(mask2)
		h3(i+1)=count(mask3)
		h4(i+1)=count(mask4)
	end do

	!Print histogram data
	print*,"Printing data"
  open(unit=1,file="data_histpsi.1")
  open(unit=2,file="data_histphi.2")
  open(unit=3,file="data_histpsi_dot.3")
  open(unit=4,file="data_histphi_dot.4")
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

end program histogram
