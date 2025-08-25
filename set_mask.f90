subroutine set_mask(mask,diag,nelem,nterm,ndep,ntdof,nem, &
                      ibcw,ibce,ibcs,ibcn)
      implicit none
      integer, intent(in) :: nelem,nterm,ndep,ntdof,nem
      integer, intent(in) :: ibcw(nem),ibce(nem),ibcs(nem),ibcn(nem)
      integer, intent(out) :: mask(ntdof,nem)
      real(8), intent(out) :: diag(ntdof,nem)
      integer :: i,j,k,ne
      diag=1.
      mask=1
      do ne=1,nelem
        if( ibcw(ne).eq.1 .or. ibcw(ne).eq.2 .or. ibcw(ne).eq.3 ) then
          do j=1,nterm
            k = (j-1)*nterm + 1
            mask((k-1)*ndep+1,ne) = 0
            mask((k-1)*ndep+2,ne) = 0
            mask((k-1)*ndep+3,ne) = 0
            mask((k-1)*ndep+4,ne) = 0
          enddo
        endif
        if( ibce(ne).eq.1 .or. ibce(ne).eq.2 .or. ibce(ne).eq.3 ) then
          do j=1,nterm
            k = (j-1)*nterm + nterm
            mask((k-1)*ndep+1,ne) = 0
            mask((k-1)*ndep+2,ne) = 0
            mask((k-1)*ndep+3,ne) = 0
            mask((k-1)*ndep+4,ne) = 0
          enddo
        endif
        if( ibcs(ne).eq.1 .or. ibcs(ne).eq.2 .or. ibcs(ne).eq.3 ) then
          do i=1,nterm
            k = i
            mask((k-1)*ndep+1,ne) = 0
            mask((k-1)*ndep+2,ne) = 0
            mask((k-1)*ndep+3,ne) = 0
            mask((k-1)*ndep+4,ne) = 0
          enddo
        endif
        if( ibcn(ne).eq.1 .or. ibcn(ne).eq.2 .or. ibcn(ne).eq.3 ) then
          do i=1,nterm
            k = (nterm-1)*nterm + i
            mask((k-1)*ndep+1,ne) = 0
            mask((k-1)*ndep+2,ne) = 0
            mask((k-1)*ndep+3,ne) = 0
            mask((k-1)*ndep+4,ne) = 0
          enddo
        endif
      enddo
      return
      end
