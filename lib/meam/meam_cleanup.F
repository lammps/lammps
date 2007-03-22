c Declaration in pair_meam.h:
c
c void meam_cleanup()
c
c Call from PairMEAM destructor
c
c meam_cleanup()
c

      subroutine meam_cleanup
      use meam_data
      implicit none

      integer dealloc_error

      deallocate(phir,STAT=dealloc_error)
      deallocate(phirar,STAT=dealloc_error)
      deallocate(phirar1,STAT=dealloc_error)
      deallocate(phirar2,STAT=dealloc_error)
      deallocate(phirar3,STAT=dealloc_error)
      deallocate(phirar4,STAT=dealloc_error)
      deallocate(phirar5,STAT=dealloc_error)
      deallocate(phirar6,STAT=dealloc_error)

      return
      end
