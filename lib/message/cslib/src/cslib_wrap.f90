! ISO_C_binding wrapper on CSlib C interface

module cslib_wrap

interface
  subroutine cslib_open_fortran(csflag,mode,str,pcomm,ptr) bind(c)
    use iso_c_binding
    integer(c_int), value :: csflag
    character(c_char) :: mode(*),str(*)
    type(c_ptr), value :: pcomm
    type(c_ptr) :: ptr
  end subroutine cslib_open_fortran

  subroutine cslib_open_fortran_mpi_one(csflag,mode,pboth,pcomm,ptr) bind(c)
    use iso_c_binding
    integer(c_int), value :: csflag
    character(c_char) :: mode(*)
    type(c_ptr), value :: pboth,pcomm
    type(c_ptr) :: ptr
  end subroutine cslib_open_fortran_mpi_one

  subroutine cslib_close(ptr) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
  end subroutine cslib_close

  subroutine cslib_send(ptr,msgID,nfield) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    integer(c_int), value :: msgID,nfield
  end subroutine cslib_send

  subroutine cslib_pack_int(ptr,id,value) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    integer(c_int), value :: id
    integer(c_int), value :: value
  end subroutine cslib_pack_int

  subroutine cslib_pack_int64(ptr,id,value) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    integer(c_int), value :: id
    integer(c_int64_t), value :: value
  end subroutine cslib_pack_int64

  subroutine cslib_pack_float(ptr,id,value) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    integer(c_int), value :: id
    real(c_float), value :: value
  end subroutine cslib_pack_float

  subroutine cslib_pack_double(ptr,id,value) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    integer(c_int), value :: id
    real(c_double), value :: value
  end subroutine cslib_pack_double

  subroutine cslib_pack_string(ptr,id,value) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    integer(c_int), value :: id
    character(c_char) :: value(*)
  end subroutine cslib_pack_string

  subroutine cslib_pack(ptr,id,ftype,flen,data) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    integer(c_int), value :: id,ftype,flen
    type(c_ptr), value :: data
  end subroutine cslib_pack

  subroutine cslib_pack_parallel(ptr,id,ftype,nlocal,ids,nper,data) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    integer(c_int), value :: id,ftype,nlocal,nper
    type(c_ptr), value :: ids,data
  end subroutine cslib_pack_parallel

  function cslib_recv(ptr,nfield,fieldID,fieldtype,fieldlen) bind(c)
    use iso_c_binding
    integer(c_int) :: cslib_recv
    type(c_ptr), value :: ptr
    integer(c_int) :: nfield
    type(c_ptr) :: fieldID,fieldtype,fieldlen
  end function cslib_recv

  function cslib_unpack_int(ptr,id) bind(c)
    use iso_c_binding
    integer(c_int) :: cslib_unpack_int
    type(c_ptr), value :: ptr
    integer(c_int), value :: id
  end function cslib_unpack_int

  function cslib_unpack_int64(ptr,id) bind(c)
    use iso_c_binding
    integer(c_int64_t) :: cslib_unpack_int64
    type(c_ptr), value :: ptr
    integer(c_int), value :: id
  end function cslib_unpack_int64

  function cslib_unpack_float(ptr,id) bind(c)
    use iso_c_binding
    real(c_float) :: cslib_unpack_float
    type(c_ptr), value :: ptr
    integer(c_int), value :: id
  end function cslib_unpack_float

  function cslib_unpack_double(ptr,id) bind(c)
    use iso_c_binding
    real(c_double) :: cslib_unpack_double
    type(c_ptr), value :: ptr
    integer(c_int), value :: id
  end function cslib_unpack_double

  function cslib_unpack_string(ptr,id) bind(c)
    use iso_c_binding
    type(c_ptr) :: cslib_unpack_string
    type(c_ptr), value :: ptr
    integer(c_int), value :: id
  end function cslib_unpack_string

  function cslib_unpack(ptr,id) bind(c)
    use iso_c_binding
    type(c_ptr) :: cslib_unpack
    type(c_ptr), value :: ptr
    integer(c_int), value :: id
  end function cslib_unpack

  subroutine cslib_unpack_parallel(ptr,id,nlocal,ids,nper,data) bind(c)
    use iso_c_binding
    type(c_ptr), value :: ptr
    integer(c_int), value :: id,nlocal,nper
    type(c_ptr), value :: ids,data
  end subroutine cslib_unpack_parallel

  function cslib_extract(ptr,flag) bind(c)
    use iso_c_binding
    integer(c_int) :: cslib_extract
    type(c_ptr), value :: ptr
    integer(c_int), value :: flag
  end function cslib_extract
end interface

end module cslib_wrap
