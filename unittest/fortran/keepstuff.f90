MODULE keepstuff
  USE liblammps
  IMPLICIT NONE
  TYPE(LAMMPS) :: lmp
  INTEGER :: mycomm
  CHARACTER(len=40), DIMENSION(3), PARAMETER :: demo_input = &
      [ CHARACTER(len=40) ::                                &
      'region       box block 0 $x 0 2 0 2',                &
      'create_box 1 box',                                   &
      'create_atoms 1 single 1.0 1.0 ${zpos}' ]
  CHARACTER(LEN=40), DIMENSION(3), PARAMETER :: big_input = &
      [ CHARACTER(len=40) ::                                &
      'region       box block 0 $x 0 3 0 4',                &
      'create_box 1 box',                                   &
      'create_atoms 1 single 1.0 1.0 ${zpos}' ]
  CHARACTER(len=40), DIMENSION(2), PARAMETER :: cont_input = &
      [ CHARACTER(len=40) ::                                &
      'create_atoms 1 single &',                            &
      ' 0.2 0.1 0.1' ]
  CHARACTER(LEN=40), DIMENSION(3), PARAMETER :: pair_input = &
      [ CHARACTER(LEN=40) ::                                &
      'pair_style lj/cut 2.5',                              &
      'pair_coeff 1 1 1.0 1.0',                             &
      'mass 1 1.0' ]
END MODULE keepstuff

