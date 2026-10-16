! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to set the precision of real variables.

MODULE realtype_rd

  USE, INTRINSIC :: iso_fortran_env, ONLY: real32, real64

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: RealK, RealExt

  ! Internal Real precision within SOCRATES
  ! (set with -DSOCRATES_PRECISION=32, default 64)
#if (SOCRATES_PRECISION == 32)
  INTEGER, PARAMETER :: RealK=real32
#else
  INTEGER, PARAMETER :: RealK=real64
#endif

  ! External Real precision for variables passed through the Runes interface
  INTEGER, PARAMETER :: RealExt=real64

END MODULE realtype_rd
