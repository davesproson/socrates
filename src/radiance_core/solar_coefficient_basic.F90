! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate the basic coefficients for the solar beam.
!
!- ---------------------------------------------------------------------
MODULE solar_coefficient_basic_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SOLAR_COEFFICIENT_BASIC_MOD'
CONTAINS
SUBROUTINE solar_coefficient_basic(control                              &
    , n_profile, i_layer_first, i_layer_last                            &
    , omega, asymmetry, sec_0, path_div                                 &
    , sum, diff, lambda, lambda_r64                                     &
    , gamma_up, gamma_down                                              &
    , nd_profile, id_lt, id_lb, id_trs_lt, id_trs_lb                    &
    )


  USE, INTRINSIC :: iso_fortran_env, ONLY: real64
  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: i_err_fatal, ip_2s_test, ip_discrete_ord,          &
                     ip_eddington, ip_elsasser, ip_hemi_mean, ip_pifm80,&
                     ip_pifm85
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE def_control, ONLY: StrCtrl
  USE ereport_mod, ONLY: ereport
  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl), INTENT(IN) :: control

! Sizes of dummy arrays.
  INTEGER                                                               &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , id_lt                                                             &
!       Topmost declared layer
    , id_lb                                                             &
!       Bottom declared layer
    , id_trs_lt                                                         &
!       Topmost declared layer for transmission coefficients
    , id_trs_lb
!       Bottom declared layer for transmission coefficients


! Dummy variables.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , i_layer_first                                                     &
!       First layer to consider
    , i_layer_last
!       First layer to consider

  REAL (RealK), INTENT(IN) ::                                           &
      omega(nd_profile, id_lt: id_lb)                                   &
!       Albedo of single scattering
    , asymmetry(nd_profile, id_lt: id_lb)                               &
!       Asymmetry
    , sec_0(nd_profile)                                                 &
!       Secant of solar zenith angle
    , path_div(nd_profile, id_trs_lt: id_trs_lb)
!       Scaled path for spherical geometry

! Basic two-stream coefficients:
  REAL (RealK), INTENT(INOUT) ::                                        &
      sum(nd_profile, id_lt: id_lb)                                     &
!       Sum of two-stream coefficients
    , diff(nd_profile, id_lt: id_lb)                                    &
!       Difference of two-stream coefficients
    , lambda(nd_profile, id_lt: id_lb)
!       Lambda
! The coefficients below are held in real64 whatever RealK is. They
! scale as 1/(lambda-sec_0), and trans_source_coeff subtracts them
! from each other, so near lambda=sec_0 rounding in RealK=real32
! would swamp the source terms.
  REAL (real64), INTENT(OUT) ::                                         &
      lambda_r64(nd_profile, id_lt: id_lb)                              &
!       Lambda, consistent with any perturbation of sum and diff
    , gamma_up(nd_profile, id_lt: id_lb)                                &
!       Coefficient for upward radiation
    , gamma_down(nd_profile, id_lt: id_lb)
!       Coefficient for downwad radiation


! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable
  REAL (real64) ::                                                      &
      ksi_0(nd_profile, id_lt: id_lb)                                   &
!       Difference in solar scattering fractions
    , factor                                                            &
!       Temporary variable
    , sum_r64                                                           &
!       Sum of two-stream coefficients
    , diff_r64                                                          &
!       Difference of two-stream coefficients
    , sec_r64
!       Secant of the solar zenith angle or scaled path
  REAL (real64) ::                                                      &
      root_3
!       Square root of 3
  PARAMETER(                                                            &
      root_3=1.7320508075688772e+00_real64                              &
    )

! Variables related to the treatment of ill-conditioning
  REAL (RealK) ::                                                       &
      tol_perturb
!       The tolerance used to judge where the two-stream
!       expressions for the solar source become ill-conditioned

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  INTEGER :: ierr
  CHARACTER (LEN=errormessagelength) :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'SOLAR_COEFFICIENT_BASIC'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the tolerances used in avoiding ill-conditioning, testing
! on any variable.
  tol_perturb=3.2e+01_RealK*EPSILON(sec_0(1))

! If LAMBDA is too close to SEC_0 it must be perturbed.
  IF (control%l_spherical_solar) THEN
    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile
        IF (ABS(lambda(l, i)-path_div(l, i)) < tol_perturb) THEN
          sum(l, i)=(1.0e+00_RealK+tol_perturb)*sum(l, i)
          diff(l, i)=(1.0e+00_RealK+tol_perturb)*diff(l, i)
          lambda(l, i)=(1.0e+00_RealK+tol_perturb)*lambda(l, i)
        END IF
      END DO
    END DO
  ELSE
    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile
        IF ((ABS(lambda(l, i)-sec_0(l))) < tol_perturb) THEN
          sum(l, i)=(1.0e+00_RealK+tol_perturb)*sum(l, i)
          diff(l, i)=(1.0e+00_RealK+tol_perturb)*diff(l, i)
          lambda(l, i)=(1.0e+00_RealK+tol_perturb)*lambda(l, i)
        END IF
      END DO
    END DO
  END IF

  IF ( (control%i_2stream == ip_eddington).OR.                          &
       (control%i_2stream == ip_elsasser).OR.                           &
       (control%i_2stream == ip_pifm85).OR.                             &
       (control%i_2stream == ip_2s_test).OR.                            &
       (control%i_2stream == ip_hemi_mean).OR.                          &
       (control%i_2stream == ip_pifm80) ) THEN

    IF (control%l_spherical_solar) THEN
      DO i=i_layer_first, i_layer_last
        DO l=1, n_profile
          ksi_0(l, i)=1.5e+00_real64*asymmetry(l, i)/path_div(l, i)
        END DO
      END DO
    ELSE
      DO i=i_layer_first, i_layer_last
        DO l=1, n_profile
          ksi_0(l, i)=1.5e+00_real64*asymmetry(l, i)/sec_0(l)
        END DO
      END DO
    END IF

  ELSE IF (control%i_2stream == ip_discrete_ord) THEN

    IF (control%l_spherical_solar) THEN
      DO i=i_layer_first, i_layer_last
        DO l=1, n_profile
          ksi_0(l, i)=root_3*asymmetry(l, i)/path_div(l, i)
        END DO
      END DO
    ELSE
      DO i=i_layer_first, i_layer_last
        DO l=1, n_profile
          ksi_0(l, i)=root_3*asymmetry(l, i)/sec_0(l)
        END DO
      END DO
    END IF

  ELSE

    cmessage = '*** Error: An illegal solar two-stream scheme has '     &
      //'been selected.'
    ierr=i_err_fatal
    CALL ereport(RoutineName, ierr, cmessage)

  END IF


! Determine the basic solar coefficients for the two-stream equations.
  IF (control%l_spherical_solar) THEN
    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile
        sum_r64=REAL(sum(l, i), real64)
        diff_r64=REAL(diff(l, i), real64)
        sec_r64=REAL(path_div(l, i), real64)
        lambda_r64(l, i)=SQRT(sum_r64*diff_r64)
        factor=0.5e+00_real64*omega(l, i)*sec_r64                       &
          /((lambda_r64(l, i)-sec_r64)*(lambda_r64(l, i)+sec_r64))
        gamma_up(l, i)=factor*(sum_r64-sec_r64                          &
          -ksi_0(l, i)*(diff_r64-sec_r64))
        gamma_down(l, i)=factor*(sum_r64+sec_r64                        &
          +ksi_0(l, i)*(diff_r64+sec_r64))
      END DO
    END DO
  ELSE
    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile
        sum_r64=REAL(sum(l, i), real64)
        diff_r64=REAL(diff(l, i), real64)
        sec_r64=REAL(sec_0(l), real64)
        lambda_r64(l, i)=SQRT(sum_r64*diff_r64)
        factor=0.5e+00_real64*omega(l, i)*sec_r64                       &
          /((lambda_r64(l, i)-sec_r64)*(lambda_r64(l, i)+sec_r64))
        gamma_up(l, i)=factor*(sum_r64-sec_r64                          &
          -ksi_0(l, i)*(diff_r64-sec_r64))
        gamma_down(l, i)=factor*(sum_r64+sec_r64                        &
          +ksi_0(l, i)*(diff_r64+sec_r64))
      END DO
    END DO
  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE solar_coefficient_basic
END MODULE solar_coefficient_basic_mod
