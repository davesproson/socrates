! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate transmission and reflection coefficients.
!
!- ---------------------------------------------------------------------
MODULE trans_source_coeff_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'TRANS_SOURCE_COEFF_MOD'
CONTAINS
SUBROUTINE trans_source_coeff(control, n_profile                        &
     , i_layer_first, i_layer_last                                      &
     , tau_dir, tau, sum, diff, lambda, lambda_r64, sec_0, path_div     &
     , gamma_up, gamma_down                                             &
     , trans, reflect, trans_0_dir, trans_0, source_coeff               &
     , nd_profile                                                       &
     , id_op_lt, id_op_lb, id_trs_lt, id_trs_lb                         &
     , nd_source_coeff                                                  &
     )


  USE, INTRINSIC :: iso_fortran_env, ONLY: real64
  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: ip_infra_red, ip_solar, ip_scf_solar_up,           &
                     ip_scf_solar_down, ip_direct_noscaling,            &
                     ip_direct_csr_scaling, ip_scf_ir_1d, ip_scf_ir_2d
  USE svectlib_mod, ONLY : exp_v
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE def_control, ONLY: StrCtrl

  IMPLICIT NONE

! Control options:
  TYPE(StrCtrl), INTENT(IN) :: control

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , id_op_lt                                                          &
!       Topmost declared layer for optical properties
    , id_op_lb                                                          &
!       Bottom declared layer for optical properties
    , id_trs_lt                                                         &
!       Topmost declared layer for transmission coefficients
    , id_trs_lb                                                         &
!       Bottom declared layer for transmission coefficients
    , nd_source_coeff
!       Size allocated for source coefficients


! Dummy variables.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , i_layer_first                                                     &
!       First layer to consider
    , i_layer_last
!       Last layer to consider

! Optical properties of the layer
  REAL (RealK), INTENT(IN) ::                                           &
      tau(nd_profile, id_op_lt: id_op_lb)                               &
!       Optical depths of layers
    , tau_dir(nd_profile, id_op_lt: id_op_lb)                           &
!       Unscaled optical depths of layers
    , sum(nd_profile, id_op_lt: id_op_lb)                               &
!       Sum of alpha_1 and alpha_2
    , diff(nd_profile, id_op_lt: id_op_lb)                              &
!       Difference of alpha_1 and alpha_2
    , lambda(nd_profile, id_op_lt: id_op_lb)                            &
!       Lambda
    , sec_0(nd_profile)                                                 &
!       Secant of solar zenith angle
    , path_div(nd_profile, id_trs_lt: id_trs_lb)
!       Path scaling for spherical geometry

! Set only in the solar case, and held in real64 whatever RealK is:
! see solar_coefficient_basic.
  REAL (real64), INTENT(IN) ::                                          &
      lambda_r64(nd_profile, id_op_lt: id_op_lb)                        &
!       Lambda
    , gamma_up(nd_profile, id_op_lt: id_op_lb)                          &
!       Basic solar coefficient for upward radiation
    , gamma_down(nd_profile, id_op_lt: id_op_lb)
!       Basic solar coefficient for downward radiation

! Transmission and reflection coefficients and coefficients for
! source terms.
  REAL (RealK), INTENT(OUT) ::                                          &
      trans(nd_profile, id_trs_lt: id_trs_lb)                           &
!       Diffuse transmission coefficient
    , reflect(nd_profile, id_trs_lt: id_trs_lb)                         &
!       Diffuse reflection coefficient
    , trans_0(nd_profile, id_trs_lt: id_trs_lb)                         &
!       Direct transmission coefficient
    , trans_0_dir(nd_profile, id_trs_lt: id_trs_lb)                     &
!       Direct transmission coefficient without scaling
    , source_coeff(nd_profile, id_trs_lt: id_trs_lb                     &
        , nd_source_coeff)
!       Source coefficients


! Local variables
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable
  REAL (RealK) ::                                                       &
      gamma                                                             &
!       Gamma
    , exponential                                                       &
!       Exponential of scaled optical depth
    , gamma2                                                            &
!       Gamma squared
    , exponential2
!       Exponential squared

! The solar source coefficients are small differences of large terms
! near lambda=sec_0, so all of these terms are calculated in real64.
! The transmission and reflection coefficients returned to the solver
! remain in RealK.
  REAL (real64) ::                                                      &
      trans_r64(nd_profile, id_trs_lt: id_trs_lb)                       &
!       Diffuse transmission coefficient
    , reflect_r64(nd_profile, id_trs_lt: id_trs_lb)                     &
!       Diffuse reflection coefficient
    , trans_0_r64                                                       &
!       Direct transmission coefficient
    , sum_r64                                                           &
!       Sum of alpha_1 and alpha_2
    , gamma_r64                                                         &
!       Gamma
    , gamma2_r64                                                        &
!       Gamma squared
    , exponential_r64                                                   &
!       Exponential of scaled optical depth
    , exponential2_r64                                                  &
!       Exponential squared
    , tmp_inv_r64
!       Temporary variable

! Variables related to the treatment of ill-conditioning
  REAL (RealK) ::                                                       &
      sq_eps_r                                                          &
!       The square root of a real number that is negligible
!       compared to 1.0.
    , tmp_inv                                                           &
!       Temporary work variable
    , tol
!       The tolerance used for switching to the asymptotic form
!       of the quadratic source function term.

  REAL (RealK) ::                                                       &
      temp(nd_profile),                                                 &
      temp_in(nd_profile)


  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='TRANS_SOURCE_COEFF'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the tolerances used in avoiding ill-conditioning, testing
! on any variable.
  sq_eps_r=SQRT(EPSILON(sq_eps_r))
  tol=SQRT(sq_eps_r)

! Determine the diffuse transmission and reflection coefficients.

  DO i=i_layer_first, i_layer_last
    DO l=1, n_profile
      temp_in(l)=-lambda(l,i)*tau(l,i)
    END DO
    CALL exp_v(n_profile, temp_in, temp)
    DO l=1, n_profile
      exponential=temp(l)
      exponential2=exponential*exponential
      gamma=(sum(l, i)-lambda(l, i)) / (sum(l, i)+lambda(l, i))
      gamma2=gamma*gamma
      tmp_inv=1.0_RealK / ( 1.0_RealK - exponential2*gamma2 )
      trans(l, i)=exponential*(1.0_RealK-gamma2)*tmp_inv
      reflect(l, i)=gamma*(1.0_RealK-exponential2)*tmp_inv
    END DO
  END DO



  IF (control%isolir == ip_solar) THEN

!   Calculate the direct transmission and the source coefficients
!   for the solar beam: in the solar case these are
!   the coefficients which will multiply the direct flux at the
!   top of the layer to give the source terms for the upward
!   diffuse flux and the total downward flux.

    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile
        sum_r64=REAL(sum(l, i), real64)
        exponential_r64=EXP(-lambda_r64(l, i)*REAL(tau(l, i), real64))
        exponential2_r64=exponential_r64*exponential_r64
        gamma_r64=(sum_r64-lambda_r64(l, i))/(sum_r64+lambda_r64(l, i))
        gamma2_r64=gamma_r64*gamma_r64
        tmp_inv_r64=1.0_real64/(1.0_real64-exponential2_r64*gamma2_r64)
        trans_r64(l, i)=exponential_r64*(1.0_real64-gamma2_r64)         &
          *tmp_inv_r64
        reflect_r64(l, i)=gamma_r64*(1.0_real64-exponential2_r64)       &
          *tmp_inv_r64
      END DO
    END DO

    IF (control%l_spherical_solar) THEN
      DO i=i_layer_first, i_layer_last
        DO l=1, n_profile
          temp(l) = -tau(l,i)*path_div(l,i)
        END DO
        CALL exp_v(n_profile,temp,trans_0(1,i))
        DO l=1, n_profile
          trans_0_r64=EXP(-REAL(tau(l, i), real64)                      &
            *REAL(path_div(l, i), real64))
          source_coeff(l, i, ip_scf_solar_up)=REAL(                     &
            gamma_up(l, i)-reflect_r64(l, i)*gamma_down(l, i)           &
            -gamma_up(l, i)*trans_r64(l, i)*trans_0_r64, RealK)
          source_coeff(l, i, ip_scf_solar_down)=REAL(trans_0_r64        &
            *(gamma_down(l, i)-gamma_up(l, i)*reflect_r64(l, i))        &
            -gamma_down(l, i)*trans_r64(l, i), RealK)
        END DO
      END DO
    ELSE
      DO i=i_layer_first, i_layer_last
        DO l=1, n_profile
          temp(l) = -tau(l,i)*sec_0(l)
        END DO
        CALL exp_v(n_profile,temp,trans_0(1,i))
        IF (control%i_direct_tau == ip_direct_noscaling .OR.            &
            control%i_direct_tau == ip_direct_csr_scaling) THEN
          DO l=1, n_profile
             temp(l) = -tau_dir(l,i)*sec_0(l)
          END DO
          CALL exp_v(n_profile,temp,trans_0_dir(1,i))
        END IF
        DO l=1, n_profile
          trans_0_r64=EXP(-REAL(tau(l, i), real64)                      &
            *REAL(sec_0(l), real64))
          source_coeff(l, i, ip_scf_solar_up)=REAL(                     &
            (gamma_up(l, i)-reflect_r64(l, i)                           &
            *(1.0_real64+gamma_down(l, i)))                             &
            -gamma_up(l, i)*trans_r64(l, i)*trans_0_r64, RealK)
          source_coeff(l, i, ip_scf_solar_down)=REAL(trans_0_r64        &
            *(1.0_real64+gamma_down(l, i)                               &
            -gamma_up(l, i)*reflect_r64(l, i))                          &
            -(1.0_real64+gamma_down(l, i))*trans_r64(l, i), RealK)
        END DO
      END DO
    END IF

  ELSE IF (control%isolir == ip_infra_red) THEN

!   In the case of infra-red radiation, the first source
!   coefficient holds the multiplier for the first difference
!   of the Planckian function across the layer, and the second
!   that for the second difference.

    DO i=i_layer_first, i_layer_last
      DO l=1, n_profile

!       A tolerance is added to the numerator and the denomiator
!       to avoid ill-conditioning at small optical depths.

        source_coeff(l, i, ip_scf_ir_1d)=(1.0_RealK-trans(l, i)         &
          +reflect(l, i)+sq_eps_r)                                      &
          /(sq_eps_r+tau(l, i)*sum(l, i))

      END DO
    END DO


    IF (control%l_ir_source_quad) THEN

!     Quadratic correction to source function.
!     This correction is very ill-conditioned for
!     small optical depths so the asymptotic form is then used.

      DO i=i_layer_first, i_layer_last
        DO l=1, n_profile
          IF (tau(l, i) > tol) THEN
            source_coeff(l, i, ip_scf_ir_2d)                            &
              =-2.0_RealK                                               &
              *(1.0_RealK-trans(l, i)-reflect(l, i)+sq_eps_r)           &
              /(diff(l, i)*tau(l, i)+sq_eps_r)
          ELSE
            source_coeff(l, i, ip_scf_ir_2d)                            &
              =-2.0_RealK+diff(l, i)*tau(l, i)
          END IF
          source_coeff(l, i, ip_scf_ir_2d)                              &
            =-(1.0_RealK+reflect(l, i)+trans(l, i)                      &
            +source_coeff(l, i, ip_scf_ir_2d))                          &
            /(sum(l, i)*tau(l, i)+sq_eps_r)
        END DO
      END DO

    END IF

  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE trans_source_coeff
END MODULE trans_source_coeff_mod
