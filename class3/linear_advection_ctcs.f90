PROGRAM linear_advection
    USE netcdf
    IMPLICIT NONE
    REAL(8) :: u, c, dt, dx, x1, x0, t0, t1, tp, alpha_ra, alpha_raw, beta
    REAL(8), DIMENSION(:), ALLOCATABLE :: x, phi_new, phi_now, phi_old, d_ra, d_raw
    REAL(8), DIMENSION(:), ALLOCATABLE :: phi_now_ra, phi_now_raw, phi_old_ra, phi_old_raw, phi_new_ra, phi_new_raw
    INTEGER :: n, t, nx, nt, itp, snap_idx
    INTEGER :: ncid, var_x, var_time, dim_x, dim_time, ierr
    INTEGER :: var_phi_nonfilt, var_phi_rafilt, var_phi_rawfilt

    !Initial conditions
    u = -0.31
    dx = 0.1
    x0 = 0.0
    x1 = 500.0
    t0 = 0.0
    t1 = 1001.0
    tp = 200.0
    alpha_ra = 0.1
    alpha_raw = 0.05
    beta = 0.53

    dt = 0.1
    c = u * dt / dx
    nx = INT((x1 - x0) / dx) + 1
    nt = INT((t1 - t0) / dt) + 1
    itp = NINT(tp / dt)

    snap_idx = 0

    ALLOCATE(x(nx), phi_now(nx), phi_old(nx), phi_new(nx), d_ra(nx), d_raw(nx))
    ALLOCATE(phi_now_ra(nx), phi_now_raw(nx), phi_old_ra(nx), phi_old_raw(nx), phi_new_ra(nx), phi_new_raw(nx))

    DO n = 1, nx
        x(n) = x0 + (n-1)*dx
        phi_old(n) = phi0(x(n))
        phi_old_ra(n) = phi0(x(n))
        phi_old_raw(n) = phi0(x(n))
    END DO

    ! Initial conditions for arrays
    phi_now(:) = phi_old(:)
    phi_now_ra(:) = phi_old_ra(:)
    phi_now_raw(:) = phi_old_raw(:)
    d_ra(:) = phi_old_ra(:)
    d_raw(:) = phi_old_raw(:)


    IF (u > 0) THEN
        phi_now = ftbs(phi_old, c)
        phi_now_ra = ftbs(phi_old_ra, c)
        phi_now_ra = ftbs(phi_old_ra, c)
    ELSE
        phi_now = ftfs(phi_old, c)
        phi_now_ra = ftfs(phi_old_ra, c)
        phi_now_raw = ftfs(phi_old_raw, c)
    END IF

    ! create nc
    ierr = nf90_create("advection_ctcs.nc", NF90_CLOBBER, ncid)

    ! dimensions
    ierr = nf90_def_dim(ncid, "x", nx, dim_x)
    ierr = nf90_def_dim(ncid, "time", NF90_UNLIMITED, dim_time)

    ! variables
    ierr = nf90_def_var(ncid, "x", NF90_DOUBLE, (/dim_x/), var_x)
    ierr = nf90_def_var(ncid, "time", NF90_DOUBLE, (/dim_time/), var_time)
    ierr = nf90_def_var(ncid, "phi nonfilter", NF90_DOUBLE, (/dim_x, dim_time/), var_phi_nonfilt)
    ierr = nf90_def_var(ncid, "phi RA filter", NF90_DOUBLE, (/dim_x, dim_time/), var_phi_rafilt)
    ierr = nf90_def_var(ncid, "phi RAW filter", NF90_DOUBLE, (/dim_x, dim_time/), var_phi_rawfilt)

    ! attributtes
    ierr = nf90_put_att(ncid, var_x, "units", "m")
    ierr = nf90_put_att(ncid, var_time, "units", "s")
    ierr = nf90_put_att(ncid, var_phi_nonfilt, "units", "1")
    ierr = nf90_put_att(ncid, var_phi_rafilt, "units", "1")
    ierr = nf90_put_att(ncid, var_phi_rawfilt, "units", "1")
    ierr = nf90_put_att(ncid, NF90_GLOBAL, "dx", dx)
    ierr = nf90_put_att(ncid, NF90_GLOBAL, "dt", dt)
    ierr = nf90_put_att(ncid, NF90_GLOBAL, "u", u)

    ierr = nf90_enddef(ncid)

    ! write nc
    ierr = nf90_put_var(ncid, var_x, x)

    snap_idx = snap_idx + 1
    call save_snapshot(ncid, var_time, var_phi_nonfilt, snap_idx, 0.0d0, phi_now)
    call save_snapshot(ncid, var_time, var_phi_rafilt, snap_idx, 0.0d0, phi_now_ra)
    call save_snapshot(ncid, var_time, var_phi_rawfilt, snap_idx, 0.0d0, phi_now_raw)

    DO t = 2, nt
        phi_new = ctcs(phi_old, phi_now, c)
        phi_new_ra = ctcs(phi_old_ra, phi_now_ra, c)
        phi_new_raw = ctcs(phi_old_raw, phi_now_raw, c)

        ! Non filter
        phi_old(:) = phi_now(:)
        phi_now(:) = phi_new(:)

        ! RA filter
        d_ra = alpha_ra*(phi_old_ra + phi_new_ra - 2.0*phi_now_ra)
        phi_old_ra(:) = phi_now_ra + d_ra
        phi_now_ra(:) = phi_new_ra

        ! RAW filter
        d_raw = alpha_raw*(phi_old_raw + phi_new_raw - 2.0*phi_now_raw)
        phi_old_raw(:) = phi_now_raw + beta*d_raw
        phi_now_raw(:) = phi_new_raw + (1 - beta)*d_raw

        ! save every tp seconds
        IF (MOD(t-1, itp) == 0) THEN
            snap_idx = snap_idx + 1
            call save_snapshot(ncid, var_time, var_phi_nonfilt, snap_idx, (snap_idx-1)*tp, phi_now)
            call save_snapshot(ncid, var_time, var_phi_rafilt, snap_idx, (snap_idx-1)*tp, phi_now_ra)
            call save_snapshot(ncid, var_time, var_phi_rawfilt, snap_idx, (snap_idx-1)*tp, phi_now_raw)
        END IF

    END DO

    ierr = nf90_close(ncid)

CONTAINS

    SUBROUTINE save_snapshot(ncid, var_time, var_phi, k, tsec, phi_vec)
        INTEGER, INTENT(IN) :: ncid, var_time, var_phi, k
        REAL(8), INTENT(IN) :: tsec
        REAL(8), DIMENSION(:), INTENT(IN) :: phi_vec
        INTEGER :: ierr, start_time(1), count_time(1), start_phi(2), count_phi(2)

        ! write time
        start_time = (/k/)
        count_time = (/1/)
        ierr = nf90_put_var(ncid, var_time, (/tsec/), start=start_time, count=count_time)

        ! write phi
        start_phi = (/1, k/)
        count_phi = (/SIZE(phi_vec), 1/)
        ierr = nf90_put_var(ncid, var_phi, phi_vec, start=start_phi, count=count_phi)

    END SUBROUTINE save_snapshot

    REAL(8) FUNCTION phi0(x) RESULT(phi)

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: x

        IF (x >= 200 .AND. x <= 250) THEN
            phi = 2.0
        ELSE IF (x >= 250 .AND. x <= 300) THEN
            phi = 1.0
        ELSE
            phi = 0.1
        END IF

    END FUNCTION phi0

    FUNCTION ftfs(phi_now, c) RESULT(phi_out)

        IMPLICIT NONE
        REAL(8), DIMENSION(:), INTENT(IN) :: phi_now
        REAL(8), DIMENSION(SIZE(phi_now)) :: phi_out
        REAL(8) :: c
        INTEGER :: i, n

        n = SIZE(phi_now)

        DO i = 1, n-1
            phi_out(i) = (1 + c) * phi_now(i) - c * phi_now(i+1)
        END DO

        phi_out(n) = (1 + c) * phi_now(n) - c * phi_now(1)

    END FUNCTION ftfs

    FUNCTION ftbs(phi_now, c) RESULT(phi_out)

        IMPLICIT NONE
        REAL(8), DIMENSION(:), INTENT(IN) :: phi_now
        REAL(8), DIMENSION(SIZE(phi_now)) :: phi_out
        REAL(8) :: c
        INTEGER :: i, n

        n = SIZE(phi_now)
        phi_out(1) = (1 - c) * phi_now(1) + c * phi_now(n)

        DO i = 2, n
            phi_out(i) = (1 - c) * phi_now(i) + c * phi_now(i-1)
        END DO

    END FUNCTION ftbs


    FUNCTION ctcs(phi_old, phi_now, c) RESULT(phi_out)

        IMPLICIT NONE
        REAL(8), DIMENSION(:), INTENT(IN) :: phi_now, phi_old
        REAL(8), DIMENSION(SIZE(phi_now)) :: phi_out
        REAL(8) :: c
        INTEGER :: i, n

        n = SIZE(phi_now)
        phi_out(1) = phi_old(1) - c*(phi_now(2) - phi_now(n))

        DO i = 2, n-1
            phi_out(i) = phi_old(i) - c*(phi_now(i+1) - phi_now(i-1))
        END DO

        phi_out(n) = phi_old(n) - c*(phi_now(1) - phi_now(n-1))

    END FUNCTION ctcs

END PROGRAM linear_advection
