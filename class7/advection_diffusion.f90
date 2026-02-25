PROGRAM advection_diffusion

    USE netcdf
    IMPLICIT NONE

    REAL(8) :: u, c, dx, dt, x1, x0, t0, t1, tp, alpha, beta, k, dc, dt_adv, dt_diff
    REAL(8), DIMENSION(:), ALLOCATABLE :: x, phi_new, phi_now, phi_old, d
    INTEGER :: n, t, nx, nt, itp, snap_idx, i
    INTEGER :: ncid, var_phi, var_x, var_time, dim_x, dim_time, ierr

    !Initial conditions
    u = 0.95
    dx = 0.2
    x0 = 0.0
    x1 = 1000.0
    t0 = 0.0
    t1 = 2001.0
    tp = 500.0
    k = 0.029
    alpha = 0.1
    beta = 0.53

    ! CFL
    dt_adv = dx / u
    dt_diff = (dx**2) / (2.0 * k)
    dt = 0.4 * MIN(dt_adv, dt_diff)

    c = u * dt / dx
    dc = k * dt / (dx**2)

    ! parameters
    nx = INT((x1 - x0) / dx) + 1
    nt = INT((t1 - t0) / dt) + 1
    itp = NINT(tp / dt)
    snap_idx = 0

    ALLOCATE(x(nx), phi_new(nx), phi_now(nx), phi_old(nx), d(nx))

    ! first step
    DO n = 1, nx

        x(n) = x0 + (n-1)*dx
        phi_old(n) = phi0(x(n))
        phi_now(n) = phi0(x(n))

    END DO

    ! create nc
    ierr = nf90_create("advection_diffusion.nc", NF90_CLOBBER, ncid)

    ! dimensions
    ierr = nf90_def_dim(ncid, "x", nx, dim_x)
    ierr = nf90_def_dim(ncid, "time", NF90_UNLIMITED, dim_time)

    ! variables
    ierr = nf90_def_var(ncid, "x", NF90_DOUBLE, (/dim_x/), var_x)
    ierr = nf90_def_var(ncid, "time", NF90_DOUBLE, (/dim_time/), var_time)
    ierr = nf90_def_var(ncid, "phi", NF90_DOUBLE, (/dim_x, dim_time/), var_phi)


    ! attributtes
    ierr = nf90_put_att(ncid, var_x, "units", "m")
    ierr = nf90_put_att(ncid, var_time, "units", "s")
    ierr = nf90_put_att(ncid, var_phi, "units", "1")
    ierr = nf90_put_att(ncid, NF90_GLOBAL, "dx", dx)
    ierr = nf90_put_att(ncid, NF90_GLOBAL, "dt", dt)

    ierr = nf90_enddef(ncid)

    ! write nc
    ierr = nf90_put_var(ncid, var_x, x)

    snap_idx = snap_idx + 1
    call save_snapshot(ncid, var_time, var_phi, snap_idx, 0.0d0, phi_now)

    ! ctcs scheme with diffusion
    phi_now = ftcs_diff(phi_now, c, dc)

    DO t = 2, nt

        phi_new = ctcs_diff(phi_old, phi_now, c, dc)
        d = alpha*(phi_old + phi_new - 2.0*phi_now)

        phi_old = phi_now + beta*d
        phi_now = phi_new + (1-beta)*d

        ! save every tp seconds
        IF (MOD(t-1, itp) == 0) THEN
            snap_idx = snap_idx + 1
            call save_snapshot(ncid, var_time, var_phi, snap_idx, (snap_idx-1)*tp, phi_now)
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

        IF (x >= 400.0 .AND. x <= 500.0) THEN
            phi = 0.01*(x - 400.0)
        ELSE IF (x >= 500.0 .AND. x <= 600) THEN
            phi = 2.0 - 0.01*(x - 400.0)
        ELSE
            phi = 0.0
        END IF

    END FUNCTION phi0

    FUNCTION ftcs_diff(phi_now, c, dc) RESULT(phi_out)

        IMPLICIT NONE
        REAL(8), DIMENSION(:), INTENT(IN) :: phi_now
        REAL(8), DIMENSION(SIZE(phi_now)) :: phi_out
        REAL(8), INTENT(IN) :: c, dc
        INTEGER :: i, n

        n = SIZE(phi_now)

        phi_out(1) = phi_now(1) - 0.5*c*(phi_now(2) - phi_now(n)) + dc*(phi_now(2) - 2.0*phi_now(1) + phi_now(n))

        DO i = 2, n-1
            phi_out(i) = phi_now(i) - 0.5*c*(phi_now(i+1) - phi_now(i-1)) + dc*(phi_now(i+1) - 2.0*phi_now(i) + phi_now(i-1))
        END DO

        phi_out(n) = phi_now(n) - 0.5*c*(phi_now(1) - phi_now(n-1)) + dc*(phi_now(1) - 2.0*phi_now(n) + phi_now(n-1))

    END FUNCTION ftcs_diff

    FUNCTION ctcs_diff(phi_old, phi_now, c, dc) RESULT(phi_out)

        IMPLICIT NONE
        REAL(8), DIMENSION(:), INTENT(IN) :: phi_now, phi_old
        REAL(8), DIMENSION(SIZE(phi_now)) :: phi_out
        REAL(8), INTENT(IN) :: c, dc
        INTEGER :: i, n

        n = SIZE(phi_now)

        phi_out(1) = phi_old(1) - c*(phi_now(2) - phi_now(n)) + 2.0*dc*(phi_old(2) - 2.0*phi_old(1) + phi_old(n))

        DO i = 2, n-1
            phi_out(i) = phi_old(i) - c*(phi_now(i+1) - phi_now(i-1)) + 2.0*dc*(phi_old(i+1) - 2.0*phi_old(i) + phi_old(i-1))
        END DO

        phi_out(n) = phi_old(n) - c*(phi_now(1) - phi_now(n-1)) + 2.0*dc*(phi_old(1) - 2.0*phi_old(n) + phi_old(n-1))

    END FUNCTION ctcs_diff


END PROGRAM advection_diffusion
