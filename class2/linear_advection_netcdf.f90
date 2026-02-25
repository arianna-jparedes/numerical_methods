PROGRAM linear_advection_netcdf
    USE netcdf
    IMPLICIT NONE

    REAL(8) :: u, c, dt, dx, x1, x0, t0, t1, tp
    REAL(8), DIMENSION(:), ALLOCATABLE :: x, phi, phi_new
    INTEGER :: n, t, nx, nt, itp, snap_idx
    INTEGER :: ncid, dim_x, dim_time, var_x, var_time, var_phi, ierr

    ! initial conditions
    u = 0.087
    dx = 0.1
    x0 = 0.0
    x1 = 100.0
    t0 = 0.0
    t1 = 1002.0
    tp = 200.0
    dt = 1.1
    c = u * dt / dx

    ! parameters
    nx = INT((x1 - x0) / dx) + 1
    nt = INT((t1 - t0) / dt) + 1
    itp = NINT(tp / dt)
    snap_idx = 0

    ALLOCATE(x(nx), phi(nx), phi_new(nx))

    ! first step
    DO n = 1, nx
        x(n) = x0 + (n-1)*dx
        phi(n) = phi0(x(n))
    END DO

    ! create nc
    ierr = nf90_create("advection.nc", NF90_CLOBBER, ncid)

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
    ierr = nf90_put_att(ncid, NF90_GLOBAL, "u", u )

    ierr = nf90_enddef(ncid)

    ! write nc
    ierr = nf90_put_var(ncid, var_x, x)

    snap_idx = snap_idx + 1
    call save_snapshot(ncid, var_time, var_phi, snap_idx, 0.0d0, phi)

    ! original loop
    DO t = 2, nt

        IF (u > 0) THEN
            phi_new = ftbs(phi, c)
        ELSE
            phi_new = ftfs(phi, c)
        END IF
        phi = phi_new

        ! save every tp seconds
        IF (MOD(t-1, itp) == 0) THEN
            snap_idx = snap_idx + 1
            call save_snapshot(ncid, var_time, var_phi, snap_idx, (snap_idx-1)*tp, phi)
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

    REAL(8) FUNCTION phi0(x) RESULT(phi_out)

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: x
        REAL(8) :: pi

        pi = 4.0 * ATAN(1.0)
        phi_out = 0.0

        IF (x >= 40.0 .AND. x <= 70.0) phi_out = SIN(pi*(x - 40.0)/30.0)**2

    END FUNCTION phi0

    FUNCTION ftfs(phi_now, c) RESULT(phi_out)

        IMPLICIT NONE
        REAL(8), DIMENSION(:), INTENT(IN) :: phi_now
        REAL(8), DIMENSION(SIZE(phi_now)) :: phi_out
        REAL(8), INTENT(IN) :: c
        INTEGER :: i, n

        n = SIZE(phi_now)

        DO i = 1, n-1
            phi_out(i) = (1.0 + c) * phi_now(i) - c * phi_now(i+1)
        END DO

        phi_out(n) = (1.0 + c) * phi_now(n) - c * phi_now(1)

    END FUNCTION ftfs

    FUNCTION ftbs(phi_now, c) RESULT(phi_out)

        IMPLICIT NONE
        REAL(8), DIMENSION(:), INTENT(IN) :: phi_now
        REAL(8), DIMENSION(SIZE(phi_now)) :: phi_out
        REAL(8), INTENT(IN) :: c
        INTEGER :: i, n

        n = SIZE(phi_now)
        phi_out(1) = (1.0 - c) * phi_now(1) + c * phi_now(n)

        DO i = 2, n
            phi_out(i) = (1.0 - c) * phi_now(i) + c * phi_now(i-1)
        END DO

    END FUNCTION ftbs

END PROGRAM linear_advection_netcdf
