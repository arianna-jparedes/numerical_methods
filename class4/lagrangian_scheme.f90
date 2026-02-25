PROGRAM lagrangian

    USE netcdf
    IMPLICIT NONE

    REAL(8) :: u, dx, x0, x1, t0, t1, tp, dt, c, alpha
    REAL(8), DIMENSION(:), ALLOCATABLE :: x, xdep, phi_lin, phi_lin_new, phi_cub, phi_cub_new
    INTEGER :: n, t, nx, nt, itp, snap_idx, m, mp
    INTEGER :: ncid, dim_x, dim_time, var_x, var_time, var_phi_lin, var_phi_cub, ierr

    !Initial conditions
    u = 0.75
    dx = 0.5
    x0 = 0.0
    x1 = 1000.0
    t0 = 0.0
    t1 = 2001.0
    tp = 250.0
    dt = 9.0
    c = u * dt / dx

    ! parameters
    nx = INT((x1 - x0) / dx) + 1
    nt = INT((t1 - t0) / dt) + 1
    itp = NINT(tp / dt)
    snap_idx = 0

    ALLOCATE(x(nx), xdep(nx), phi_lin(nx), phi_cub(nx), phi_lin_new(nx), phi_cub_new(nx))

    ! first step
    DO n = 1, nx

        x(n) = x0 + (n-1)*dx
        phi_lin(n) = phi0(x(n))
        phi_cub(n) = phi0(x(n))

    END DO

    ! create nc
    ierr = nf90_create("advection_lagrangian.nc", NF90_CLOBBER, ncid)

    ! dimensions
    ierr = nf90_def_dim(ncid, "x", nx, dim_x)
    ierr = nf90_def_dim(ncid, "time", NF90_UNLIMITED, dim_time)

    ! variables
    ierr = nf90_def_var(ncid, "x", NF90_DOUBLE, (/dim_x/), var_x)
    ierr = nf90_def_var(ncid, "time", NF90_DOUBLE, (/dim_time/), var_time)
    ierr = nf90_def_var(ncid, "phi_linear", NF90_DOUBLE, (/dim_x, dim_time/), var_phi_lin)
    ierr = nf90_def_var(ncid, "phi_cubic",  NF90_DOUBLE, (/dim_x, dim_time/), var_phi_cub)


    ! attributtes
    ierr = nf90_put_att(ncid, var_x, "units", "m")
    ierr = nf90_put_att(ncid, var_time, "units", "s")
    ierr = nf90_put_att(ncid, var_phi_lin, "units", "1")
    ierr = nf90_put_att(ncid, var_phi_cub, "units", "1")
    ierr = nf90_put_att(ncid, NF90_GLOBAL, "dx", dx)
    ierr = nf90_put_att(ncid, NF90_GLOBAL, "dt", dt)
    ierr = nf90_put_att(ncid, NF90_GLOBAL, "u", u )

    ierr = nf90_enddef(ncid)

    ! write nc
    ierr = nf90_put_var(ncid, var_x, x)

    snap_idx = snap_idx + 1
    call save_snapshot(ncid, var_time, var_phi_lin, snap_idx, 0.0d0, phi_lin)
    call save_snapshot(ncid, var_time, var_phi_cub, snap_idx, 0.0d0, phi_cub)


    ! Lagrangian schemes
    DO t = 2, nt

        DO n = 1, nx

            xdep(n) = x0 + MODULO((x(n) - u*dt - x0), (x1 - x0))
            phi_lin_new(n) = lagrangian_linear(phi_lin, xdep(n), dx, nx)
            phi_cub_new(n) = lagrangian_cubic(phi_cub, xdep(n), dx, nx)

        END DO

        phi_lin = phi_lin_new
        phi_cub = phi_cub_new

        ! save every tp seconds
        IF (MOD(t-1, itp) == 0) THEN

            snap_idx = snap_idx + 1
            call save_snapshot(ncid, var_time, var_phi_lin, snap_idx, (snap_idx-1)*tp, phi_lin)
            call save_snapshot(ncid, var_time, var_phi_cub, snap_idx, (snap_idx-1)*tp, phi_cub)

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
            phi = 0.1*(x - 400.0)
        ELSE IF (x >= 500.0 .AND. x <= 600.0) THEN
            phi = 20.0 - 0.1*(x - 400.0)
        ELSE
            phi = 0
        END IF

    END FUNCTION phi0

    REAL(8) FUNCTION lagrangian_linear(phi, xdep, dx, nx) RESULT(phi_out)

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: phi(:), xdep, dx
        INTEGER, INTENT(IN) :: nx
        INTEGER :: m, mp
        REAL(8) :: alpha

        m = INT(FLOOR(xdep/dx))
        alpha = (xdep/dx) - m
        m = m + 1
        mp = m + 1

        IF (mp > nx) mp = 2

        phi_out = (1.0d0 - alpha)*phi(m) + alpha*phi(mp)

    END FUNCTION lagrangian_linear

    REAL(8) FUNCTION lagrangian_cubic(phi, xdep, dx, nx) RESULT(phi_out)

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: phi(:), xdep, dx
        INTEGER, INTENT(IN) :: nx
        INTEGER :: m, im1, mp, mp2
        REAL(8) :: alpha, w1, w2, w3, w4

        m = INT(FLOOR(xdep/dx))
        alpha = (xdep/dx) - m
        m = m + 1

        im1 = m - 1
        mp  = m + 1
        mp2 = m + 2

        IF (im1 < 1 ) im1 = nx
        IF (mp  > nx) mp  = 2
        IF (mp2 > nx) mp2 = 3

        w1 = -alpha*(alpha-1.0)*(alpha-2.0)/6.0
        w2 =  (alpha+1.0)*(alpha-1.0)*(alpha-2.0)/2.0
        w3 = -(alpha+1.0)*alpha*(alpha-2.0)/2.0
        w4 =  (alpha+1.0)*alpha*(alpha-1.0)/6.0

        phi_out = w1*phi(im1) + w2*phi(m) + w3*phi(mp) + w4*phi(mp2)

    END FUNCTION lagrangian_cubic

END PROGRAM lagrangian

! Analysis
! We do not need to worry about the CFL condition for stability, because the solution at the new time
! is obtained by tracing the trajectory backward to its departure point and interpolating the value there
! rather than using an explicit update based on neighboring grid points, like the previous methods
! so, the method remains stable even when the CFL number is larger than one (I tested it!), but the CFL
! is important for solution quality, with a big CFL, we are going to have a bigger interpolation error.
