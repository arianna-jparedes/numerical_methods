PROGRAM diffusion_bc

    USE netcdf
    IMPLICIT NONE

    REAL(8) :: dx, x0, x1, t0, t1, dt, tp, k, bc
    REAL(8), DIMENSION(:), ALLOCATABLE :: x, phi_new, phi_now
    INTEGER :: ncid, var_x, var_time, dim_x, dim_time, ierr, var_phi
    INTEGER :: n, t, nx, nt, itp, snap_idx, m, mp

    !Initial conditions
    dx = 0.01
    x0 = 0.0
    x1 = 1.0
    t0 = 0.0
    t1 = 21600.0
    tp = 3600.0
    k = 2.9e-5
    bc = 273.15
    dt = (dx**2)/(2.0*k)

    ! parameters
    nx = INT((x1 - x0) / dx) + 1
    nt = INT((t1 - t0) / dt) + 1
    itp = NINT(tp / dt)
    snap_idx = 0

    ALLOCATE(x(nx), phi_new(nx), phi_now(nx))

    ! first step
    DO n = 1, nx

        x(n) = x0 + (n-1)*dx
        phi_now(n) = phi0(x(n))

    END DO

    ! create nc
    ierr = nf90_create("diffusion.nc", NF90_CLOBBER, ncid)

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

    ! diffusion loop
    DO t = 2, nt

        phi_new = diffusion(phi_now, dx, dt, k, bc)
        phi_now = phi_new

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
        REAL(8):: pi

        pi = 4.0 * ATAN(1.0)

        IF (x >= 0 .AND. x <= 0.5) THEN
            phi = 273.15 + 20.0*x + SIN(50*pi*x)
        ELSE IF (x >= 0.5 .AND. x <= 1) THEN
            phi = 273.15 + 20.0 - 20.0*x + SIN(50*pi*x)
        END IF

    END FUNCTION phi0

    FUNCTION diffusion(phi_now, dx, dt, k, bc) RESULT(phi_new)

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: dx, dt, k, bc
        REAL(8), INTENT(IN) :: phi_now(:)
        REAL(8), ALLOCATABLE :: phi_new(:)
        REAL(8) :: m
        INTEGER :: n

        n = SIZE(phi_now)
        ALLOCATE(phi_new(n))

        m = (k*dt)/(dx*dx)

        ! boundary conditions
        phi_new(1) = bc
        phi_new(n) = bc

        phi_new(2:n-1) = phi_now(2:n-1) + m * (phi_now(3:n) - 2.0d0*phi_now(2:n-1) + phi_now(1:n-2))

    END FUNCTION diffusion


END PROGRAM diffusion_bc
