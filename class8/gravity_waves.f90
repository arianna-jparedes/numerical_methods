PROGRAM gravity_wave

    USE netcdf
    IMPLICIT NONE

    REAL(8) :: u0, x0, x1, dx, t0, t1, tp, dt, dt_c, c, dtdx, p
    REAL(8), DIMENSION(:), ALLOCATABLE :: x, p_now, u_now, p_old, u_old, p_new, u_new
    REAL(8), DIMENSION(:), ALLOCATABLE :: u_new2, p_new2, u_now2, p_now2
    INTEGER :: n, t, nx, nt, itp, snap_idx, i
    INTEGER :: ncid, var_p, var_u, var_x, var_time, dim_x, dim_time, ierr, var2_p, var2_u

    !Initial conditions
    u0 = 0.0
    dx = 0.5
    x0 = 0.0
    x1 = 1000.0
    t0 = 0.0
    t1 = 2001.0
    tp = 200.0
    p = 1.0

    ! CFL
    dt = 0.1
    dtdx = dt/dx
    c = SQRT(p) * dtdx


    ! parameters
    nx = INT((x1 - x0) / dx) + 1
    nt = INT((t1 - t0) / dt) + 1
    itp = NINT(tp / dt)
    snap_idx = 0

    ALLOCATE(x(nx), p_now(nx), u_now(nx), p_old(nx), u_old(nx), p_new(nx), u_new(nx))

    ! first step
    DO n = 1, nx
        x(n) = x0 + (n-1)*dx
        p_old(n) = phi0(x(n))
        u_old(n) = u0
    END DO

    ! create nc
    ierr = nf90_create("gravity_wave.nc", NF90_CLOBBER, ncid)

    ! dimensions
    ierr = nf90_def_dim(ncid, "x", nx, dim_x)
    ierr = nf90_def_dim(ncid, "time", NF90_UNLIMITED, dim_time)

    ! variables
    ierr = nf90_def_var(ncid, "x", NF90_DOUBLE, (/dim_x/), var_x)
    ierr = nf90_def_var(ncid, "time", NF90_DOUBLE, (/dim_time/), var_time)
    ierr = nf90_def_var(ncid, "p", NF90_DOUBLE, (/dim_x, dim_time/), var_p)
    ierr = nf90_def_var(ncid, "u", NF90_DOUBLE, (/dim_x, dim_time/), var_u)
    ierr = nf90_def_var(ncid, "p2", NF90_DOUBLE, (/dim_x, dim_time/), var2_p)
    ierr = nf90_def_var(ncid, "u2", NF90_DOUBLE, (/dim_x, dim_time/), var2_u)

    ! attributtes
    ierr = nf90_put_att(ncid, var_x, "units", "m")
    ierr = nf90_put_att(ncid, var_time, "units", "s")
    ierr = nf90_put_att(ncid, var_p, "units", "m2 s-2")
    ierr = nf90_put_att(ncid, var_u, "units", "m2 s-1")
    ierr = nf90_put_att(ncid, var2_p, "units", "m2 s-2")
    ierr = nf90_put_att(ncid, var2_u, "units", "m2 s-1")
    ierr = nf90_put_att(ncid, NF90_GLOBAL, "dx", dx)
    ierr = nf90_put_att(ncid, NF90_GLOBAL, "dt", dt)

    ierr = nf90_enddef(ncid)

    ! write nc
    ierr = nf90_put_var(ncid, var_x, x)

    snap_idx = snap_idx + 1
    call save_snapshot(ncid, var_time, var_p, snap_idx, 0.0d0, p_old)
    call save_snapshot(ncid, var_time, var_u, snap_idx, 0.0d0, u_old)
    call save_snapshot(ncid, var_time, var2_p, snap_idx, 0.0d0, p_old)
    call save_snapshot(ncid, var_time, var2_u, snap_idx, 0.0d0, u_old)

    p_now = peq_ft(p_old, u_old, c)
    u_now = ueq_ft(u_old, p_old, dtdx)

    p_now2 = p_old
    u_now2 = u_old

    DO t = 2, nt

        ! First method
        u_new = ueq(u_old, p_now, dtdx)
        p_new = peq(p_old, u_now, c)

        p_old = p_now
        p_now = p_new
        u_old = u_now
        u_now = u_new

        ! Second method
        u_new2 = ueq_fo(u_now2, p_now2, dtdx)
        p_new2 = peq_fo(p_now2, u_new2, c)

        p_now2 = p_new2
        u_now2 = u_new2

        ! save every tp seconds
        IF (MOD(t-1, itp) == 0) THEN
            snap_idx = snap_idx + 1
            call save_snapshot(ncid, var_time, var_p, snap_idx, (snap_idx-1)*tp, p_now)
            call save_snapshot(ncid, var_time, var_u, snap_idx, (snap_idx-1)*tp, u_now)
            call save_snapshot(ncid, var_time, var2_p, snap_idx, (snap_idx-1)*tp, p_now2)
            call save_snapshot(ncid, var_time, var2_u, snap_idx, (snap_idx-1)*tp, u_now2)
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
        REAL(8) :: pi

        pi = 4.0*ATAN(1.0)

        IF (x < 400.0) THEN
            phi = 0.0
        ELSE IF (x >= 400.0 .AND. x < 600) THEN
            phi = SIN(((X - 400.0)/200)*pi)**2
        ELSE
            phi = 0.0
        END IF

    END FUNCTION phi0


    FUNCTION ueq_ft(u_now, p_now, dtdx) RESULT(u_out)

        IMPLICIT NONE
        REAL(8), DIMENSION(:), INTENT(IN) :: u_now, p_now
        REAL(8), DIMENSION(SIZE(u_now)) :: u_out
        REAL(8), INTENT(IN) :: dtdx
        INTEGER :: n, i

        n = SIZE(u_now)

        u_out(1) = u_now(1) - 0.5*dtdx*(p_now(2) - p_now(n))

        DO i = 2, n - 1

            u_out(i) = u_now(i) - 0.5*dtdx*(p_now(i+1) - p_now(i-1))

        END DO

        u_out(n) = u_now(n) - 0.5*dtdx*(p_now(1) - p_now(n-1))


    END FUNCTION ueq_ft

    FUNCTION peq_ft(p_now, u_now, c) RESULT(p_out)

        IMPLICIT NONE
        REAL(8), DIMENSION(:), INTENT(IN) :: u_now, p_now
        REAL(8), DIMENSION(SIZE(p_now)) :: p_out
        REAL(8), INTENT(IN) :: c
        INTEGER :: n, i

        n = SIZE(p_now)

        p_out(1) = p_now(1) - 0.5*c*(u_now(2) - u_now(n))

        DO i = 2, n - 1

            p_out(i) = p_now(i) - 0.5*c*(u_now(i+1) - u_now(i-1))

        END DO

        p_out(n) = p_now(n) - 0.5*c*(u_now(1) - u_now(n-1))


    END FUNCTION peq_ft

        FUNCTION ueq(u_old, p_now, dtdx) RESULT(u_new)

        IMPLICIT NONE
        REAL(8), DIMENSION(:), INTENT(IN) :: u_old, p_now
        REAL(8), DIMENSION(SIZE(u_old)) :: u_new
        REAL(8), INTENT(IN) :: dtdx
        INTEGER :: n, i

        n = SIZE(u_old)

        u_new(1) = u_old(1) - dtdx*(p_now(2) - p_now(n))

        DO i = 2, n - 1

            u_new(i) = u_old(i) - dtdx*(p_now(i+1) - p_now(i-1))

        END DO

        u_new(n) = u_old(n) - dtdx*(p_now(1) - p_now(n-1))


    END FUNCTION ueq

    FUNCTION peq(p_old, u_now, c) RESULT(p_new)

        IMPLICIT NONE
        REAL(8), DIMENSION(:), INTENT(IN) :: u_now, p_old
        REAL(8), DIMENSION(SIZE(p_old)) :: p_new
        REAL(8), INTENT(IN) :: c
        INTEGER :: n, i

        n = SIZE(p_old)

        p_new(1) = p_old(1) - c*(u_now(2) - u_now(n))

        DO i = 2, n - 1

            p_new(i) = p_old(i) - c*(u_now(i+1) - u_now(i-1))

        END DO

        p_new(n) = p_old(n) - c*(u_now(1) - u_now(n-1))


    END FUNCTION peq

    FUNCTION ueq_fo(u_now, p_now, dtdx) RESULT(u_new)

        IMPLICIT NONE
        REAL(8), DIMENSION(:), INTENT(IN) :: u_now, p_now
        REAL(8), DIMENSION(SIZE(u_now)) :: u_new
        REAL(8), INTENT(IN) :: dtdx
        INTEGER :: n, i

        n = SIZE(u_now)

        u_new(1) = u_now(1) - 0.5*dtdx*(p_now(2) - p_now(n))

        DO i=2, n-1

            u_new(i) = u_now(i) - 0.5*dtdx*(p_now(i+1) - p_now(i-1))

        END DO

        u_new(n) = u_now(n) - 0.5*dtdx*(p_now(1) - p_now(n-1))

    END FUNCTION ueq_fo

    FUNCTION peq_fo(p_now, u_fut, c) RESULT(p_new)

        IMPLICIT NONE
        REAL(8), DIMENSION(:), INTENT(IN) :: p_now, u_fut
        REAL(8), DIMENSION(SIZE(p_now)) :: p_new
        REAL(8), INTENT(IN) :: c
        INTEGER :: n, i

        n = SIZE(p_now)

        p_new(1) = p_now(1) - 0.5*c*(u_fut(2) - u_fut(n))

        DO i = 2, n - 1

            p_new(i) = p_now(i) - 0.5*c*(u_fut(i+1) - u_fut(i-1))

        END DO

        p_new(n) = p_now(n) - 0.5*c*(u_fut(1) - u_fut(n-1))


    END FUNCTION peq_fo



END PROGRAM gravity_wave

! Both solutions are the same, we should have select different dt for each method, but as we choose a
! small enough dt, it works for both.
