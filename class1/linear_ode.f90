PROGRAM linear_ode
    ! Here we solve a linear ode on a discrete grid, and compare the accuracy of
    ! two methods: Euler forward(1st order) and Heun(2nd order).

    USE netcdf
    IMPLICIT NONE

    REAL(8) :: x0, x1, dx, y0
    INTEGER :: nx, n, stats, lode, dim_x, var_euler, var_heun, var_x
    REAL(8), ALLOCATABLE :: x(:), err_eul(:), err_heu(:)
    REAL(8) :: yeul, yheu, yexc

    ! Parameters
    x0 = 0.0
    x1 = 10.0
    dx = 0.1
    y0 = 0.0
    nx = INT((x1 - x0)/dx) + 1

    ALLOCATE(x(nx), err_eul(nx), err_heu(nx))

    ! Build x grid
    DO n = 1, nx
        x(n) = x0 + DBLE(n-1)*dx
    END DO

    ! Initial conditions
    yeul = y0
    yheu = y0
    yexc = sol_func(x(1))

    ! Errors
    err_eul(1)  = yeul - yexc
    err_heu(1)  = yheu - yexc

    ! Create nc
    stats = nf90_create("euler_heun.nc", NF90_CLOBBER, lode); CALL handle_err(stats)
    stats = nf90_def_dim(lode, "x", nx, dim_x); CALL handle_err(stats)

    ! variables
    stats = nf90_def_var(lode, "x", NF90_DOUBLE, (/dim_x/), var_x); CALL handle_err(stats)
    stats = nf90_def_var(lode, "euler_error", NF90_DOUBLE, (/dim_x/), var_euler); CALL handle_err(stats)
    stats = nf90_def_var(lode, "heun_error", NF90_DOUBLE, (/dim_x/), var_heun); CALL handle_err(stats)

    ! attributtes
    stats = nf90_put_att(lode, var_x, "units", "m"); CALL handle_err(stats)
    stats = nf90_put_att(lode, var_euler, "units", "m"); CALL handle_err(stats)
    stats = nf90_put_att(lode, var_heun, "units", "m"); CALL handle_err(stats)
    stats = nf90_put_att(lode, NF90_GLOBAL, "dx", dx); CALL handle_err(stats)

    stats = nf90_enddef(lode); CALL handle_err(stats)

    ! write nc
    stats = nf90_put_var(lode, var_x, x); CALL handle_err(stats)
    stats = nf90_put_var(lode, var_euler, err_eul); CALL handle_err(stats)
    stats = nf90_put_var(lode, var_heun, err_heu); CALL handle_err(stats)

    ! Time step
    DO n = 2, nx
        ! Use methods
        yeul = euler_forward(x(n-1), yeul, dx)
        yheu = heun(x(n-1), yheu, dx)

        ! Calculate errors
        yexc = sol_func(x(n))
        err_eul(n) = yeul - yexc
        err_heu(n) = yheu - yexc

        ! Save variables in netcdf
        stats = nf90_put_var(lode, var_euler, err_eul); CALL handle_err(stats)
        stats = nf90_put_var(lode, var_heun, err_heu); CALL handle_err(stats)

    END DO

    ! Close and deallocate
    stats = nf90_close(lode); CALL handle_err(stats)
    DEALLOCATE(x, err_eul, err_heu)

CONTAINS

    SUBROUTINE handle_err(stats)

        INTEGER, INTENT(IN) :: stats
        IF(stats/=nf90_noerr) THEN
            PRINT *, TRIM(nf90_strerror(stats))
            STOP "Stopped"
        END IF

    END SUBROUTINE handle_err

    REAL(8) FUNCTION func(x, y) RESULT(f)

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: x, y
        f = -0.5*y + 4.0*EXP(-0.5*x)*COS(4.0*x)

    END FUNCTION func

    REAL(8) FUNCTION sol_func(x) RESULT(exact_sol)

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: x
        exact_sol = EXP(-0.5*x)*SIN(4.0*x)

    END FUNCTION sol_func

    REAL(8) FUNCTION euler_forward(x, y, dx) RESULT(ynew)

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: x, y, dx
        ynew = y + dx*func(x, y)

    END FUNCTION euler_forward

    REAL(8) FUNCTION heun(x, y, dx) RESULT(ynew)

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: x, y, dx
        REAL(8) :: f0, yp, xp
        f0 = func(x, y)
        xp = x + dx
        yp = y + dx*f0
        ynew  = y + 0.5*dx*(f0 + func(xp, yp))

    END FUNCTION heun

END PROGRAM linear_ode
