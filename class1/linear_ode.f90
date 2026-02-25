PROGRAM linear_ode

    IMPLICIT NONE
    REAL(8) :: x0, x1, dx, y0
    INTEGER :: nx, n, u, ios
    REAL(8), ALLOCATABLE :: x(:), err_eu(:), err_he(:)
    REAL(8) :: y_eu, y_he, y_sol

    ! Parameters
    x0 = 0.0
    x1 = 10.0
    dx = 0.1
    y0 = 0.0

    nx = INT((x1 - x0)/dx) + 1

    ALLOCATE(x(nx), err_eu(nx), err_he(nx))

    ! Build x grid
    DO n = 1, nx
        x(n) = x0 + DBLE(n-1)*dx
    END DO

    ! Initial conditions
    y_eu = y0
    y_he = y0
    y_sol = sol_func(x(1))
    err_eu(1)  = y_eu - y_sol
    err_he(1)  = y_he - y_sol

    ! Time step
    DO n = 2, nx
        y_eu = euler_forward(x(n-1), y_eu, dx)
        y_he = heun(x(n-1), y_he, dx)
        y_sol = sol_func(x(n))
        err_eu(n) = y_eu - y_sol
        err_he(n) = y_he - y_sol
    END DO

    ! Save file
    OPEN(NEWUNIT=u, IOSTAT=ios, FILE='euler_heun.txt', &
      STATUS='new', ACTION='write')

        IF (ios == 0) THEN
            WRITE(u,'(A)') "# x  err_euler  err_heun"
            DO n = 1, nx
                WRITE(u,'(3F20.12)') x(n), err_eu(n), err_he(n)
            END DO
            CLOSE(u)
        END IF

    DEALLOCATE(x, err_eu, err_he)

CONTAINS

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
