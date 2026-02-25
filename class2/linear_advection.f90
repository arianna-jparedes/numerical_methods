PROGRAM linear_advection

    IMPLICIT NONE
    REAL(8) :: u, c, dt, dx, x1, x0, t0, t1, tp
    REAL(8), DIMENSION(:), ALLOCATABLE :: x, phi, phi_new
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: phi_plot
    INTEGER :: n, t, i, nx, nt, ios, w, itp, itp_n, itp_s

    !Initial conditions
    u = 0.087
    dx = 0.1
    x0 = 0.0
    x1 = 100.0
    t0 = 0.0
    t1 = 1001.0
    tp = 200.0

    !dt = dx / u
    dt = 0.5
    c = u * dt / dx
    nx = INT((x1 - x0) / dx) + 1
    nt = INT((t1 - t0) / dt) + 1

    itp = INT(tp / dt)
    itp_n = INT(t1 / tp) + 1
    itp_s = 1

    ALLOCATE(x(nx), phi_new(nx), phi(nx), phi_plot(nx, itp_n))

    DO n = 1, nx
        x(n) = x0 + (n-1)*dx
        phi(n) = phi0(x(n))
        phi_plot(n, 1) = phi(n)
    END DO

    ! To check
    !DO n = 1, nx
    !    PRINT *, x(n), phi(n, 1)
    !END DO

    DO t = 2, nt

        IF (u > 0) THEN
            phi_new = ftbs(phi, c)
        ELSE
            phi_new = ftfs(phi, c)
        END IF
        phi(:) = phi_new(:)

        IF (MOD(t-1, itp) == 0) THEN
            itp_s = itp_s + 1
            IF (itp_s <= itp_n) THEN
                phi_plot(:, itp_s) = phi(:)
            END IF
        END IF


    END DO

    ! Save file
    OPEN(NEWUNIT=w, IOSTAT=ios, FILE='advection.txt', &
      STATUS='new', ACTION='write')

        IF (ios == 0) THEN
            DO n = 1, nx
                WRITE(w, *) phi_plot(n, 1:itp_s)
            END DO
            CLOSE(w)
        END IF


CONTAINS

    REAL(8) FUNCTION phi0(x) RESULT(phi)

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: x
        REAL(8) :: pi

        pi = 4.0 * ATAN(1.0)

        phi = 0.0
        IF (x >= 40 .AND. x <= 70) phi = SIN(pi*(x - 40)/30)**2

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

END PROGRAM linear_advection
