PROGRAM linear_advection

    IMPLICIT NONE
    REAL(8) :: u, c, dt, dx, x1, x0, t0, t1, tp, alpha, beta
    REAL(8), DIMENSION(:), ALLOCATABLE :: x, phi_new, phi_now, phi_old, d
    REAL(8), DIMENSION(:, :), ALLOCATABLE :: phi_plot
    INTEGER :: n, t, i, nx, nt, ios, w, itp, itp_n, itp_s

    !Initial conditions
    u = -0.31
    dx = 0.1
    x0 = 0.0
    x1 = 500.0
    t0 = 0.0
    t1 = 1001.0
    tp = 200.0
    alpha = 0.05
    beta = 0.53

    dt = 0.1
    c = u * dt / dx
    nx = INT((x1 - x0) / dx) + 1
    nt = INT((t1 - t0) / dt) + 1

    itp = NINT(tp / dt)
    itp_n = INT(t1 / tp) + 1
    itp_s = 1

    ALLOCATE(x(nx), phi_now(nx), phi_old(nx), phi_new(nx), phi_plot(nx, itp_n), d(nx))

    DO n = 1, nx
        x(n) = x0 + (n-1)*dx
        phi_old(n) = phi0(x(n))
    END DO

    ! Initial conditions for arrays
    phi_now(:) = phi_old(:)
    d(:) = phi_old(:)
    phi_plot(:, 1) = phi_old(:)

    ! To check
    !DO n = 1, nx
    !    PRINT *, x(n), phi(n, 1)
    !END DO

    IF (u > 0) THEN
        phi_now = ftbs(phi_old, c)
    ELSE
        phi_now = ftfs(phi_old, c)
    END IF

    DO t = 2, nt
        phi_new = ctcs(phi_old, phi_now, c)

        ! Non filter
        !phi_old(:) = phi_now(:)
        !phi_now(:) = phi_new(:)

        d = alpha*(phi_old + phi_new - 2.0*phi_now)

        ! RA filter
        !phi_old(:) = phi_now + d
        !phi_now(:) = phi_new

        ! RAW filter
        phi_old(:) = phi_now + beta*d
        phi_now(:) = phi_new + (1 - beta)*d

        IF (MOD(t-1, itp) == 0) THEN
            itp_s = itp_s + 1
            IF (itp_s <= itp_n) THEN
                phi_plot(:, itp_s) = phi_now(:)
            END IF
        END IF


    END DO

    ! Save file
    OPEN(NEWUNIT=w, IOSTAT=ios, FILE='advection_ctcs_RAW_filter.txt', &
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
