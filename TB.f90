program tight_binding
    use,intrinsic :: iso_fortran_env
    implicit none
    !parameter
    integer, parameter :: N = 153
    integer, parameter :: allhop = 280
    double precision, parameter :: tau = 1. + sqrt(2.)
    ! make_Hamiltonianのedge potentialの部分を領域D_nに応じて変更して.
    double precision, parameter :: edge_potential = 1.
    double precision, parameter :: pi = 4*atan(1.)
    character(7), parameter :: hop_file = "hop.txt"
    character(15), parameter :: coord_file = "coordinates.txt"
    character(24), parameter :: E_file = "eigval.txt"
    character(10), parameter :: result_file = "result.txt"
    !variable
    integer :: unit_write_result
    double precision :: Hamiltonian(N, N)
    ! for ZHEEVD of lapack
    integer :: INFO, IWORK(3 + 5 * N)
    double precision :: W(N)
    double precision :: WORK(1 + 6*N + 2*N**2)
    ! to clock time
    integer(int64) :: time_begin_c,time_end_c, CountPerSec, CountMax

    call system_clock(time_begin_c, CountPerSec, CountMax)

    call make_Hamiltonian()
    call DSYEVD("V", "U", N, Hamiltonian, N, W, WORK, 1 + 6*N + 2*N**2, IWORK, 3 + 5*N, INFO)

    call system_clock(time_end_c)

    call write_files()
contains
    subroutine make_Hamiltonian()
        integer :: i, j, k, l, m
        double precision :: x, y, Xs(N), Ys(N), H_ij

        Hamiltonian = 0
        open(newunit = unit_write_result, file = coord_file)
        do k = 1, N
            read(unit_write_result, *) x, y
            Xs(k) = x
            Ys(k) = y
        end do
        close(unit_write_result)

        ! edge potentials
        x = maxval(Xs)
        do k = 1, N
            !if ( sqrt(Xs(k)**2 + Ys(k)**2) > x - 10.**(-6) ) then
            if ( sqrt(Xs(k)**2 + Ys(k)**2) > sqrt((x - 1. / sqrt(2.))**2 + tau**2) + 10.**(-6) ) then
                i = k
                j = k
                H_ij = edge_potential

                Hamiltonian(i, j) = H_ij
            end if
        end do

        ! hopping elements
        open(newunit = unit_write_result, file = hop_file)
        do k = 1, allhop
            read(unit_write_result, *) l, m
            l = l + 1
            m = m + 1

            ! spin up
            i = l
            j = m
            H_ij = -1.

            Hamiltonian(i, j) = H_ij
            Hamiltonian(j, i) = H_ij
        end do
        close (unit_write_result)
    end subroutine make_Hamiltonian

    pure function int32_to_string(num) result(str)
        use, intrinsic :: iso_fortran_env
        implicit none

        integer(int32), intent(in) :: num
        character(:), allocatable :: str  ! retval

        integer(int32), parameter :: Extra_Digits = 0 ! 追加する桁数
        integer(int32) :: num_digits                 ! 整数numの桁数
        integer(int32) :: dgt_digits                 ! 整数numの桁数num_digitsの桁数
        character(:), allocatable :: str_digits
        character(:), allocatable :: fmt

        num_digits = get_digits_of(num) + Extra_Digits ! 整数numの桁数を取得
        dgt_digits = get_digits_of(num_digits)         ! 整数numの桁数num_digitsの桁数を取得

        ! 整数の桁数を文字列に変換
        allocate(character(dgt_digits)::str_digits)
        write(str_digits,'(I0)') num_digits

        ! 書式を作成
        fmt = "(I"//str_digits//"."//str_digits//")"

        ! 整数numを文字列に変換
        allocate(character(num_digits)::str)
        write(str,fmt) num
    end function int32_to_string

    pure function get_digits_of(num) result(num_digit)
        implicit none
        integer(int32),intent(in) :: num
        integer(int32) :: num_digit

        num_digit = int(log10(dble(num)))+1
    end function get_digits_of

    subroutine write_files()
        integer :: i, j

        do i = 1, N
            if ( abs(W(i)) < 10.**(-10) ) then
                open(newunit=unit_write_result, file="wave_functions/WF_" // int32_to_string(i) // ".txt")
                    do j = 1, N
                        write(unit_write_result, *) Hamiltonian(j, i)
                    end do
                close(unit_write_result)
            end if
        end do

        ! eigenenergies
        open(newunit=unit_write_result, file=E_file)
            do i = 1, N
                write(unit_write_result, *) W(i)
            end do
        close(unit_write_result)
    end subroutine write_files
end program tight_binding