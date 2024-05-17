program optimize_IPR
    use,intrinsic :: iso_fortran_env
    implicit none
    integer :: loop
    integer :: unit_write_result
    integer, parameter :: max_loop = 100
    integer, parameter :: N = 33
    integer, parameter :: N_c = 2
    integer, parameter :: number = 2
    logical, parameter :: is_to_update = .true.
    integer, parameter :: WF_levels(N_c) = [(13 + (loop - 1), loop = 1, N_c)]
    double precision, parameter :: eta = 0.1
    double precision :: WFs(N, N_c)
    double precision :: IPR, preIPR
    double precision :: gradients(N_c)
    double precision :: psis(N_c) = [((-1) / sqrt(real(N_c, kind(0d0))), loop = 1, N_c)]
    double precision :: old_psis(N_c, N_c) = 0

    WFs = 0
    old_psis = 0
    ! is_to_update = .true. または number > 1 なら必須.
    call read_psis()
    do loop = 1, N_c
        call read_wave_functions(loop)
    end do

    preIPR = 0
    do loop = 1, max_loop
        call get_IPR()
        call update_variants()
        print *, (IPR - preIPR) / IPR
        preIPR = IPR
    end do

    call write_files()
contains
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

    subroutine read_wave_functions(column)
        integer :: i, column
        double precision :: phi
        character(50) :: filename

        filename = "wave_functions/AB" // int32_to_string(N) // "/WF_" // int32_to_string(WF_levels(column)) // ".txt"
        open(newunit = unit_write_result, file = filename)
            do i = 1, N
                read(unit_write_result, *) phi
                WFs(i, column) = phi
            end do
        close(unit_write_result)
    end subroutine read_wave_functions

    subroutine get_IPR()
        integer :: i, j
        double precision :: sums(N)
        double precision :: temp(N)
        double precision :: gradient

        IPR = 0
        sums = 0
        do i = 1, N_c
            temp = WFs(1:N, i)
            sums = sums + psis(i) * temp
        end do
        if ( number > 1 ) then
            do i = 1, number - 1
                do j = 1, N_c
                    temp = WFs(1:N, j)
                    sums = sums - old_psis(j, i)
                end do
            end do
        end if

        IPR = sum(sums**4)

        do i = 1, N_c
            temp = WFs(1:N, i)
            gradient = 4 * sum(sums**3 * temp)
            gradients(i) = gradient
        end do
    end subroutine get_IPR

    subroutine update_variants()
        integer :: i
        
        do i = 1, N_c
            psis(i) = psis(i) + eta * gradients(i)
        end do

        ! normalization
        psis = psis * (1. / sqrt(sum(psis**2)))
    end subroutine update_variants

    subroutine write_files()
        integer :: i
        character(50) :: filename

        filename = "confined_states/AB" // int32_to_string(N) // "/psis_" // int32_to_string(number) // ".txt"
        open(newunit=unit_write_result, file=filename)
            do i = 1, N_c
                write(unit_write_result, *) psis(i)
            end do
        close(unit_write_result)
    end subroutine write_files

    subroutine read_psis()
        integer :: i, j
        double precision :: psi
        character(50) :: filename

        if ( is_to_update ) then
            filename = "confined_states/AB" // int32_to_string(N) // "/psis_" // int32_to_string(number) // ".txt"
            open(newunit=unit_write_result, file=filename)
                do i = 1, N_c
                    read(unit_write_result, *) psi
                    psis(i) = psi
                end do
            close(unit_write_result)
        end if

        if ( number > 1 ) then
            do i = 1, number - 1
                filename = "confined_states/AB" // int32_to_string(N) // "/psis_" // int32_to_string(i) // ".txt"
                open(newunit=unit_write_result, file=filename)
                do j = 1, N_c
                    read(unit_write_result, *) psi
                    old_psis(j, i) = psi
                end do
                close(unit_write_result)
            end do
        end if
    end subroutine read_psis
end program optimize_IPR