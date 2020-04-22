!  Interpolationspolynom.f90 
!
!  FUNCTIONS:
!  Interpolationspolynom - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Interpolationspolynom
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program Interpolationspolynom

    implicit none
    double precision, dimension(:,:), allocatable :: fa
    double precision, dimension(:), allocatable :: a, basis
    double precision :: x, fx = 0
    double precision, parameter :: x_step = 1d-2, x_min = -6d0
    integer :: n = 0, stat, i, j, k
    integer, parameter :: steps = 1300
    open(unit=1, file="int-pol.txt", action="read")
    do
        read(1,*,iostat = stat)
        if (stat /= 0) exit
        n = n +1      
    end do
    close(1)
    allocate(fa(n,n))
    allocate(a(n))
    allocate(basis(n))
    basis(1) = 1d0
    !Koeffizienten berechnen
    open(unit=1, file="int-pol.txt", action="read") 
    do i = 1, n
        read(1,*) a(i), fa(i,1)
        do j = 2, i
            fa(i,j) = (fa(i-1, j-1) - fa(i, j-1)) / (a(i - j + 1) - a(i))
        end do
    end do
    close(1)
    
    do i = 1, size(fa, dim=1)
            do j= 1, size(fa, dim=2)
                write(*,'(F12.6)', advance="no") fa(i,j)
            end do
            print*,""
    end do
    
    open(unit = 2, file="polynom.dat", action="write")
    do k = 0, steps !gibt uns den x Wert
        x = k * x_step + x_min
        fx = fa(1,1)
        do i = 2, n !geht durch die Anzahl der Daten, gleichbedeutend mit der Zeilenanzahl (runter)
            !do j = 2, i !geht ín der Spalte bis zur aktuellen Zeilenzahl               
            basis(i) = basis(i-1) * (x - a(i-1))                            
            fx = fx + basis(i) *fa(i,i)
            !end do
        end do
        write (2,*) x, fx
        fx = 0d0
        basis(1) = 1d0
    end do
    deallocate(fa)
    deallocate(a)
    deallocate(basis)
    end program Interpolationspolynom
