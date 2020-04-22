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
    
    call bonusReihe()
    end program Interpolationspolynom
    
    function f(x)
        double precision :: f
        double precision, intent(in) :: x
        
        f = 3d0 / (x**2 + 1d0)
    end function
    subroutine bonusReihe()
        implicit none
        double precision :: x, treihe
        double precision, dimension(4) :: a
        double precision, dimension(3) :: b
        double precision, parameter :: x_step = 0.01d0, x_min = -1d1
        integer, parameter :: steps = 2000
        integer :: i, j
        double precision :: f
        
        
        open(unit=2, file="koeffizients.dat", action="write")
        !Berechnen der Koeffizienten mit der Trapezregel fürs Integral
        do j = 2, 4
            x = x_min
            do i = 1, steps - 1
                x = x_min + i * x_step
                a(j) = a(j) + (f(x)*cos(x * (j-1))+f(x+x_step)*cos(x * (j-1)))
                b(j-1) = b(j-1) + (f(x)*sin(x * (j-1))+f(x+x_step)*sin(x * (j-1)))
                if (j == 2) then
                    a(1) = a(1) + (f(x)+f(x+x_step)) / (2d0 * acos(-1d0)) * x_step
                end if
            end do
            if (j == 2) write(2,*) "a0 " , a(1) 
            a(j) = (a(j) / (2d0* acos(-1d0)))  * x_step
            b(j-1) = (b(j-1) / (2d0* acos(-1d0))) * x_step
            write(2,*) "a", (j-1) , " mit ", a(j), " b" , (j-1), " mit " , b(j-1)
        end do
        
    
        
        open(unit=1, file="TrigReihe.dat", action="write")
        !Bestimmen der Funktion
        do i = 1, steps !x werte
            x = x_min + i * x_step
            treihe = a(1) / 2d0
            do j = 2, 4 !n Werte
                treihe = treihe + a(j) * cos(x * (j-1)) + b(j-1) * sin(x* (j-1))
            end do
            write(1,*) x, treihe
        end do
        close(1)
    end subroutine bonusReihe
    
    
