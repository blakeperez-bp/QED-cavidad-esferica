Module funciones_index

implicit none

contains
    integer function Index_photonic(tupla, nmax) result(indice)
        implicit none
        integer, intent(in) :: tupla(:)
        integer, intent(in) :: nmax
        integer :: m,k
        indice = 0
        m = size(tupla)
        do k = 1, m
            indice = indice + tupla(k) * (nmax+1)**(k-1)
        end do
    end function Index_photonic

    function fock_photonic(indice, nmax, vector_size) result(tupla)  !!hay que pasarle el indice, el nro max de fotones y la cantidad de modos
        implicit none
        integer, intent(in) :: indice, nmax, vector_size
        integer :: tupla(vector_size)
        integer :: b, N, r, j, digit_count
        
        ! Inicializar
        tupla = 0
        N = indice
        b = nmax + 1
        digit_count = vector_size
        
        ! Conversión CORRECTA a base b
        do j = vector_size, 1, -1
            r = mod(N, b)
            tupla(j) = r
            N = N / b
            if (N == 0) exit
        end do
        
    end function fock_photonic


    integer function encode_tuple(xs, ranges, nvars) result(ind)
        implicit none
        integer, intent(in) :: xs(:)
        integer, intent(in) :: ranges(:,:)
        integer, intent(in) :: nvars
        integer :: i
        integer :: a, b, d, y
        integer :: mult

        ind = 0
        mult = 1
        do i = 1, nvars
            a = ranges(i,1)
            b = ranges(i,2)
            d = b - a + 1
            y = xs(i) - a

            if (xs(i) < a .or. xs(i) > b) then
                print*, "Error en la variable ", i
                stop
            end if

            ind = ind + y * mult
            mult = mult * d  ! actualizar mult para la siguiente variable
        end do
    end function encode_tuple

    subroutine decode_index(ind, ranges, nvars, xs)
    implicit none
    integer, intent(in) :: ind
    integer, intent(in) :: ranges(:,:)
    integer, intent(in) :: nvars
    integer, intent(out) :: xs(:)
    integer :: i, a, b, d, y, idx, prod_d

    ! Copiar el índice
    idx = ind

    ! Calcular el producto de los rangos
    prod_d = 1
    do i = 1, nvars
        prod_d = prod_d * (ranges(i,2) - ranges(i,1) + 1)
    end do

    ! Verificar rango
    if (ind < 0 .or. ind >= prod_d) then
        print*, "Error: indice fuera de rango"
        stop
    end if

    ! Decodificar
    do i = 1, nvars
        a = ranges(i,1)
        b = ranges(i,2)
        d = b - a + 1
        y = mod(idx,d)
        xs(i) = y + a
        idx = idx / d
    end do
    end subroutine decode_index


    !!funcion que checkea si dos vectores cumplen que en el lugar j-esimo difieren en 1 y si el resto son iguales
    !!va a servir para hacer el elemento <N\a_j\N'> , si esas dos condiciones no se cumplen, todo el elemento es cero
    !!da el resultado apenas ve que falla una condicion, es mas eficiente
    !!IMPORTANTE a va por derecha y b por izquierda para el operador destruccion
    !!devuelve la raiz del numero de ocupacion

    real*8 function delta_a(a, b, j) result(res)
        implicit none
        integer, intent(in) :: a(:), b(:)   ! vectores de enteros
        integer, intent(in) :: j
        integer :: i

        ! Condición en el índice j
        if (a(j) - b(j) /= 1) then
            res = 0.d0
            return
        end if

        ! Condición en el resto
        do i = 1, size(a)
            if (i /= j) then
                if (a(i) - b(i) /= 0) then
                    res = 0.d0
                    return
                end if
            end if
        end do

        res = sqrt(a(j)*1.d0)
    end function delta_a

    real*8 function delta_adaga(a, b, j) result(res) !!para el operador creacion
        implicit none
        integer, intent(in) :: a(:), b(:)   ! vectores de enteros
        integer, intent(in) :: j
        integer :: i

        ! Condición en el índice j
        if (b(j) - a(j) /= 1) then
            res = 0.d0
            return
        end if

        ! Condición en el resto
        do i = 1, size(a)
            if (i /= j) then
                if (a(i) - b(i) /= 0) then
                    res = 0.d0
                    return
                end if
            end if
        end do

        res = sqrt(a(j)*1.d0 + 1.d0)
    end function delta_adaga



end module funciones_index


