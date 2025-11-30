Module ModuleState
use ModuleParameters
use funciones_index
use sph_bessel_roots

implicit none


integer, allocatable :: state_matrix_m(:,:)
integer, allocatable :: state_matrix_ph(:,:)
integer, allocatable :: modes_matrix(:,:)

contains
    subroutine swap(a, b)
        integer, intent(inout) :: a, b
        integer :: tmp
        tmp = a
        a = b
        b = tmp
    end subroutine swap



    subroutine InitConfSpace() !--> devuelvo las matrices de estados
        use ModuleParameters
        use funciones_index
        implicit none

        integer, allocatable :: te(:), tm(:)
        integer :: i, ii , kk, jj, m, j
        integer :: inMA, ilMA, imMA, alpha, il_ph, im_ph, in_ph
        integer :: nstates, ncols
        integer, allocatable :: total_photons(:), idx(:)
        integer, allocatable :: temp_matrix(:,:)



        !!idea: puedo poner una matriz de Nsystem filas y con 4 columnas
        !la 1era sera un numero que codifique el estado de fock
        !las otras 3 seran nlm materiales
        !!--> de esta forma la fila i de la matriz tiene toda la informacion


        write(*,*) "InitConfSpace"

        allocate(state_matrix_m(NsysMA,3))
        allocate(state_matrix_ph(NsysPH, M_TE + M_TM ))
        allocate(modes_matrix(M_TE + M_TM,4))
        allocate(te(M_TE), tm(M_TM))


        write(*,'(A)') "ESTADOS MATERIALES"
        !!!matriz de estados materiales
        i = 1
        Do ilMA = 0, Lmax  !!parte material
            Do imMA = -ilMA, ilMA
                Do inMA = 1, Nmax

                    ! write(*,*) i
                    write(*,'(A)') "----------------------------------"
                    write(*,*) ilMA, imMA, inMA

                    state_matrix_m(i,1) = ilMA
                    state_matrix_m(i,2) = imMA
                    state_matrix_m(i,3) = inMA

                    i = i + 1
                End Do
            End Do
        End Do



        write(*,'(A)') "ESTADOS FOTONICOS"

        Do ii = 1, NsysPH_TE
            Do jj = 1, NsysPH_TM
                te =  fock_photonic(ii-1 , Nph_TE, M_TE)
                tm =  fock_photonic(jj-1 , Nph_TM, M_TM)
                m = NsysPH_TM*(ii-1) + jj !!indice de la matriz
                state_matrix_ph(m,:) = [te,tm]  !!concateno configuraciones de Fock, fijo una TE y barro por las TM
!                 write(*,*) state_matrix_ph(m,:)
                ! write(*,*) (state_matrix_ph(m,j), j=1,size(state_matrix_ph,2))
                ! write(*,'(A)') "-------------------------------------------------------------------------"
            End Do
        End Do

        !!!!!!!!!!!!!!!!!!! ordeno la matriz segun numero de fotones !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        nstates = size(state_matrix_ph, 1)
        ncols   = size(state_matrix_ph, 2)

        allocate(total_photons(nstates))
        allocate(idx(nstates))

        ! Calculo el numero total de fotones en cada estado (suma de la fila)
        do i = 1, nstates
            total_photons(i) = sum(state_matrix_ph(i, :))
            idx(i) = i
        end do

        ! Ordeno indices idx segun total_photons
        do i = 1, nstates - 1
            do j = i + 1, nstates
                if (total_photons(j) < total_photons(i)) then
                    call swap(total_photons(i), total_photons(j))
                    call swap(idx(i), idx(j))
                end if
            end do
        end do

        ! Reordeno la matriz segun los indices ordenados
        allocate(temp_matrix(nstates, ncols))
        do i = 1, nstates
            temp_matrix(i, :) = state_matrix_ph(idx(i), :)
        end do

        state_matrix_ph = temp_matrix

        !!para verfificar        
        write(*,'(A)') "Matriz ordenada por numero total de fotones:"
        do i = 1, nstates
            write(*,'(A)') "-------------------------------"
            write(*,'(100I2)') state_matrix_ph(i,:)
        end do


        deallocate(temp_matrix, total_photons, idx)




        !!Ahora tengo dos matrices --> una con los estados de materia y otra con las configuraciones de fock
        !!es como si tuviera una matriz de NsysPH X NsysMA con Nsystem en total
        !!filas --> fotones
        !!columna --> materia
        !!Indice j se puede construir como J = NsysMA*(i-1) + ii


        !!!forma en la que ordeno los estados materiales
        !! indice que recorre los estados desde 1 hasta el numero de modos fotonicos

        ! kk = 1
        ! Do alpha = 1,2
        !     if (alpha == 1) then !TE
        !         Do il_ph = 1, Lmodes_TE 
        !             Do im_ph = -il_ph, il_ph
        !                 Do in_ph = 1, Rmodes_TE

        !                 ! write(*,*) kk
        !                 ! write(*,'(A)') "----------------------------------"
        !                 write(*,*) alpha, il_ph, im_ph, in_ph

        !                 modes_matrix(kk,1) = alpha
        !                 modes_matrix(kk,2) = il_ph
        !                 modes_matrix(kk,3) = im_ph
        !                 modes_matrix(kk,4) = in_ph
        !                 kk = kk+1
        !                 End Do
        !             End Do
        !         End Do

        !     else  !TM
        !         Do il_ph = 1, Lmodes_TM 
        !             Do im_ph = -il_ph, il_ph
        !                 Do in_ph = 1, Rmodes_TM

        !                 ! write(*,*) kk
        !                 ! write(*,'(A)') "----------------------------------"
        !                 write(*,*) alpha, il_ph, im_ph, in_ph

        !                 modes_matrix(kk,1) = alpha
        !                 modes_matrix(kk,2) = il_ph
        !                 modes_matrix(kk,3) = im_ph
        !                 modes_matrix(kk,4) = in_ph
        !                 kk =kk+1
        !                 End Do
        !             End Do
        !         End Do
        !     End If
        ! End Do
        !!!con modificacion de m0
        kk = 1
        Do alpha = 1,2
            if (alpha == 1) then !TE
                Do il_ph = 1, Lmodes_TE 
                    Do im_ph = -min(il_ph,m0), min(il_ph,m0)
                        Do in_ph = 1, Rmodes_TE

                        ! write(*,*) kk
                        ! write(*,'(A)') "----------------------------------"
                        write(*,*) alpha, il_ph, im_ph, in_ph

                        modes_matrix(kk,1) = alpha
                        modes_matrix(kk,2) = il_ph
                        modes_matrix(kk,3) = im_ph
                        modes_matrix(kk,4) = in_ph
                        kk = kk+1
                        End Do
                    End Do
                End Do

            else  !TM
                Do il_ph = 1, Lmodes_TM 
                    Do im_ph = -min(il_ph,m0), min(il_ph,m0)
                        Do in_ph = 1, Rmodes_TM

                        ! write(*,*) kk
                        ! write(*,'(A)') "----------------------------------"
                        write(*,*) alpha, il_ph, im_ph, in_ph

                        modes_matrix(kk,1) = alpha
                        modes_matrix(kk,2) = il_ph
                        modes_matrix(kk,3) = im_ph
                        modes_matrix(kk,4) = in_ph
                        kk =kk+1
                        End Do
                    End Do
                End Do
            End If
        End Do
        !!estan ordenados primero todos los TE y despues los TM

    End subroutine InitConfSpace



!!subrutina que dado un valor de k devuelve los estados de Fock y materiales
    subroutine get_state(k,NsysMAT, conf_fock, est_material)!!usa las matrices de estados que ya estan en module parameters
        use ModuleParameters
        implicit none

        integer, intent(in) :: k, NsysMAT
        integer, intent(out), allocatable :: conf_fock(:), est_material(:)
        integer :: fila, columna

        !!ya reservo la memoria de salida
        allocate(conf_fock(M_TE+M_TM))
        allocate(est_material(3))

        !    fila = ((k - mod(k,NsysMAT))/NsysMAT) + 1 !!--> hacer esto me dice en que fila cae
        !    conf_fock = state_matrix_ph(fila,:)


        !    if (mod(k,NsysMAT) == 0) then          !!--> hacer esto me dice en que columna cae
        !       columna = NsysMAT             
        !    else 
        !       columna = mod(k,NsysMAT)
        !    end if
        !    est_material = state_matrix_m(columna,:)
   
        fila    = (k - 1) / NsysMAT + 1
        columna = mod(k - 1, NsysMAT) + 1
        conf_fock = state_matrix_ph(fila,:)
        est_material = state_matrix_m(columna,:)


    end subroutine get_state


    !!para calcular la energia
    real*8 function E_N(vector_fock) result(energia_N)
        integer, intent(in):: vector_fock(:) !!esto va a ser una de las filas de state_matrix_ph, tiene los nros de ocupacion
        integer :: i,l,n
        energia_N = 0.d0
        do i = 1,M_TE
            l = modes_matrix(i,2)
            n = modes_matrix(i,4)
            energia_N = energia_N + 137.d0*(1.d0/Rmax)*bessel_j_root(l,n)*vector_fock(i)
        End do
        
        do i = M_TE+1, M_TE + M_TM
            l = modes_matrix(i,2)
            n = modes_matrix(i,4)
            energia_N = energia_N + 137.d0*(1.d0/Rmax)*bessel_jprime_root(l,n)*vector_fock(i)
        End do

    end function E_N

End Module ModuleState
