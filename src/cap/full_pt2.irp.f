subroutine full_pt2()

    implicit none
    BEGIN_DOC
    ! Computes the second-order perturbation correction (PT2) to the CAP-CI energy.
    !
    ! Iterates over all generator determinants and their double-excitation spaces,
    ! evaluating the PT2 contribution using the CAP Hamiltonian H + i*eta*W.
    ! The total PT2 correction is accumulated in val and printed at the end.
    END_DOC

    complex*16 :: pt2(mo_num,mo_num), aha(mo_num,mo_num), val
    integer(bit_kind) :: Dg(N_int,2), Gh1h2(N_int,2), Di(N_int,2), Gh1h2p1p2(N_int,2), tmp(N_int,2)
    integer :: occ(N_int*bit_kind_size,2), vir(N_int*bit_kind_size,2)
    integer :: exc(0:2,2,2)
    double precision :: iha, phase
    integer :: i,j,k,l,g,h,p,r,s,idx, n_diff,nb
    integer :: n_occ(2), n_vir(2), degree, start_h2, start_p2
    integer :: N_g, N_s, h1,p1,h2,p2,s1,s2,spin,n_simple,n_double
    integer :: h1bis,p1bis,h2bis,p2bis,s1bis,s2bis
    logical :: ok, b(mo_num, mo_num)
    integer, allocatable :: iorder(:)
    double precision, external :: diag_H_mat_elem, diag_H_mat_elem_cap
    logical, external :: is_in_wavefunction
    double precision :: wij
    integer(bit_kind), allocatable :: alpha(:,:,:)

    N_g = N_det
    N_s = N_det
    nb = 0
    n_simple = 0
    n_double = 0
 
    val = (0d0,0d0)

    ! Generators
    do g = 1, N_g 
        Dg = psi_det(:,:,g)
        !print*,g
        ! Spin
        do s1 = 1, 2
            do s2 = s1, 2
                occ = 0
                call bitstring_to_list(Dg(1,1), occ(1,1), n_occ(1), N_int)
                call bitstring_to_list(Dg(1,2), occ(1,2), n_occ(2), N_int)
 
                ! hole 1 and 2
                do h1 = 1, n_occ(s1)
                    if (s1 == s2) then
                        start_h2 = h1 + 1
                    else
                        start_h2 = 1
                    endif
                    do h2 = start_h2, n_occ(s2)
                        if (h1 == h2 .and. s1 == s2) cycle
                        call apply_holes(Dg, s1, occ(h1,s1), s2, occ(h2,s2), Gh1h2, ok, N_int)
                        if (.not. ok) cycle
                        
                        ! List of virtual orbital for each spin of Gh1h2p1p2
                        vir = 0
                        do spin = 1, 2
                            l = 0
                            do k = 1, mo_num
                                call apply_particle(Gh1h2, spin, k, tmp, ok, N_int)
                                if (ok) then
                                    l += 1
                                    vir(l,spin) = k
                                endif
                                n_vir(spin) = l
                            enddo
                        enddo
  
                        pt2 = (0d0,0d0)
                        b = .False.

                        do p1 = 1, n_vir(s1)
                            if (s1 == s2) then
                                start_p2 = p1 + 1
                            else
                                start_p2 = 1
                            endif
                            do p2 = start_p2, n_vir(s2)
                                if (occ(h1,s1) == vir(p1,s1) .and. s1 /= s2 .and. min(n_occ(s1),n_occ(s2)) > 1) then
                                    b(p1,p2) = .True.
                                    cycle
                                endif
                                if (occ(h1,s1) == vir(p1,s1) .and. s1 /= s2 .and. min(n_occ(s1),n_occ(s2)) == 1 .and. occ(h1,s1) /= 1) then
                                    b(p1,p2) = .True.
                                    cycle
                                endif
                                if (occ(h2,s2) == vir(p2,s2) .and. s1 /= s2) then
                                    b(p1,p2) = .True.
                                    cycle
                                endif

                                if (occ(h1,s1) == vir(p1,s1) .and. s1 == s2 .and. h1 /= 1) then
                                    b(p1,p2) = .True.
                                    cycle
                                endif
                                if (occ(h2,s2) == vir(p1,s1) .and. s1 == s2 .and. h2 /= 2) then
                                    b(p1,p2) = .True.
                                    cycle
                                endif
                                if (occ(h1,s1) == vir(p2,s2) .and. s1 == s2 .and. h1 /= 1) then
                                    b(p1,p2) = .True.
                                    cycle
                                endif
                                if (occ(h2,s2) == vir(p2,s2) .and. s1 == s2 .and. h2 /= 2) then
                                    b(p1,p2) = .True.
                                    cycle
                                endif

                                if (p1 == p2 .and. s1 == s2) then
                                    b(p1,p2) = .True.
                                    cycle
                                endif
                            enddo
                        enddo

                        do i = 1, N_s
                            Di = psi_det(:,:,i)

                            n_diff = 0
                            do k = 1, 2
                                do j = 1, N_int
                                    n_diff += popcnt(IEOR(Di(j,k),Gh1h2(j,k)))
                                enddo
                            enddo

                            if (n_diff > 6) then
                                cycle
                            endif

                            do p1 = 1, n_vir(s1)
                                if (s1 == s2) then
                                    start_p2 = p1 + 1
                                else
                                    start_p2 = 1
                                endif
                                do p2 = start_p2, n_vir(s2)
                                    if (b(p1,p2)) cycle 

                                    call apply_particles(Gh1h2, s1, vir(p1,s1), s2, vir(p2,s2), Gh1h2p1p2, ok, N_int)
                                    if (.not. ok) cycle
                                    call get_excitation_degree(Di, Gh1h2p1p2, degree, N_int)

                                    if (degree == 0 .or. (degree <= 2 .and. i < g)) then
                                        ! if degree <= 2 and i < N_g, Ghp have been generated
                                        ! by another generator
                                        b(p1,p2) = .True.
                                        cycle
                                    else if (degree <= 2)  then
                                        ! Gh1h2p1p2 is connected to Di
                                        if (degree == 1) then
                                            call get_excitation(Di, Gh1h2p1p2, exc, degree, phase, N_int)
                                            call decode_exc(exc,degree,h1bis,p1bis,h2bis,p2bis,s1bis,s2bis)
                                            call i_W_j_single_spin_cap(Di,Gh1h2p1p2,N_int,s1bis,wij)
                                        else
                                            wij = 0d0
                                        endif
                                        !wij = 0d0
                                        call i_H_j(Di, Gh1h2p1p2, N_int, iha)
                                        pt2(p1,p2) += psi_cap_coef(i,1) * dcmplx(iha,wij)
                                    endif
                                enddo
                            enddo
                        enddo

                        do p1 = 1, n_vir(s1)
                            if (s1 == s2) then
                                start_p2 = p1 + 1
                            else
                                start_p2 = 1
                            endif
                            do p2 = start_p2, n_vir(s2)
                                if (b(p1,p2)) cycle
                                call apply_particles(Gh1h2, s1, vir(p1,s1), s2, vir(p2,s2), Gh1h2p1p2, ok, N_int)
                                aha(p1,p2) = dcmplx(diag_H_mat_elem(Gh1h2p1p2, N_int),diag_H_mat_elem_cap(Gh1h2p1p2, N_int))
                            enddo
                        enddo

                        do p1 = 1, n_vir(s1)
                            if (s1 == s2) then
                                start_p2 = p1 + 1
                            else
                                start_p2 = 1
                            endif
                            do p2 = start_p2, n_vir(s2)
                                if (b(p1,p2)) cycle
                                call apply_particles(Gh1h2, s1, vir(p1,s1), s2, vir(p2,s2), Gh1h2p1p2, ok, N_int)
                                pt2(p1,p2) = pt2(p1,p2)**2 / (psi_energy_cap(1) - aha(p1,p2))
                                val += pt2(p1,p2)
                                nb += 1
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
    print*,'val', val,nb

end
