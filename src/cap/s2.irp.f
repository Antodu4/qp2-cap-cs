subroutine u_0_S2_u_0_complex(e_0,u_0,n,keys_tmp,Nint,N_st,sze_8)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes e_0 = <u_0|S2|u_0>/<u_0|u_0>
  !
  ! n : number of determinants
  !
  END_DOC
  integer, intent(in)            :: n,Nint, N_st, sze_8
  complex*16, intent(out)  :: e_0(N_st)
  complex*16, intent(in)   :: u_0(sze_8,N_st)
  integer(bit_kind),intent(in)   :: keys_tmp(Nint,2,n)

  complex*16, allocatable  :: v_0(:,:)
  complex*16 :: uv,uu
  integer :: i,j
  allocate (v_0(sze_8,N_st))

  call S2_u_0_nstates_complex(v_0,u_0,n,keys_tmp,Nint,N_st,sze_8)
  do i=1,N_st
    call inner_product_complex(u_0(1,i),v_0(1,i),n,uv)
    call inner_product_complex(u_0(1,i),u_0(1,i),n,uu)
    e_0(i) = uv/uu
  enddo
end

subroutine S2_u_0_nstates_complex(v_0,u_0,n,keys_tmp,Nint,N_st,sze_8)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v_0  = S^2|u_0>
  !
  ! n : number of determinants
  !
  END_DOC
  integer, intent(in)            :: N_st,n,Nint, sze_8
  complex*16, intent(out)  :: v_0(sze_8,N_st)
  complex*16, intent(in)   :: u_0(sze_8,N_st)
  integer(bit_kind),intent(in)   :: keys_tmp(Nint,2,n)
  double precision               :: s2_tmp
  complex*16, allocatable  :: vt(:,:)
  integer                        :: i,j,k,l, jj,ii
  integer                        :: i0, j0

  integer, allocatable           :: shortcut(:,:), sort_idx(:,:)
  integer(bit_kind), allocatable :: sorted(:,:,:), version(:,:,:)
  integer(bit_kind)              :: sorted_i(Nint)

  integer                        :: sh, sh2, ni, exa, ext, org_i, org_j, endi, istate


  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (n>0)
  PROVIDE ref_bitmask_energy

  allocate (shortcut(0:n+1,2), sort_idx(n,2), sorted(Nint,n,2), version(Nint,n,2))
  v_0 = dcmplx(0.d0,0d0)

  call sort_dets_ab_v(keys_tmp, sorted(1,1,1), sort_idx(1,1), shortcut(0,1), version(1,1,1), n, Nint)
  call sort_dets_ba_v(keys_tmp, sorted(1,1,2), sort_idx(1,2), shortcut(0,2), version(1,1,2), n, Nint)

  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP PRIVATE(i,s2_tmp,j,k,jj,vt,ii,sh,sh2,ni,exa,ext,org_i,org_j,endi,sorted_i,istate)&
      !$OMP SHARED(n,u_0,keys_tmp,Nint,v_0,sorted,shortcut,sort_idx,version,N_st,sze_8)
  allocate(vt(sze_8,N_st))
  vt = dcmplx(0.d0,0d0)

  do sh=1,shortcut(0,1)
    !$OMP DO SCHEDULE(static,1)
    do sh2=sh,shortcut(0,1)
      exa = 0
      do ni=1,Nint
        exa = exa + popcnt(xor(version(ni,sh,1), version(ni,sh2,1)))
      end do
      if(exa > 2) then
        cycle
      end if

      do i=shortcut(sh,1),shortcut(sh+1,1)-1
        org_i = sort_idx(i,1)
        if(sh==sh2) then
          endi = i-1
        else
          endi = shortcut(sh2+1,1)-1
        end if
        do ni=1,Nint
          sorted_i(ni) = sorted(ni,i,1)
        enddo

        do j=shortcut(sh2,1),endi
          org_j = sort_idx(j,1)
          ext = exa
          do ni=1,Nint
            ext = ext + popcnt(xor(sorted_i(ni), sorted(ni,j,1)))
          end do
          if(ext <= 4) then
            call get_s2(keys_tmp(1,1,org_i),keys_tmp(1,1,org_j),Nint,s2_tmp)
            do istate=1,N_st
              vt (org_i,istate) = vt (org_i,istate) + dcmplx(s2_tmp,0d0)*u_0(org_j,istate)
              vt (org_j,istate) = vt (org_j,istate) + dcmplx(s2_tmp,0d0)*u_0(org_i,istate)
            enddo
          endif
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT
  enddo

  do sh=1,shortcut(0,2)
    !$OMP DO
    do i=shortcut(sh,2),shortcut(sh+1,2)-1
      org_i = sort_idx(i,2)
      do j=shortcut(sh,2),i-1
        org_j = sort_idx(j,2)
        ext = 0
        do ni=1,Nint
          ext = ext + popcnt(xor(sorted(ni,i,2), sorted(ni,j,2)))
        end do
        if(ext == 4) then
          call get_s2(keys_tmp(1,1,org_i),keys_tmp(1,1,org_j),Nint,s2_tmp)
          do istate=1,N_st
            vt (org_i,istate) = vt (org_i,istate) + dcmplx(s2_tmp,0d0)*u_0(org_j,istate)
            vt (org_j,istate) = vt (org_j,istate) + dcmplx(s2_tmp,0d0)*u_0(org_i,istate)
          enddo
        end if
      end do
    end do
    !$OMP END DO NOWAIT
  enddo
  !$OMP BARRIER

  do istate=1,N_st
    do i=n,1,-1
      !$OMP ATOMIC
      v_0(i,istate) = v_0(i,istate) + vt(i,istate)
    enddo
  enddo

  deallocate(vt)
  !$OMP END PARALLEL

  do i=1,n
    call get_s2(keys_tmp(1,1,i),keys_tmp(1,1,i),Nint,s2_tmp)
    do istate=1,N_st
      v_0(i,istate) += dcmplx(s2_tmp,0d0) * u_0(i,istate)
    enddo
  enddo

  deallocate (shortcut, sort_idx, sorted, version)
end

