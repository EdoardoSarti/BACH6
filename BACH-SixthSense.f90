program tt
implicit none

  integer, parameter :: NLINE_max=50000
  integer, parameter :: NCHAIN_max=25
  integer, parameter :: NRES_max=1500
  integer, parameter :: DBIAS=4
  integer, parameter :: CBIAS=1
  real*8, parameter  :: R_PROBE=1.4d0
  real*8, parameter  :: RESDIST = 17.d0
  real*8, parameter  :: threshold = 8.d0
  real*8, parameter  :: pi = 3.1415926
  real*8, parameter  :: CLASH_THR = 6.d0
  real*8, parameter  :: VERLET_RADIUS= 20.d0
  real*8, parameter  :: INTERFACE_THRESHOLD= 10.d0
  real*8, parameter  :: R_interaction= 4.5d0

  integer                  :: NXRES_max = 0 , NXCHAIN_max=10, NXNEIGH_max = 200
  integer                  :: argcount,IARGC,nrestot,namelen
  integer                  :: NumRes(NCHAIN_max),holes(NCHAIN_max)
  integer                  :: NCHAIN,MAX_CODE
  integer, allocatable     :: numresneigh(:),nneigh(:,:),resneigh(:,:),neighborhood(:,:,:,:)
  integer                  :: apolar(20,14),nexc(21,14)
  integer, allocatable     :: ncmap(:,:,:,:),nbridgecmap(:,:,:,:),nhor_cp(:),nti_cp(:),strictlist(:,:)
  integer                  :: i,j,ii,NL,NM,ii_old,j_old
  real*8, allocatable      :: h_pos(:,:,:),atom_pos(:,:,:,:),ap_cp(:,:,:)
  real*8                   :: vGUT(20,20,0:5),vSAS(20,2),eneSOLV(2),es_int,resareamax(20)
  real*8                   :: param(4),PARAM_ALL(4,10),atomcode(21,10),sasaparam(21,14,4)
  real*8                   :: eneGUT(0:5),eneGUT_int(0:5)
  real*8, allocatable      :: resarea(:),sum_1(:,:),sum_2(:,:),sum_3(:,:),area(:,:),radius(:,:),clash_param(:,:),thr_couple(:)
  real*8                   :: normt,normaa,en_clashes(2)
  real*8                   :: ntot(20),nconct(20,2)
  character*3, allocatable :: ResName(:,:)
  character*3              :: order(21,14),order_p(21,14),exceptions(21,14,10),exception
  character*5              :: check
  character*45             :: pdbname45
  character*150            :: pdbname,wq_char,filepar,fileparat,filepdb,outputfile,line_r(NLINE_max),line,prefix
  logical                  :: ONLY_GUT=.false., QUIET=.false.,QUIET_W=.false., INTF=.false., S_INTF=.false., COMPUTE_PARM=.false.
  logical                  :: TRAJ=.false., PDBLIST=.false., READ_MODE=.false., COMPUTE_ENE=.false., SEPARATE_TERMS=.false.

  param(:)=0.d0
  nexc(:,:)=0
  argcount = IARGC()
  filepdb='XXX'
  filepar='BSS.par'
  fileparat='ATOMIC_PARAMETERS_BSS'
  prefix='mod_'
  outputfile='output.bss'
  if(argcount==0)then
      write(6,*)'USAGE:'
      write(6,*)'bach++.x -COMPUTE_PAR xor -COMPUTE_ENE : learns parameters XOR computes eneregies'
      write(6,*)'         -PDBLIST <file_pdb_list>      : specifies file with the structure names'
      write(6,*)'                  xor'
      write(6,*)'         -TRAJ <file_pdb_traj>         : specifies the trajectory file with the structure pdb MODELs'
      write(6,*)'         [-FILE_PAR <file_parameters>] : specifies filename for parameter file. Default: bach++.par'
      write(6,*)'         [-FILE_PAR_AT <file_atomic_parameters>] : specifies filename for atomic parameter file. Default: ATOMIC_PARAMETERS_BSS'
      write(6,*)'         [-READ_MODE]                  : just reads the PDB file and outputs the lines considered as good'
      write(6,*)'         [-PREFIX <output_pdb_prefix>] : adds this prefix to output structures in READ_MODE option. Default: mod_'
      write(6,*)'         [-NO_SOLV]                    : speeds up calculations by considering only the pairwise potential'
      write(6,*)'         [-INTERFACE]                  : learns/computes only using interfaces of complexes. Disables -ONLY_GUT if present'
      write(6,*)'         [-STRICT_INTERFACE]           : learns/computes only using interfaces of complexes (strict def. of interface)'
      write(6,*)'         [-SEPARATE_TERMS]             : writes pairwise, solvation and clash contributions separately'
      write(6,*)'         [-o <output_file>]            : specifies output file name. Default: output.bv0'
      write(6,*)'         [-q]                          : quiet mode (writes output only on file, warnings suppressed)'
      write(6,*)'         [-qwarn]                      : quiet warnings mode (warnings suppressed)'
  endif
  j=0
  do i=1,argcount
      call get_command_argument(i,wq_char)
      if (trim(wq_char)=='-COMPUTE_PAR') then
          j=j+1
          COMPUTE_PARM=.true.
      endif
      if (trim(wq_char)=='-READ_MODE') then
          j=j+1
          READ_MODE=.true.
      endif
      if (trim(wq_char)=='-PREFIX') then
          call get_command_argument(i+1,prefix)
      endif
      if (trim(wq_char)=='-COMPUTE_ENE') then
          j=j+1
          COMPUTE_ENE=.true.
      endif
      if (trim(wq_char)=='-FILE_PAR') then
          call get_command_argument(i+1,filepar)
      endif
      if (trim(wq_char)=='-FILE_PAR_AT') then
          call get_command_argument(i+1,fileparat)
      endif
      if (trim(wq_char)=='-PDBLIST') then
          call get_command_argument(i+1,filepdb)
          PDBLIST=.true.
      endif
      if (trim(wq_char)=='-TRAJ') then
          call get_command_argument(i+1,filepdb)
          TRAJ=.true.
      endif
      if (trim(wq_char)=='-NO_SOLV')ONLY_GUT=.true.
      if (trim(wq_char)=='-q')QUIET=.true.
      if (trim(wq_char)=='-qwarn')QUIET_W=.true.
      if (trim(wq_char)=='-INTERFACE')then
          INTF=.true.
          ONLY_GUT=.false.
      endif
      if (trim(wq_char)=='-STRICT_INTERFACE') then
          S_INTF=.true.
          ONLY_GUT=.true.
      endif
      if (trim(wq_char)=='-o') then
          call get_command_argument(i+1,outputfile)
      endif
      if (trim(wq_char)=='-SEPARATE_TERMS')SEPARATE_TERMS=.true.
  enddo
  
  if(QUIET.eqv..true.) QUIET_W=.true.

  if(j==0)stop 'choose either -COMPUTE_PAR or -COMPUTE_ENE'
  if(j>1)stop 'sure you know what you are doing?'
  if(PDBLIST.and.TRAJ) stop 'choose either -TRAJ <file_pdb_traj> or -PDBLIST <file_pdb_list>'
  if(trim(filepdb)=="XXX")stop 'provide a file with a list of pdb structure files'


  open(unit=23,file=trim(fileparat),status='old')
  read(23,*)
  read(23,*)
  do i=1,10
      read(23,*) wq_char,PARAM_ALL(:,i)
  enddo
  read(23,*)
  do ii=1,20
      do j=1,nheavy_of_res(i_toname(ii))
          read(23,*) wq_char,order(ii,j),atomcode(ii,j),sasaparam(ii,j,:)
          apolar(ii,j+4)=fill_apolar(ii,j)
      enddo
  enddo
  do i=1,6
      read(23,*) wq_char,order(21,i),atomcode(21,i),sasaparam(21,i,:)
      if (i<5) then
          do j=1,20
              apolar(j,i)=fill_apolar(21,i)
          enddo
      endif
  enddo
  read(23,*)
  ii_old=0
  j_old=0
  do
      read(23,*,end=800) ii,j,exception
      if(ii==ii_old.or.j==j_old) then
          nexc(ii,j)=nexc(ii,j)+1
      else
          nexc(ii,j)=1
          ii_old=ii
          j_old=j
      endif
      exceptions(ii,j,nexc(ii,j))=exception
  enddo
800   continue
  close(23)
  order_p(:,:)=order(:,:)


  if (COMPUTE_ENE) then
      open(unit=22,file=trim(filepar),status='old')
      read(22,*)
      do i=1,20
          do j=1,20
              read(22,*) wq_char,wq_char,vGUT(i,j,:)
          enddo
      enddo
      read(22,*)
      do i=1,20
          read(22,*) wq_char,vSAS(i,:)
      enddo
      close(22)
  else if (COMPUTE_PARM) then
      vGUT(:,:,:)=0.d0
      vSAS(:,:)=0.d0
      nconct(:,:)=0.d0
      ntot(:)=0
      normaa=0.d0
      normt=0.d0
  endif

  if(.not.TRAJ)then
      open(22,file=filepdb,status='old')
      ii=0
      if (COMPUTE_ENE) then
          open(33,file=trim(outputfile),status='unknown')
          if (.not.SEPARATE_TERMS) then 
              write(33,'(44x,"PDB name",5x,"# res",12x,"TOT")')
          else
              write(33,'(44x,"PDB name",5x,"# res",7x,"PAIRWISE",10x,"CLASH",6x,"SOLVATION")')
          endif
      else if(READ_MODE) then
         open(33,file=trim(outputfile),status='unknown')
         write(33,'(44x,"PDB name",3x,"# progr",5x,"# res")')
      endif
      do 
          read(22,'(a150)',end=200) pdbname
          namelen=len(trim(pdbname))
          if (namelen>45) then
              pdbname45=pdbname(namelen-44:)
          else 
              pdbname45=pdbname(:45)
          endif
          ii=ii+1
          if (COMPUTE_PARM.and.(.not.QUIET)) write(6,'(i6,a100)')ii,trim(pdbname45)
          call read_pdb()
          if(READ_MODE) then
              write(33,'("E ",a50,2i10)')trim(pdbname45),ii,sum(numres(:))
              cycle
          endif
          if(.not.ONLY_GUT) then
              call compute_surf()
              call solvate()
          endif
          call compute_cmap()
          call compute_cmap_hb()
          call compute_GUTene()
          if (COMPUTE_ENE) then
              if((.not.S_INTF).and.SEPARATE_TERMS) then
                  if(.not.QUIET) write(6,'(3x,"ENERGIES",a50,i10,3F15.6)')trim(pdbname45),sum(numres(:)),&
                                        &0.6*sum(eneGUT(:)),0.6*0.03*sum(en_clashes(:)),sum(eneSOLV(:))
                  write(33,'("E ",a155,i10,3F15.6)')trim(pdbname),sum(numres(:)),0.6*sum(eneGUT(:)),0.6*0.03*sum(en_clashes(:)),sum(eneSOLV(:))
              else if (S_INTF.and.SEPARATE_TERMS) then
                  if(.not.QUIET) write(6,'(3x,"ENERGIES",a50,i10,2F15.6)')trim(pdbname45),sum(numres(:)),&
                                        &0.6*(sum(eneGUT_int(:))-eneGUT_int(3)),0.6*0.03*en_clashes(2)
                  write(33,'("E ",a155,i10,2F15.6)')trim(pdbname),sum(numres(:)),0.6*(sum(eneGUT_int(:))-eneGUT_int(3)),0.6*0.03*en_clashes(2)
              else if (S_INTF.and.(.not.SEPARATE_TERMS)) then
                  if(.not.QUIET) write(6,'(3x,"ENERGIES",a50,i10,F15.6)')trim(pdbname45),sum(numres(:)),&
                                        &0.6*(sum(eneGUT_int(:))-eneGUT_int(3))+0.6*0.03*en_clashes(2)
                  write(33,'("E",a155,i10,F15.6)')trim(pdbname),sum(numres(:)),0.6*(sum(eneGUT_int(:))-eneGUT_int(3))+0.6*0.03*en_clashes(2)
              else
                  if(.not.QUIET) write(6,'(3x,"ENERGIES",a50,i10,F15.6)')trim(pdbname45),sum(numres(:)),&
                                            &0.6*sum(eneGUT(:))+0.6*0.03*sum(en_clashes(:))+sum(eneSOLV(:))
                  write(33,'("E ",a155,i10,F15.6)')trim(pdbname),sum(numres(:)),0.6*sum(eneGUT(:))+0.6*0.03*sum(en_clashes(:))+sum(eneSOLV(:))
              endif
          endif
          deallocate(h_pos,atom_pos,ResName,thr_couple,clash_param)
          deallocate(ncmap,nbridgecmap)
          deallocate(resarea)
          deallocate(strictlist)
      enddo
200   continue
      close(22)
      if (COMPUTE_PARM) then
          open(unit=23,file=trim(filepar),status='unknown')
          write (23,*)'GUT POTENTIAL'
          do i=1,20
              do j=1,20
                  do ii=0,5
                      if(isnan(vGUT(i,j,ii))) vGUT(i,j,ii)=0.d0
                  enddo
                  write (23,'(2i3,6F12.6)') i,j,vGUT(i,j,:)
              enddo
          enddo
          write (23,*)'SOLVATION POTENTIAL'
          do i=1,20
              write(23,'(i8,2F16.8)')i,vSAS(i,1),vSAS(i,2)
          enddo
          close(23)   
      else
          close(33)
      endif
  else
      open(22,file=filepdb,status='old')
      open(33,file=trim(outputfile),status='unknown')  
      NM=0
      NL=0
      do
          NL=NL+1
          read(22,'(a100)',end=900)line
          if(NM>0) line_r(NL)=line
          read(line(1:5),'(a5)')check
          if(check=="MODEL") then
              if(NM>0) then
                  call read_pdb()
                  if(.not.READ_MODE) then
                      if(.not.ONLY_GUT) then
                          call compute_surf()
                          call solvate()
                      endif
                      call compute_cmap()
                      call compute_cmap_hb()
                      call compute_GUTene()
                      if(ONLY_GUT) then
                          if(.not.QUIET) write(6,'(3x,"ENERGIES ","MODEL",2i10,F15.6)')NM,sum(numres(:)),&
                                                &0.6*sum(eneGUT(:))+0.6*0.03*sum(en_clashes(:))+sum(eneSOLV(:))
                          write(33,'("E ","MODEL",2i10,F15.6)')NM,sum(numres(:)),0.6*sum(eneGUT(:))+0.6*0.03*sum(en_clashes(:))+sum(eneSOLV(:))
                      else 
                          if(.not.QUIET) write(6,'(3x,"ENERGIES ","MODEL",2i10,F15.6)')NM,sum(numres(:)),&
                                                &0.6*sum(eneGUT(:))+0.6*0.03*sum(en_clashes(:))+sum(eneSOLV(:))
                          write(33,'("E ","MODEL",2i10,F15.6)')NM,sum(numres(:)),0.6*sum(eneGUT(:))+0.6*0.03*sum(en_clashes(:))+sum(eneSOLV(:))
                      endif
                      deallocate(h_pos,atom_pos,ResName,thr_couple,clash_param)
                      deallocate(ncmap,nbridgecmap)
                      deallocate(resarea)
                      deallocate(strictlist)
                  endif
              endif
              NM=NM+1
              NL=0
          endif
      enddo
      close(22)
      close(33)
  endif
900  continue

contains

  subroutine compute_surf()
    implicit none

    integer a,aa,b,bb,c,cc,i,ii,j,jj,l,ll,l3,k,kk,m,mm,n,nn,p,pp,ic1,counter,ic
    real*8  dr(3),r,s2

    allocate(numresneigh(nrestot),nneigh(14,nrestot),resneigh(nrestot,nrestot),nhor_cp(nrestot),ap_cp(3,14,nrestot),nti_cp(nrestot))
    allocate(sum_1(14,nrestot),sum_2(14,nrestot),sum_3(14,nrestot),area(14,nrestot),radius(14,nrestot))

    numresneigh(:) = 0
    resarea(:) = 0.d0
    sum_1(:,:) = 0.d0
    sum_2(:,:) = 0.d0
    sum_3(:,:) = 0.d0
    nneigh(:,:) = 0
    resneigh(:,:) = 0
    n = 0

    do ic=1,NCHAIN
        do i=1,NumRes(ic)
            n = n+1        
            nhor_cp(n) = nheavy_of_res(ResName(i,ic))
            nti_cp(n) = name_toi(ResName(i,ic))
            do k=1,4
                radius(k,n) = radparam(order(21,k))
                ap_cp(1:3,k,n) = atom_pos(1:3,k,i,ic)
            enddo
            if(nti_cp(n)==11) cycle
            do k=5,4+nhor_cp(n)
                radius(k,n) = radparam(order(nti_cp(n),k-4))
                ap_cp(1:3,k,n) = atom_pos(1:3,k,i,ic)
            enddo
        enddo
    enddo

  !Filling of a matrix of neighboring RESIDUES (neighbor # of residue #)

    do i=1,nrestot
        do j=i,nrestot !counts itself as a neighbor in order to count atoms of the same residue
            dr(1:3)=ap_cp(1:3,1,i)-ap_cp(1:3,1,j)
            r=sqrt(sum(dr(1:3)**2)) 
            if (r > RESDIST) cycle 
            numresneigh(i) = numresneigh(i)+1
            resneigh(numresneigh(i),i) = j
            if (i/=j) then
                numresneigh(j) = numresneigh(j)+1
                resneigh(numresneigh(j),j) = i
            endif
        enddo
    enddo
  
    NXNEIGH_max = MAXVAL(numresneigh(:))*14
    allocate(neighborhood(2,NXNEIGH_max,14,nrestot))
    neighborhood(:,:,:,:) = 0

  !Filling of a matrix of neighboring ATOMS (neighbor <m',l'> of atom <m,l>)

    do l=1,nrestot !on the residues
        do m=1,4+nhor_cp(l)
            do ll=1,numresneigh(l) !sui residui vicini (ATT: NON E' IN NOTAZIONE GIUSTA)
                l3=resneigh(ll,l)
                do mm = 1,4+nhor_cp(l3)
                    if (l==l3.and.m==mm) cycle
                    if (radius(m,l)+radius(mm,l3) > distance(m,mm,l,l3)) then
                        nneigh(m,l) = nneigh(m,l)+1
                        neighborhood(1,nneigh(m,l),m,l) = mm
                        neighborhood(2,nneigh(m,l),m,l) = resneigh(ll,l)
                    endif
                enddo
            enddo
        enddo
    enddo

  !Sums computation

    do l=1,nrestot
        do m=1,4+nhor_cp(l)
            do n=1,nneigh(m,l)
                s2=0.d0
                a = neighborhood(1,n,m,l)
                b = neighborhood(2,n,m,l)
                sum_1(m,l) = sum_1(m,l) + Aij(m,a,l,b,radius(m,l),radius(a,b))
                do k=1,nneigh(a,b)
                    c = neighborhood(1,k,a,b)
                    ll = neighborhood(2,k,a,b)
                    if (ll==l.and.c==m) cycle
                    if (radius(m,l)+radius(c,ll)<distance(m,c,l,ll)) cycle
                    s2 = s2 + Aij(a,c,b,ll,radius(a,b),radius(c,ll))
                end do
                sum_2(m,l) = sum_2(m,l) + s2
                sum_3(m,l) = sum_3(m,l) + s2*Aij(m,a,l,b,radius(m,l),radius(a,b)) 
            end do
        end do
    end do

  !Residues' SASA computation

    do l=1,nrestot
        do m=1,4
            param(:)=sasaparam(21,5,:)
            area(m,l) = param(1)*4*pi*(radius(m,l)**2) + param(2)*sum_1(m,l) + param(3)*sum_2(m,l) + param(4)*sum_3(m,l)
            resarea(l) = resarea(l) + area(m,l)
        enddo
        if(nti_cp(l)/=11) then
            do m=5,4+nhor_cp(l)
                param(:)=sasaparam(21,5,:)
                area(m,l) = param(1)*4*pi*(radius(m,l)**2) + param(2)*sum_1(m,l) + param(3)*sum_2(m,l) + param(4)*sum_3(m,l)
                resarea(l) = resarea(l) + area(m,l)
            end do
        endif
    end do

    deallocate(neighborhood)
    deallocate(numresneigh,nneigh,resneigh,nhor_cp,ap_cp)
    deallocate(sum_1,sum_2,sum_3,area,radius)

    return
  end subroutine compute_surf

  
!SOLVATION ENERGY  
  
  subroutine solvate()
    implicit none

    ! intent in/out
    
    integer ii,nn
    integer nr2(20)
    real*8 coeff, func

    eneSOLV(:)=0.d0
    do i=1,nrestot
        ii = nti_cp(i)
        if(COMPUTE_ENE) then
            if(resarea(i)>threshold) then
                eneSOLV(1)=eneSOLV(1)+vSAS(ii,1)
            else
                eneSOLV(2)=eneSOLV(2)+vSAS(ii,2)
            endif
        else
            if(resarea(i)>threshold) then
                nconct(ii,1) = nconct(ii,1)+1
                normaa = normaa+1
            else
                nconct(ii,2) = nconct(ii,2)+1
            endif
            ntot(ii)=ntot(ii)+1
            normt = normt+1
        endif
    enddo

    deallocate(nti_cp)

    if(COMPUTE_PARM)then
       do i=1,20
          vSAS(i,1)=-log(dble(nconct(i,1)))+log(dble(ntot(i)))+log(normaa)-log(normt)
          vSAS(i,2)=-log(dble(nconct(i,2)))+log(dble(ntot(i)))+log(normt-normaa)-log(normt)
       enddo
    endif

  return
  end subroutine solvate


  subroutine compute_cmap() !computes a matrix of neighbor chains
    implicit none

    integer i,j,m,n,m_max,n_max,ic1,ic2,ii,jj,nmin_1,nmin_2
    real*8 pos(3,50,NRES_max),r,dr(3),x(3),rmin,norm,clash_res,rmin_ini,atomcode_copy(20,14)
    character*4 a_1,a_2

    ncmap(:,:,:,:)=0
    en_clashes(:)=0.d0
    nmin_1=1
    nmin_2=1
    rmin_ini=R_interaction
    
    do i=1,20
        do j=1,4+nheavy_of_res(i_toname(i))
            if (j<5) then
                atomcode_copy(i,j)=atomcode(21,j)
            else
                atomcode_copy(i,j)=atomcode(i,j-4)
            endif
        enddo
    enddo

    do ic1=1,NCHAIN
        do ic2=ic1,NCHAIN  
            do i=1,NumRes(ic1)
                ii=name_toi(ResName(i,ic1))
                m_max=nheavy_of_res(ResName(i,ic1))
                do j=1,NumRes(ic2)
                    if(ic1==ic2.and.j<=i+CBIAS)cycle
                    dr(1:3)=atom_pos(1:3,2,i,ic1)-atom_pos(1:3,2,j,ic2)
                    if (sqrt(sum(dr(1:3)**2))>VERLET_RADIUS) cycle
                    if (sqrt(sum(dr(1:3)**2))<INTERFACE_THRESHOLD) then
                        strictlist(i,ic1)=ic2
                        strictlist(j,ic2)=ic1
                    endif
                    jj=name_toi(ResName(j,ic2))
                    n_max=nheavy_of_res(ResName(j,ic2))
                    clash_res=0.d0
                    rmin=rmin_ini
                    do m=1,4+m_max
                        do n=1,4+n_max
                            dr(1:3)=atom_pos(1:3,m,i,ic1)-atom_pos(1:3,n,j,ic2)
                            r=sqrt(sum(dr(1:3)**2))
                            if (r>CLASH_THR) cycle
                            clash_res = clash_res+clash_function(r,atomcode_copy(ii,m)*atomcode_copy(jj,n))
                            if((ic1==ic2.and.j<=i+DBIAS).or.m<5.or.n<5) cycle
                            norm=thr_couple(atomcode_copy(ii,m)*atomcode_copy(jj,n))
                            if(r/norm<rmin.and.r<R_interaction) then
                                rmin=r/norm
                                nmin_1=m
                                nmin_2=n
                                ncmap(i,j,ic1,ic2)=1
                            endif
                        enddo
                    enddo
                    if(ic1==ic2) then
                        en_clashes(1) = en_clashes(1) + clash_res
                    else
                        en_clashes(2) = en_clashes(2) + clash_res
                    endif
                    ncmap(i,j,ic1,ic2)=ncmap(i,j,ic1,ic2)*(ncmap(i,j,ic1,ic2)+apolar(ii,nmin_1)*apolar(jj,nmin_2))
                enddo
            enddo
        enddo
    enddo

    return
  end subroutine compute_cmap


  subroutine compute_cmap_hb()

    implicit none

    integer i,j,ic1,ic2,sc
    real*8 ENEHBTH,dr(3),dist
    real*8 eneij,eneji,eneimjp,enejmip,eneimj,eneijp,enejmi,enejip

    nbridgecmap(:,:,:,:)=3
    ENEHBTH=-1.0d0
     
    do ic1=1,NCHAIN
        do i=1,NumRes(ic1)-1
            dr(:)=atom_pos(:,3,i,ic1)-atom_pos(:,4,i,ic1)
            dist=sqrt(sum(dr(:)**2))
            h_pos(:,i+1,ic1)=atom_pos(:,1,i+1,ic1)+(atom_pos(:,3,i,ic1)-atom_pos(:,4,i,ic1))/dist
        enddo
    enddo

    !       Computation of DSSP bridges; given the residue pair i,j
    !       I need to consider residues i-1,i+1,j-1,j+1 as well
    !       alpha-turns are possible ONLY if j-i=3

    do ic1=1,NCHAIN
        do ic2=ic1,NCHAIN
            if(ic1==ic2)then
                sc=1
            else
                sc=0
            endif
            do i=2,NumRes(ic1)-1
                do j=2+sc*(i+1),NumRes(ic2)-1
                    enejip=hbonds_geometry(i+1,j,j,i+1,ic1,ic2)
                    eneimj=hbonds_geometry(j,i-1,i-1,j,ic2,ic1)
                    enejmi=hbonds_geometry(i,j-1,j-1,i,ic1,ic2)
                    eneijp=hbonds_geometry(j+1,i,i,j+1,ic2,ic1)
                    eneji=hbonds_geometry(i,j,j,i,ic1,ic2)
                    eneij=hbonds_geometry(j,i,i,j,ic2,ic1)
                    enejmip=hbonds_geometry(i+1,j-1,j-1,i+1,ic1,ic2)
                    eneimjp=hbonds_geometry(j+1,i-1,i-1,j+1,ic2,ic1)
                    if((enejip.lt.ENEHBTH.and.eneimj.lt.ENEHBTH).or.(eneijp.lt.ENEHBTH.and.enejmi.lt.ENEHBTH))then
                        nbridgecmap(i,j,ic1,ic2)=0  ! PARALLEL BRIDGE
                    else if ((enejmip.lt.ENEHBTH.and.  eneimjp.lt.ENEHBTH).or.(eneij.lt.ENEHBTH .and.eneji.lt.ENEHBTH)) then
                        nbridgecmap(i,j,ic1,ic2)=1 ! ANTIPARALLEL BRIDGE
                    else if ((eneimj.lt.ENEHBTH.and.eneijp.lt.ENEHBTH.and.(j==i+3)).or.&
                      (enejmi.lt.ENEHBTH.and.enejip.lt.ENEHBTH .and.(i==j+3)))then
                         nbridgecmap(i,j,ic1,ic2)=2
                    endif
                enddo
            enddo
        enddo
    enddo

    return
  end subroutine compute_cmap_hb


  subroutine compute_GUTene() 
    implicit none

    integer i,j,l,ii,jj,ic1,ic2,h,mm
    real*8 add,prob_all,den_all,num_all(0:5),den_seq(20,20)
    integer, save ::  nGUT(20,20,0:5)=0

    if (COMPUTE_ENE) then
        eneGUT(:)=0.d0
        eneGUT_int(:)=0.d0
    endif

    if(NCHAIN==1.and.(INTF.or.S_INTF)) write(6,'(3x,a77)')"WARNING: cannot calculate energy at interface because there is just one chain"


    do ic1=1,NCHAIN
        do ic2=ic1,NCHAIN
            do i=1,NumRes(ic1)
                ii=name_toi(ResName(i,ic1))
                do j=1,NumRes(ic2)
                    if(ic1==ic2.and.j<=i+2.and.DBIAS>1)cycle
                    jj=name_toi(ResName(j,ic2))
                    l=nbridgecmap(i,j,ic1,ic2)+floor(nbridgecmap(i,j,ic1,ic2)/3.d0)*ncmap(i,j,ic1,ic2)
                    if (COMPUTE_ENE) then
                        !!!!!!!!!write(*,*) ic1, ic2, i, j, l, ResName(i,ic1), " ",  ResName(j,ic2), vGUT(ii,jj,l), nbridgecmap(i,j,ic1,ic2), ncmap(i,j,ic1,ic2)
                        eneGUT(l)=eneGUT(l)+vGUT(ii,jj,l)
                    else if (COMPUTE_PARM) then 
                        nGUT(ii,jj,l)=nGUT(ii,jj,l)+1
                        if (ii.ne.jj) nGUT(jj,ii,l)=nGUT(jj,ii,l)+1
                    endif
                enddo 
            enddo
        enddo
    enddo

    if(NCHAIN>1.and.COMPUTE_ENE)then
        do ic1=1,NCHAIN
            do i=1,NumRes(ic1)
                do ic2=ic1+1,NCHAIN
                    do j=1,NumRes(ic2)
                        if(strictlist(i,ic1)/=ic2) cycle
                        l=nbridgecmap(i,j,ic1,ic2)+floor(nbridgecmap(i,j,ic1,ic2)/3.d0)*ncmap(i,j,ic1,ic2)
                        ii=name_toi(ResName(i,ic1))
                        jj=name_toi(ResName(j,ic2))
                        eneGUT_int(l) = eneGUT_int(l)+vGUT(ii,jj,l)
                    enddo
                enddo
            enddo
        enddo
    else if (INTF.eqv..true..and.NCHAIN>1) then
        do ic1=1,NCHAIN
            do i=1,NumRes(ic1)
                if(resarea(i)<=threshold) cycle
                do ic2=ic1+1,NCHAIN
                    do j=1,NumRes(ic2)
                        if(strictlist(j,ic2)/=ic1.and.resarea(j)<=threshold) cycle
                        l=nbridgecmap(i,j,ic1,ic2)+floor(nbridgecmap(i,j,ic1,ic2)/3.d0)*ncmap(i,j,ic1,ic2)
                        ii=name_toi(ResName(i,ic1))
                        jj=name_toi(ResName(j,ic2))
                        if(COMPUTE_ENE)then
                            eneGUT_int(l) = eneGUT_int(l)+vGUT(ii,jj,l)
                        else if (COMPUTE_PARM) then
                            nGUT(ii,jj,l)=nGUT(ii,jj,l)+1
                            if (ii.ne.jj) nGUT(jj,ii,l)=nGUT(jj,ii,l)+1
                        endif
                    enddo
                enddo
            enddo
        enddo
    endif
                    

    if (COMPUTE_PARM) then
        den_all=0.d0
        num_all(:)=0.d0
        den_seq(:,:)=0.d0
        do l=0,5
            do ii=1,20
                do jj=ii,20
                    den_all=den_all+nGUT(ii,jj,l)
                    num_all(l)=num_all(l)+nGUT(ii,jj,l)
                    den_seq(ii,jj)=den_seq(ii,jj)+nGUT(ii,jj,l)
                    if(jj/=ii) den_seq(jj,ii)=den_seq(jj,ii)+nGUT(ii,jj,l)
                    vGUT(ii,jj,l)=nGUT(ii,jj,l)
                    if(jj/=ii) vGUT(jj,ii,l)=nGUT(jj,ii,l)
                enddo
            enddo
        enddo

        do l=0,5
            prob_all=num_all(l)/den_all
            if (isnan(prob_all)) prob_all=0.d0
            do ii=1,20
                do jj=ii,20
                    if(isnan((vGUT(ii,jj,l)/den_seq(ii,jj))/prob_all)) then
                        vGUT(ii,jj,l)=0.d0
                        if(jj/=ii) vGUT(jj,ii,l)=0.d0
                    else
                        vGUT(ii,jj,l)=-log((vGUT(ii,jj,l)/den_seq(ii,jj))/prob_all)
                        if(jj/=ii) vGUT(jj,ii,l)=-log((vGUT(jj,ii,l)/den_seq(jj,ii))/prob_all)
                    endif
                enddo
            enddo
        enddo
    endif

    return
  end subroutine compute_GUTene


  subroutine read_pdb()
    implicit none

    integer                   :: goodline(NLINE_max), usedj(14)
    integer                   :: i,is,ic,ii,is_old,nres,c,n,jj,nl_old,nn
    integer                   :: iH,iheavy,nheavy,nh,nalt,skipup,nlast
    integer                   :: struct_holes,mainrec
    integer, allocatable      :: ind(:), invind(:), writeorder(:)
    real*8                    :: x(3),dr_CA(3),dr_CN1(3),dr_CN2(3),dist_CA,dist_CN
    real*8, allocatable       :: occ(:)
    character*1               :: chainID,altLoc,chainID_old
    character*1               :: ins_code,ins_code_old
    character*3               :: ResName_tmp,ResName_tmp_old
    character*4               :: atom_name
    character*11              :: altmode, altmode_old
    character*3               :: name_prov
    character*16              :: piece1
    character*83              :: piece2
    character*4, allocatable  :: altseq1(:)
    character*10, allocatable :: altseq2(:)
    logical                   :: warn=.false., check2=.false., ini=.true., found


    !! 0) Reads the file, counts the lines, checks for altLocs
    if(.not.TRAJ)then
        open(50,file=pdbname,status='old')
        NL=0
        do while (NL<NLINE_max)
            NL=NL+1
            read(50,'(a100)',end=100)line_r(NL)
            line=line_r(NL)
            if(index(line(1:4),"ATOM").eq.1.and.index(line(17:17)," ").ne.1) warn=.true.
        enddo
        write(6,'(3x,a50,i10)')"FATAL: number of PDB file lines exceeds NLINE_max=",NLINE_max
        return
        close(50)
    endif
100 continue
    NL = NL-1

    goodline(:)=0
    n=0

    if(warn.eqv..true.) then
        allocate(altseq1(NL),altseq2(NL),ind(NL),invind(NL),occ(NL))
        invind(:)=0
    endif
    !! 0)

    !! 1) Discards "bad lines": checks for copies of atoms with altLoc=/" ". If there is an instance with altLoc==" " it keeps that, otherwise
    !!    it keeps the instance with higher occupancy (if present, otherwise it keeps the first met)
    do i=1,NL
        line=line_r(i)
        if(index(line(1:4),"ATOM").ne.1) then
            goodline(i)=-1
            cycle
        endif
        read(line(17:17),'(a1)')altLoc
        if(altLoc/=" ")then
            n = n+1
            read(line(13:16),'(a4)')altseq1(n)
            read(line(18:27),'(a10)')altseq2(n)
            read(line(57:60),'(f5.2)')occ(n)
            ind(n)=i
            invind(i)=n
        endif
    enddo
    if(n>0)then
        j=1
        i=0
        do while(j<=n)
            i = i+1
            if(goodline(i)<0) cycle
            line=line_r(i)
            if(invind(i)==j) then
                check2=.true.
                j = j+1
                cycle
            endif
            if(index(line(13:16),altseq1(j)).eq.1.and.index(line(18:27),altseq2(j)).eq.1) then
                goodline(ind(j))=-1
                i = ind(j)
                j = j+1
            endif
        end do
        if(check2.eqv..true.)then
            do j=1,n-1
                do jj=j+1,n
                    if(altseq1(j)==altseq1(jj).and.altseq2(j)==altseq2(jj)) then
                        if(occ(j)>=occ(jj))then
                            goodline(ind(jj))=-1
                        else
                            goodline(ind(j))=-1
                        endif
                    endif
                enddo
            enddo
        endif
    endif
    !! 1)

    if(warn.eqv..true.) deallocate(altseq1,altseq2,ind,invind,occ)
    NumRes(:)=0
    holes(:)=0
    nn=0

    !! 2) Checks if the residues have the proper atoms. If not, it returns warnings and omits the residue (if there are less atoms) or the extra atoms
    ini=.true.
    mainrec=0
    allocate(writeorder(NL))
    writeorder(:)=0
    do i=1,NL
        if(goodline(i)<0) cycle
        line=line_r(i)
        read(line(18:20),*) ResName_tmp
        if(name_toi(ResName_tmp)==-999) then
            write(6,'(3x,a26,a3)')"FATAL: resname not found: ",ResName_tmp
            return
        endif
        read(line(13:16),'(a4)')atom_name
        atom_name=trim(adjustl(atom_name))
        if(atom_name=="OXT".or.len_trim(atom_name)==0)then
            cycle
        endif
        iH=index(atom_name,"H")
        iheavy=max(index(atom_name,"C"),index(atom_name,"N"),index(atom_name,"S"),index(atom_name,"O"))
        if(.not.(iH==0.or.(iheavy>0.and.iH>iheavy))) cycle
        read(line(23:26),*)is !residue index in PDB
        read(line(27:27),'(a1)')ins_code !insertion code (if present)
        read(line(22:22),'(a1)')chainID !chain index in PDB
        read(line(31:54),'(3F8.3)')x(:)
        if(.not.ini) then
            if(is/=is_old.or.ins_code/=ins_code_old.or.chainID/=chainID_old) then
                if(c/=4+nheavy_of_res(ResName_tmp_old))then
                    if (.not.QUIET_W) write(6,'(3x,a40,a3,i5,a,a9,i5,/,12x,a28)')"WARNING: too few heavy atoms in residue ",ResName_tmp_old,&
                                  &is_old,ins_code_old," at line ",nl_old,"This residue will be omitted"
                    do j=1,4+nheavy_of_res(ResName_tmp_old)
                        if(usedj(j)/=0) goodline(usedj(j))=0
                    enddo
                    NumRes(ic)=NumRes(ic)-1
                    holes(ic)=holes(ic)+1
                else
                    mainrec = mainrec + 4+nheavy_of_res(ResName_tmp_old)
                endif
                NumRes(ic)=NumRes(ic)+1
                if(chainID/=chainID_old) then
                    if(holes(ic)/=0.and.(.not.QUIET_W)) write(6,'(3x,a9,i5,a38,a,a)')"WARNING: ",holes(ic),&
                                                        &" holes in PDB residue index of chain '",chainID_old,"'"
                    chainID_old=chainID
                    ic=ic+1
                else
                    if(is-is_old>1) holes(ic)=holes(ic)+is-is_old-1
                endif
                is_old=is
                ins_code_old=ins_code
                ResName_tmp_old=ResName_tmp
                nl_old=i
                c=0
                usedj(:)=0
            else
                if(ResName_tmp_old/=ResName_tmp) then
                    goodline(i)=0
                    cycle
                endif
            endif
        else
            ic=1
            NumRes(ic)=0
            nl_old=i
            c=0
            usedj(:)=0
            is_old=is
            ins_code_old=ins_code
            chainID_old=chainID
            ResName_tmp_old=ResName_tmp
            ini=.false.
        endif
        found=.false.
        do j=1,4
            if(atom_name==order_p(21,j)) then
                if(usedj(j)/=0) then
                    if (.not.QUIET_W) write(6,'(3x,a36,a4,a12,a3,a9,i5,/,12x,a25)')"WARNING: multiple instances of atom ",atom_name,&
                                      &" in residue ",ResName_tmp," at line ",i,"This atom will be omitted"
                    found=.true.
                    cycle
                endif
                usedj(j)=i
                c=c+1
                goodline(i)=1
                found=.true.
                writeorder(mainrec+j)=i
                exit
            endif
        enddo
        if(.not.found) then
            do j=1,nheavy_of_res(ResName_tmp)
                if(atom_name==order_p(name_toi(ResName_tmp),j)) then
                    if(usedj(4+j)/=0) then
                        if (.not.QUIET_W) write(6,'(3x,a36,a4,a12,a3,a9,i5,/,12x,a25)')"WARNING: multiple instances of atom ",atom_name,&
                                          &" in residue ",ResName_tmp," at line ",i,"This atom will be omitted"
                        found=.true.
                        cycle
                    endif
                    usedj(4+j)=i
                    c=c+1
                    goodline(i)=1
                    found=.true.
                    writeorder(mainrec+j+4)=i
                    exit
                endif
            enddo
        endif
        if(.not.found) then
            do j=1,nheavy_of_res(ResName_tmp)
                do n=1,nexc(name_toi(ResName_tmp),j)
                    if(atom_name==exceptions(name_toi(ResName_tmp),j,n)) then
                        if(usedj(4+j)/=0) then
                            if (.not.QUIET_W) write(6,'(3x,a36,a4,a12,a3,a9,i5,/,12x,a25)')"WARNING: multiple instances of atom ",atom_name,&
                                              &" in residue ",ResName_tmp," at line ",i,"This atom will be omitted"
                            found=.true.
                            cycle
                        endif
                        name_prov = exceptions(name_toi(ResName_tmp),j,n)
                        exceptions(name_toi(ResName_tmp),j,n) = order_p(name_toi(ResName_tmp),j)
                        order_p(name_toi(ResName_tmp),j) = name_prov
                        usedj(4+j)=i
                        c=c+1
                        goodline(i)=1
                        found=.true.
                        writeorder(mainrec+j+4)=i
                        exit
                    endif
                enddo
            enddo
        endif
        if(.not.found) then
            do j=1,4
                do n=1,nexc(21,j)
                    if(atom_name==exceptions(21,j,n)) then
                        if(usedj(j)/=0) then
                            if (.not.QUIET_W) write(6,'(3x,a36,a4,a12,a3,a9,i5,/,12x,a25)')"WARNING: multiple instances of atom ",atom_name,&
                                              &" in residue ",ResName_tmp," at line ",i,"This atom will be omitted"
                            found=.true.
                            cycle
                        endif
                        name_prov = exceptions(21,j,n)
                        exceptions(21,j,n) = order_p(21,j)
                        order_p(21,j) = name_prov
                        usedj(j)=i
                        c=c+1
                        goodline(i)=1
                        found=.true.
                        writeorder(mainrec+j)=i
                        exit
                    endif
                enddo
            enddo
        endif
        if((.not.found).and.((index(trim(atom_name),"O")==1.and.index(trim(atom_name),"1")==len_trim(atom_name)).or.index(trim(atom_name),"T")==len_trim(atom_name))) then
            atom_name="O"
            do j=1,4
                if(atom_name==order_p(21,j)) then
                    if(usedj(j)/=0) then
                        if (.not.QUIET_W) write(6,'(3x,a36,a4,a12,a3,a9,i5,a30,/,12x,a25)')"WARNING: multiple instances of atom ",atom_name,&
                                          &" in residue ",ResName_tmp," at line ",i,trim(pdbname45),"This atom will be omitted"
                        found=.true.
                        cycle
                    endif
                    usedj(j)=i
                    c=c+1
                    goodline(i)=1
                    found=.true.
                    writeorder(mainrec+j)=i
                    exit
                endif
            enddo
        endif
        if((.not.found).and.((index(trim(atom_name),"O")==1.and.index(trim(atom_name),"2")==len_trim(atom_name)).or.index(trim(atom_name),"T")==len_trim(atom_name))) then
            cycle
        endif
        if((.not.found).and.(.not.QUIET_W)) write(6,'(3x,a19,a4,a12,a3,a9,i5,a22,a30/,12x,a25)')"WARNING: atom name ",atom_name,&
                                          &" in residue ",ResName_tmp," at line ",i," not found in database",trim(pdbname45),"This atom will be omitted"
    enddo
    if(c/=4+nheavy_of_res(ResName_tmp_old))then
        if (.not.QUIET_W) write(6,'(3x,a40,a3,i5,a,a9,i5,/,12x,a28,i5)')"WARNING: too few heavy atoms in residue ",ResName_tmp_old,&
                                &is_old,ins_code_old," at line ",nl_old,"This residue will be omitted", c
        do j=1,4+nheavy_of_res(ResName_tmp_old)
            if(usedj(j)/=0) goodline(usedj(j))=0
            writeorder(mainrec+j)=0
        enddo
    else
        NumRes(ic)=NumRes(ic)+1
        mainrec=mainrec+4+nheavy_of_res(ResName_tmp_old)
    endif
    if(holes(ic)/=0.and.(.not.QUIET_W)) write(6,'(3x,a9,i5,a38,a,a)')"WARNING: ",holes(ic),&
                                        &" holes in PDB residue index of chain '",chainID_old,"'"
    !! 2) end

    
    nrestot=sum(NumRes(:))
    NCHAIN=ic 
    NXCHAIN_max=NCHAIN
    NXRES_max=MAXVAL(NumRes(:))
   
    if(READ_MODE) then
        open(34,file=trim(prefix)//trim(pdbname),status='unknown')
        do i=1,mainrec
            if(goodline(writeorder(i))<1) cycle
            line=line_r(writeorder(i))
            read(line(1:16),'(a16)')piece1
            read(line(18:100),'(a83)')piece2
            write(34,'(a16,x,a83)')piece1,piece2
        enddo
        close(34)
        return
    endif

    deallocate(writeorder)

    allocate(h_pos(3,NXRES_max,NXCHAIN_max),atom_pos(3,20,NXRES_max,NXCHAIN_max),ResName(NXRES_max,NXCHAIN_max))
    allocate(ncmap(NXRES_max,NXRES_max,NXCHAIN_max,NXCHAIN_max),nbridgecmap(NXRES_max,NXRES_max,NXCHAIN_max,NXCHAIN_max))
    allocate(resarea(nrestot))
    allocate(strictlist(NXRES_max,NXCHAIN_max))
    strictlist(:,:)=0

    nres=0
    ic=1
    is_old = -999
    n=0


    !! 3) fills arrays: atom and residue names, atom positions. In the 'atom' vectors the atom order is sorted with the help of the ATOMIC_PARAMETERS file
    do i=1,NL
        if(goodline(i)<1) cycle
        line=line_r(i)
        read(line(23:26),*)is
        read(line(27:27),'(a1)')ins_code
        if((is.ne.is_old).or.(ins_code_old.ne.ins_code)) then
            nheavy=0
            nres = nres+1
            if(nres>NumRes(ic)) then
                ic = ic+1
                nres = 1
            endif
        endif
        read(line(18:20),*)ResName(nres,ic)
        ii = name_toi(ResName(nres,ic))
        is_old=is
        ins_code_old=ins_code  
        read(line(13:16),*)atom_name
        read(line(31:54),'(3F8.3)')x(:)
        if(goodline(i)==1) then
            do j=1,4
                if(trim(adjustl(atom_name))==order_p(21,j)) then
                    atom_pos(:,j,nres,ic)=x(:)
                endif
            enddo
            do j=1,nheavy_of_res(ResName(nres,ic))
                if(trim(adjustl(atom_name))==order_p(ii,j)) then
                    atom_pos(:,j+4,nres,ic)=x(:)
                endif
            enddo
        endif
    enddo
    !! 3)


    !! 4) fills sparse arrays (atom codes)
    MAX_CODE=MAXVAL(PARAM_ALL(4,:))
    allocate(thr_couple(MAX_CODE),clash_param(MAX_CODE,2))
    thr_couple(:)=0.d0
    clash_param(:,:)=0.d0
    do i=1,SIZE(PARAM_ALL(4,:))
        thr_couple(PARAM_ALL(4,i))=PARAM_ALL(3,i)
        clash_param(PARAM_ALL(4,i),1)=PARAM_ALL(1,i)
        clash_param(PARAM_ALL(4,i),2)=PARAM_ALL(2,i)
    enddo
    !! 4)       

    close(50)
    return
  end subroutine read_pdb


  real*8 function hbonds_geometry(i_n,i_c,i_o,i_h,ic_nh,ic_co)

    implicit none
    integer i_n,i_c,i_o,i_h,ic_co,ic_nh
    real*8 dr(3)
    real*8 distcn,distch,diston,distoh

    if(ResName(i_n,ic_nh)=="PRO") then
        hbonds_geometry=0.d0
        return
    endif

    dr(:)=atom_pos(:,3,i_c,ic_co)-atom_pos(:,1,i_n,ic_nh)
    distcn=sqrt(sum(dr(:)**2))

    dr(:)=atom_pos(:,3,i_c,ic_co)-h_pos(:,i_h,ic_nh)
    distch=sqrt(sum(dr(:)**2))

    dr(:)=atom_pos(:,4,i_o,ic_co)-atom_pos(:,1,i_n,ic_nh)
    diston=sqrt(sum(dr(:)**2))

    dr(:)=atom_pos(:,4,i_o,ic_co)-h_pos(:,i_h,ic_nh)
    distoh=sqrt(sum(dr(:)**2))

    hbonds_geometry=(1/distch+1/diston-1/distoh-1/distcn)*332.d0*0.42d0*0.20d0
!    write(*,*) "HBOND", hbonds_geometry, i_n, ic_nh, ResName(i_n,ic_nh), i_c, ic_co, ResName(i_c,ic_co)
!    write(*,*) distch, diston, distoh, distcn
!    write(*,*) "N", atom_pos(:,1,i_n,ic_nh)
!    write(*,*) "H", h_pos(:,i_h,ic_nh)
!    write(*,*) "C", atom_pos(:,3,i_c,ic_co)
!    write(*,*) "O", atom_pos(:,4,i_o,ic_co)
    return
  end function hbonds_geometry


  integer function fill_apolar(i,m)
    implicit none

    integer i,m

    if(index(order(i,m),"C")==1.or.index(order(i,m),"S")==1) then
        fill_apolar=1
    else if(index(order(i,m),"N")==1.or.index(order(i,m),"O")==1) then
        fill_apolar=0
    else
        write(6,'(3x,a44)')"FATAL: something is wrong with atom names..."
    endif

    if(i==4.and.m==5) fill_apolar=1
    if(i==8) fill_apolar=0
    if(i==13.and.m==2) fill_apolar=0
    if(i==16.and.m==5) fill_apolar=0
    if(i==17.and.m==3) fill_apolar=0
    if(i==18.and.m==2) fill_apolar=0
    if(i==20.and.m==3) fill_apolar=0

  end function fill_apolar


  integer function nheavy_of_res(rname)
    implicit none

    character*3 rname
    nheavy_of_res=-999
    if(rname=="ALA")nheavy_of_res=1
    if(rname=="ARG")nheavy_of_res=7
    if(rname=="ASN")nheavy_of_res=4
    if(rname=="ASP")nheavy_of_res=4
    if(rname=="CYS")nheavy_of_res=2
    if(rname=="GLN")nheavy_of_res=5
    if(rname=="GLU")nheavy_of_res=5
    if(rname=="GLY")nheavy_of_res=0
    if(rname=="HIS".or. rname=="HIE")nheavy_of_res=6
    if(rname=="ILE")nheavy_of_res=4
    if(rname=="LEU")nheavy_of_res=4
    if(rname=="LYS")nheavy_of_res=5
    if(rname=="MET")nheavy_of_res=4
    if(rname=="PHE")nheavy_of_res=7
    if(rname=="PRO")nheavy_of_res=3
    if(rname=="SER")nheavy_of_res=2
    if(rname=="THR")nheavy_of_res=3
    if(rname=="TRP")nheavy_of_res=10
    if(rname=="TYR")nheavy_of_res=8
    if(rname=="VAL")nheavy_of_res=3
  end function nheavy_of_res

  integer function name_toi(res)
    implicit none
    integer iaa
    character res*3
    iaa=-999
    if(res=="CYS")iaa=1
    if(res=="PHE")iaa=2
    if(res=="LEU")iaa=3
    if(res=="TRP")iaa=4
    if(res=="VAL")iaa=5
    if(res=="ILE")iaa=6
    if(res=="MET")iaa=7
    if(res=="HIS".or.res=="HIE")iaa=8
    if(res=="TYR")iaa=9
    if(res=="ALA")iaa=10
    if(res=="GLY")iaa=11
    if(res=="PRO")iaa=12
    if(res=="ASN")iaa=13
    if(res=="THR")iaa=14
    if(res=="SER")iaa=15
    if(res=="ARG")iaa=16
    if(res=="GLN")iaa=17
    if(res=="ASP")iaa=18
    if(res=="LYS")iaa=19
    if(res=="GLU")iaa=20
    name_toi=iaa

    return
  end function name_toi

  character*3 function i_toname(iaa)
    implicit none
    integer iaa
    character*3 res
    if(iaa==1)res="CYS"
    if(iaa==2)res="PHE"
    if(iaa==3)res="LEU"
    if(iaa==4)res="TRP"
    if(iaa==5)res="VAL"
    if(iaa==6)res="ILE"
    if(iaa==7)res="MET"
    if(iaa==8)res="HIS"
    if(iaa==9)res="TYR"
    if(iaa==10)res="ALA"
    if(iaa==11)res="GLY"
    if(iaa==12)res="PRO"
    if(iaa==13)res="ASN"
    if(iaa==14)res="THR"
    if(iaa==15)res="SER"
    if(iaa==16)res="ARG"
    if(iaa==17)res="GLN"
    if(iaa==18)res="ASP"
    if(iaa==19)res="LYS"
    if(iaa==20)res="GLU"
    i_toname=res

    return
  end function i_toname

  real function distance(m1,m2,l1,l2)
    integer q
    integer m1,m2,l1,l2
    real*8 dr(3)
    
    dr(1:3)=ap_cp(1:3,m1,l1)-ap_cp(1:3,m2,l2)
    distance=sqrt(sum(dr(1:3)**2))     

  end function distance

  real function Aij(m1,m2,l1,l2,r_1,r_2)
        integer m1,m2,l1,l2
        real*8 r_1,r_2,d
        character*6 k1,k2
        
        d = distance(m1,m2,l1,l2)
        
        Aij = 2*pi*r_1*(r_1 - d/2 - (r_1**2 - r_2**2)/(2*d))
        
  end function Aij

  real function radparam(a)
    implicit none
    character*3 a

    radparam = 3.07d0
        
    return
  end function radparam

  real function norm_radii(a)
    implicit none
    character*3 a

    norm_radii = 0.d0
    if(index(a,"C")==1) norm_radii=1.70d0+R_PROBE
    if(index(a,"N")==1) norm_radii=1.65d0+R_PROBE
    if(index(a,"O")==1) norm_radii=1.60d0+R_PROBE
    if(index(a,"S")==1) norm_radii=1.90d0+R_PROBE
    if(norm_radii==0.d0) write(*,*)"ERROR"

    return
  end function norm_radii

        
  real function clash_function(x,code)
    implicit none
    real*8 x,a,b,code

    a=clash_param(code,1)
    b=clash_param(code,2)
    if(x<-b/(2*a)) then
        clash_function = a*x**2.d0 + b*x + (b**2.d0)/(4.d0*a)
    else
        clash_function = 0.d0
    end if

  return
  end function clash_function

end program tt
