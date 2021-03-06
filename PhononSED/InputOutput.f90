! -------------------------------------------------------------------------
! File I/O routines
!-------------------------------------------------------------------------
! Copyright (c) 2018 Gaurav Kumar, Daniel C. Elton
!
! This software is licensed under The MIT License (MIT)
! Permission is hereby granted, free of charge, to any person obtaining a copy of this
! software and associated documentation files (the "Software"), to deal in the Software
! without restriction, including without limitation the rights to use, copy, modify, merge,
! publish, distribute, sublicense, and/or sell copies of the Software, and to permit
! persons to whom the Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
! BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
! DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!-------------------------------------------------------------------------------------
module InputOutput
 use main_vars
 use lun_management
 implicit None

 contains

!------------------------------------------------------------
!----------- Read input file parameters --------------------
!------------------------------------------------------------
subroutine read_input_files
 use dans_timer
 implicit none
 integer :: luninp, Nlines, EoF
 real(8) :: a, c, a1, a2, a3, denom
 real(8) :: MC = 12.011000, MN = 14.007200, MO = 15.999430, MH = 1.0080000, MSi=28.0855
 real(8) :: MMg = 24.305, MHf = 178.49, MS = 32.065
 real(8)  :: n(3)


 call io_assign(luninp)
 open(luninp, file="input/PhononSED.inp", status='old', action='read')
!luninp = 5 !read from terminal standard in
 read(luninp, *) model
 read(luninp, *) a1
 read(luninp, *) a2
 read(luninp, *) a3
 read(luninp, *) fileheader
 read(luninp, *) Nunitcells
 read(luninp, *) Nk
 read(luninp, *) Neig
 read(luninp, *) fvel
 read(luninp, *) feig
 read(luninp, *) C_TYPE_EIGENVECTOR
 read(luninp, *) fcoord
 read(luninp, *) SUPERCELL_EIGENVECTOR
 read(luninp, *) Ntimesteps
 read(luninp, *) timestep
 read(luninp, *) READALL
 read(luninp, *) NPointsOut
 read(luninp, *) BTEMD
 read(luninp, *) GULPINPUT
 read(luninp, *) Ncorrptsout

 !call io_close(luninp)

 !-------------------- RDX  ------------------------------------------------
 if (model == 'RDX') then
    write(*,*) "Model is RDX"
    AtomsPerMolecule = 21
    MoleculesPerUnitCell = 8
    AtomsPerUnitCell = AtomsPerMolecule*MoleculesPerUnitCell !168
    Natoms = Nunitcells*AtomsPerUnitCell

    ! Assume equal number of PUC tiling in each lattice vector direction
    Nx = int(Nunitcells**(1./3.))
    Ny = int(Nunitcells**(1./3.))
    Nz = int(Nunitcells**(1./3.))

    allocate(MassPrefac(Natoms))
    allocate(freqs(Nk, Neig))
    allocate(eig_vecs(Nk, Neig, Natoms, 3))

    !build array of masses for ALL atoms
    idx = 1
    do i = 1, Nunitcells
        do j = 1, MoleculesPerUnitCell
            MassPrefac(idx+0:idx+2)   = MC ! Carbon
            MassPrefac(idx+3:idx+8)   = MN ! Nitrogen
            MassPrefac(idx+9:idx+14)  = MO ! Oxygen
            MassPrefac(idx+15:idx+20) = MH ! Hydrogen
            idx = idx + AtomsPerMolecule
        enddo
    enddo

    MassPrefac = sqrt(MassPrefac/real(Nunitcells))

    lattice_vector(1,:) = (/ a1, 0.d0, 0.d0 /)
    lattice_vector(2,:) = (/ 0.d0, a2, 0.d0 /)
    lattice_vector(3,:) = (/ 0.d0, 0.d0, a3 /)

     ! Tile lattice vectors such that they are commensurate with the SC used in the MD simulations (this will have knock on effects with regard to the 
     ! reciprocal lattice vectors and wherever they are used)
     lattice_vector(1,:) = Nx*lattice_vector(1,:)
     lattice_vector(2,:) = Ny*lattice_vector(2,:)
     lattice_vector(3,:) = Nz*lattice_vector(3,:)

     denom = DOT_PRODUCT(lattice_vector(1,:),CrossProd(lattice_vector(2,:),lattice_vector(3,:)))

     recip_lat_vec(1,:) = CrossProd(lattice_vector(2,:),lattice_vector(3,:))
     recip_lat_vec(2,:) = CrossProd(lattice_vector(3,:),lattice_vector(1,:))
     recip_lat_vec(3,:) = CrossProd(lattice_vector(1,:),lattice_vector(2,:))
                             
     recip_lat_vec = TwoPi*recip_lat_vec/denom

!-------------------- Magnesium Oxide (MgO) -------------------------------
else if (model == 'MgO') then
    write(*,*) "Model is MgO"
    AtomsPerMolecule = 8
    MoleculesPerUnitCell = 1
    AtomsPerUnitCell = AtomsPerMolecule*MoleculesPerUnitCell !168
    Natoms = Nunitcells*AtomsPerUnitCell

    ! Assume equal number of PUC tiling in each lattice vector direction
    Nx = int(Nunitcells**(1./3.))
    Ny = int(Nunitcells**(1./3.))
    Nz = int(Nunitcells**(1./3.))

    allocate(MassPrefac(Natoms))
    allocate(freqs(Nk, Neig))
    allocate(eig_vecs(Nk, Neig, Natoms, 3))

    !build array of masses for ALL atoms
    !idx = 1
    !do i = 1, Nunitcells
    !    do j = 1, MoleculesPerUnitCell
    !        MassPrefac(idx+0:idx+3)   = MMg ! Carbon
    !        MassPrefac(idx+4:idx+7)   = MO ! Nitrogen
    !        idx = idx + AtomsPerMolecule
    !    enddo
    !enddo
    MassPrefac(1:int(Natoms/2)) = MMg
    MassPrefac(int(Natoms/2)+1:Natoms) = MO
    MassPrefac = sqrt(MassPrefac/real(Nunitcells))

    lattice_vector(1,:) = (/ a1, 0.d0, 0.d0 /)
    lattice_vector(2,:) = (/ 0.d0, a2, 0.d0 /)
    lattice_vector(3,:) = (/ 0.d0, 0.d0, a3 /)

     ! Tile lattice vectors such that they are commensurate with the SC used in the MD simulations (this will have knock on effects with regard to the 
     ! reciprocal lattice vectors and wherever they are used)
     lattice_vector(1,:) = Nx*lattice_vector(1,:)
     lattice_vector(2,:) = Ny*lattice_vector(2,:)
     lattice_vector(3,:) = Nz*lattice_vector(3,:)

     denom = DOT_PRODUCT(lattice_vector(1,:),CrossProd(lattice_vector(2,:),lattice_vector(3,:)))

     recip_lat_vec(1,:) = CrossProd(lattice_vector(2,:),lattice_vector(3,:))
     recip_lat_vec(2,:) = CrossProd(lattice_vector(3,:),lattice_vector(1,:))
     recip_lat_vec(3,:) = CrossProd(lattice_vector(1,:),lattice_vector(2,:))
                             
     recip_lat_vec = TwoPi*recip_lat_vec/denom

!-------------------- Silicon (Si) -------------------------------
 else if (model .eq. 'silicon') then
     write(*,*) "Model is silicon"
     AtomsPerMolecule = 2
     MoleculesPerUnitCell = 1
     AtomsPerUnitCell = AtomsPerMolecule*MoleculesPerUnitCell !2
     Natoms = Nunitcells*AtomsPerUnitCell !2

     ! Assume equal number of PUC tiling in each lattice vector direction
     Nx = int(Nunitcells**(1./3.))
     Ny = int(Nunitcells**(1./3.))
     Nz = int(Nunitcells**(1./3.))

     allocate(MassPrefac(Natoms))
     allocate(freqs(Nk, Neig))
     allocate(eig_vecs(Nk, Neig, Natoms, 3))

     !build array of masses for ALL atoms
     MassPrefac = MSi
     MassPrefac = sqrt(MassPrefac/real(Nunitcells))
     
     lattice_vector(1,:) = 0.5*(/ 0.d0, a2, a3 /)
     lattice_vector(2,:) = 0.5*(/ a1, 0.d0, a3 /)
     lattice_vector(3,:) = 0.5*(/ a1, a2, 0.d0 /)

     ! Tile lattice vectors such that they are commensurate with the SC used in the MD simulations (this will have knock on effects with regard to the
     ! reciprocal lattice vectors and wherever they are used)
     lattice_vector(1,:) = Nx*lattice_vector(1,:)
     lattice_vector(2,:) = Ny*lattice_vector(2,:)
     lattice_vector(3,:) = Nz*lattice_vector(3,:)

     denom = DOT_PRODUCT(lattice_vector(1,:),CrossProd(lattice_vector(2,:),lattice_vector(3,:)))

     recip_lat_vec(1,:) = CrossProd(lattice_vector(2,:),lattice_vector(3,:))
     recip_lat_vec(2,:) = CrossProd(lattice_vector(3,:),lattice_vector(1,:))
     recip_lat_vec(3,:) = CrossProd(lattice_vector(1,:),lattice_vector(2,:))

     recip_lat_vec = TwoPi*recip_lat_vec/denom

 !-------------------- Halfnium disulphide (Si) -------------------------------
else if (model .eq. 'HfS2') then
        write(*,*) "Model is HfS2 (halfnium disulphide)"
        AtomsPerMolecule = 3
        MoleculesPerUnitCell = 1
        AtomsPerUnitCell = AtomsPerMolecule*MoleculesPerUnitCell !3
        Natoms = Nunitcells*AtomsPerUnitCell

        ! Assume equal number of PUC tiling in each lattice vector direction
        Nx = int(Nunitcells**(1./3.))
        Ny = int(Nunitcells**(1./3.))
        Nz = int(Nunitcells**(1./3.))

        allocate(MassPrefac(Natoms))
        allocate(freqs(Nk, Neig))
        allocate(eig_vecs(Nk, Neig, Natoms, 3))

        !build array of masses for ALL atoms
        MassPrefac(1:int(Natoms/3)) = MHf
        MassPrefac(int(Natoms/3)+1:Natoms) = MS
        MassPrefac = sqrt(MassPrefac/real(Nunitcells))

        !! Hexagonal closed packed cell parameters must be entered by hand here!!
        a = 3.6529
        c = 5.6544

        lattice_vector(1,:) = (/ a, 0.d0, 0.d0/)
        lattice_vector(2,:) = (/ -0.5*a, dsqrt(3.d0)/(2.d0*a), 0.d0 /)
        lattice_vector(3,:) = (/ 0.d0, 0.d0, c /)

        ! Tile lattice vectors such that they are commensurate with the SC used in the MD simulations (this will have knock on effects with regard to the
        ! reciprocal lattice vectors and wherever they are used)
        lattice_vector(1,:) = Nx*lattice_vector(1,:)
        lattice_vector(2,:) = Ny*lattice_vector(2,:)
        lattice_vector(3,:) = Nz*lattice_vector(3,:)

        denom = DOT_PRODUCT(lattice_vector(1,:),CrossProd(lattice_vector(2,:),lattice_vector(3,:)))

        recip_lat_vec(1,:) = CrossProd(lattice_vector(2,:),lattice_vector(3,:))
        recip_lat_vec(2,:) = CrossProd(lattice_vector(3,:),lattice_vector(1,:))
        recip_lat_vec(3,:) = CrossProd(lattice_vector(1,:),lattice_vector(2,:))

        recip_lat_vec = TwoPi*recip_lat_vec/denom

    else
    write(*,*) "ERROR : model not found!!"
 endif

 !------------- read length of velocities file ------------------------------
 if (READALL) then
    write(*,*) "Reading size of file ..."

    Nlines = 0

    call io_assign(lun)
    open(lun, file=fvel, status='old', action='read')

    do
        read(lun, *, iostat=EoF)
        if (.not.(EoF .eq. 0)) then
            exit
        else
            Nlines = Nlines + 1
        endif
    enddo
    Ntimesteps = floor(real(Nlines)/real(Natoms+9)) - 1
    write(*, *) "There are ", Ntimesteps, " timesteps in the file"

    call io_close(lun)

endif
!------------- allocations ----------- ------------------------------
 allocate(qdot(Ntimesteps))
 allocate(k_vectors(Nk, 3))
 allocate(r_eq(Natoms, 3))
 allocate(r(Natoms, 3))


 !read in eigenvectors we will be working with
 if (pid .eq. 0)  call start_timer("reading files")

 if (pid .eq. 0)  write(*, *) "reading eigenvector file... "

 call read_eigvector_file

 !read in necessary velocities & coordinate data
 if (pid .eq. 0)  write(*, *) "reading velocities file... "

 if (GULPINPUT) then
     call read_GULP_trajectory_file
  else
     call read_LAAMPS_files
 endif


 !!! PBC correction

 ! Vectors for PBC mapping
 P0(1,:) = (/ 0.d0, 0.d0, 0.d0 /)
 P0(2,:) = P0(1,:)
 P0(3,:) = P0(2,:)
 P0(4,:) = lattice_vector(1,:) + lattice_vector(2,:) + lattice_vector(3,:)
 P0(5,:) = P0(4,:)
 P0(6,:) = P0(5,:)

 normal(1,:) = CrossProd(lattice_vector(3,:),lattice_vector(2,:))
 normal(2,:) = CrossProd(lattice_vector(1,:),lattice_vector(3,:))
 normal(3,:) = CrossProd(lattice_vector(2,:),lattice_vector(1,:))
 normal(4,:) = -normal(1,:)
 normal(5,:) = -normal(2,:)
 normal(6,:) = -normal(3,:)

 boxmap(1,:) = lattice_vector(1,:)
 boxmap(2,:) = lattice_vector(2,:)
 boxmap(3,:) = lattice_vector(3,:)
 boxmap(4,:) = -lattice_vector(1,:)
 boxmap(5,:) = -lattice_vector(2,:)
 boxmap(6,:) = -lattice_vector(3,:)

 do ia = 1,Natoms
    do i = 1, 6
       r_eq(ia,:) = MapPBC(P0(i,:),normal(i,:),r_eq(ia,:),boxmap(i,:))
    end do
 end do

 if (C_TYPE_EIGENVECTOR) then
     write(*,*) "using equilibrium coords from GULP eig file , for use with C-type eigenvector representation ..."
     r = r_eq
     !do i = 1, Natoms
     !     write(*,*) r(i, :)
     !enddo
     !--- deprecated option to read equlibrium coordinates from an external file
     !write(*, *) "Since there is no GULP trajectory input, attempting to read equilibrium coords from external file..."
     !call io_assign(lun)
     !open(lun, file=fcoord, status='old', action='read')
     !r = one_frame_xyz(lun)
     !call io_close(lun)
 else
     if (pid .eq. 0) write(*, *) "generating unit cell coordinates for use with D-type eigenvector representation ..."
     !---- figure out equilibrium unit cell coordinates
     !---- this calculates the coordinate for the corner of the unit cell each atom is in
     !----
     !---- Here we work under the assumption that our simulation cell is the unit cell and therefore the lattice point is at 0.0
     do ia = 1, Natoms
         do ix = 1, 3
             r(ia, ix) =  floor(r_eq(ia,ix)/lattice_vector(ix,ix))*lattice_vector(ix,ix) !SC coord
         enddo
     enddo
 endif

end subroutine read_input_files


!------------------------------------------------------------
!---------------- Read eigenvector file --------------------
!------------------------------------------------------------
subroutine read_eigvector_file()
 integer  :: Natoms_file, Neig_file, my_iostat
 real(8)  :: mag
 double precision, dimension(3) :: realpart, cmplxpart
 character(len=10) :: junk, space
 character(len=12) :: junk_new


 call io_assign(luneig)
 open(luneig, file=feig, status='old', action='read')

 read(luneig, *) Natoms_file

 if ((AtomsPerUnitCell .gt. Natoms_file)) then
     write(*,*) "ERROR: N_atoms in eigenvector file ( ", Natoms_file," ) is", &
                 " less than the expected number of atoms per unit cell (", &
                 AtomsPerUnitCell, ")"
    stop
 endif
 if (SUPERCELL_EIGENVECTOR) then
     if ((Natoms .gt. Natoms_file)) then
         write(*, *) "WARNING.. I am looking for eigenvectors for a supercell representation &
                      in the input .eig file. I don't see enough atoms. Continuing anyway.."
     endif
else
    if ((AtomsPerUnitCell .lt. Natoms_file)) then
        write(*,*) "WARNING: N_atoms in eigenvector file ( ", Natoms_file," ) is", &
                 " greater than the expected number of atoms per unit cell (", &
                 AtomsPerUnitCell, "). I assume you know what you are doing  &
                 and will continue to see what happens... You may want to make sure &
                 supercell representation is set to TRUE in the input file."
    endif
 endif

 write(*,*) "reading equilibrium coordinates from GULP eigenvector file..."
 do ia = 1, Natoms_file
     read(luneig, *) junk, (r_eq(ia, ix), ix = 1,3)
 enddo

 read(luneig, *) !no of k points
 read(luneig, *) Neig_file
 write(*,'(a,i4,a)') "File contains ", Neig_file, " eigenvectors per k-point"
 if (Neig .gt. Neig_file) then
     write(*,*) "ERROR: Neig specified is larger than the number of  &
                 eigenvectors in the input file."
    stop
 endif


 do ik = 1, Nk
     read(luneig, '(a,3f10.6)') junk_new, (k_vectors(ik, ix), ix = 1,3) 
     k_vectors(ik, 1) = k_vectors(ik, 1)*recip_lat_vec(1,1)
     k_vectors(ik, 2) = k_vectors(ik, 2)*recip_lat_vec(2,2)
     k_vectors(ik, 3) = k_vectors(ik, 3)*recip_lat_vec(3,3)

     do ie = 1, Neig
        read(luneig, *) !Mode    x
        read(luneig, *) freqs(ik, ie)
        mag = 0

        do ia = 1, Natoms_file

            if (ik .eq. 1) then
                read(luneig, *) realpart(1), realpart(2), realpart(3)
                do ix = 1, 3
                    eig_vecs(ik, ie, ia, ix) = dcmplx(realpart(ix), 0)
                enddo
            else
                read(luneig, *, iostat=my_iostat) realpart(1), realpart(2), realpart(3), cmplxpart(1), cmplxpart(2), cmplxpart(3)
                if (my_iostat /= 0) then
                    write(*,*) 'WARNING: missing complex part of eigenvector in k=', ik, ' ieig = ', ie
                    backspace(luneig)
                    my_iostat = 0
                    read(luneig, *, iostat=my_iostat) realpart(1), realpart(2), realpart(3), space, &
                                                      cmplxpart(1), cmplxpart(2), cmplxpart(3)
                    if (my_iostat /= 0) then
                        write(*,'(a,i3,a,i3)') 'WARNING: error reading cmplx eigenvector for k=', ik, ' mode # = ', ie
                        read(luneig, *, iostat=my_iostat) realpart(1), realpart(2), realpart(3)
                        cmplxpart(:) = 0
                    endif
                    my_iostat = 0
                endif

                do ix = 1, 3
                    eig_vecs(ik, ie, ia, ix) = dcmplx(realpart(ix), cmplxpart(ix))
                enddo
            endif

            mag = mag + real(eig_vecs(ik, ie, ia, 1)*conjg(eig_vecs(ik, ie, ia, 1)) + &
                             eig_vecs(ik, ie, ia, 2)*conjg(eig_vecs(ik, ie, ia, 2)) + &
                             eig_vecs(ik, ie, ia, 3)*conjg(eig_vecs(ik, ie, ia, 3)) )

        enddo !do ia = 1, Natoms_file

        mag = sqrt(mag)
        eig_vecs(ik, ie, 1:AtomsPerUnitCell, :) = eig_vecs(ik, ie, 1:AtomsPerUnitCell, :)/mag !make a unit vector (normalization)

        if (.not. SUPERCELL_EIGENVECTOR) then
            !copy eigenvectors from 0,0,0 unit cell to other unit cells
            j = AtomsPerUnitCell
            do i = 2, Nunitcells
                eig_vecs(ik, ie, (i-1)*j+1 : (i-1)*j+j, :) = eig_vecs(ik, ie, 1:j, :) !fill in rest
            enddo
            !do j = 2, AtomsPerUnitCell
            !    do i = 2, Nunitcells
            !        eig_vecs(ik, ie, (i-1)*j+1 : (i-1)*j+j, :) = eig_vecs(ik, ie, 1:j, :) !fill in rest
            !enddo
         endif
     enddo  !do ie = 1, Neig

     !do ia = 1, Natoms
     !     write(*,*)  eig_vecs(1, 1, ia, :)
     !enddo

     ! read any remaining eigenvectors ls
     do ie = Neig+1, Neig_file
        read(luneig, *) !Mode    x
        read(luneig, *) !freqs
        do ia = 1, Natoms_file
            read(luneig, *)
        enddo
     enddo

 enddo !ik = 1, Nk
end subroutine read_eigvector_file

!-----------------------------------------------------------------------
!----------------- Read in necessary LAMMPS files ---------------------
!-----------------------------------------------------------------------
subroutine read_LAAMPS_files
 implicit none
 integer :: lun

     !------------- read velocities file ------
     call io_assign(lun)
     open(lun, file=fvel, status='old', action='read')

     allocate(velocities(Ntimesteps, Natoms, 3))
     do t = 1, Ntimesteps
         if (mod(t,1000).eq.0) write(*,*) t, Ntimesteps
         velocities(t, :, :) = 1000.0d0*one_frame(lun) !factor of 1000 to make velocity unit Ang/ps
     enddo

     call io_close(lun)

     !---------- read in coordinates data for BTEMD --------
     if (BTEMD) then
         write(*, *) "reading coordinates file... "

         call io_assign(lun)
         open(lun, file=fcoords, status='old', action='read')

         allocate(coordinates(Ntimesteps, Natoms, 3))
         do t = 1, Ntimesteps
             if (mod(t,1000).eq.0) write(*,*) t, Ntimesteps
             coordinates(t, :, :) = one_frame(lun)
         enddo

         call io_close(lun)
     endif !BTEMD

end subroutine read_LAAMPS_files


!-----------------------------------------------------------------------
!----------------- Read in GULP trajectory file -----------------------
!-----------------------------------------------------------------------
subroutine read_GULP_trajectory_file
 implicit none
 integer :: lun
 character(20) :: junk

     call io_assign(lun)
     open(lun, file=fvel, status='old', action='read')

     allocate(velocities(Ntimesteps, Natoms, 3))

     read(lun, *)
     read(lun, *)
     do t = 1, Ntimesteps
         if (mod(t,1000).eq.0) write(*,*) "reading", t, Ntimesteps
         do i = 1, 3
            read(lun, *) junk
         enddo
         do ia = 1, Natoms
             read(lun, *) !coordinates
         enddo
         read(lun, *)
         do ia = 1, Natoms
             read(lun, '(3ES26.16E2)') (velocities(t, ia, ix), ix=1,3)
             !if (t .eq. 1) then
             !    write(*,*) velocities(t, ia, :)
             !endif
         enddo
         do ia = 1, 2*Natoms+2
             read(lun, *) !skip derivatives & site energies data
         enddo
     enddo
     call io_close(lun)

end subroutine read_GULP_trajectory_file


!-----------------------------------------------------------------------
!----------------- Read in .vel file from Quantum Espresso ------------
!-----------------------------------------------------------------------
subroutine read_quantum_espresso
 implicit none
 integer :: lun

      !------------- read velocities file ------
      call io_assign(lun)
      open(lun, file=fvel, status='old', action='read')

      allocate(velocities(Ntimesteps, Natoms, 3))
      do t = 1, Ntimesteps
          if (mod(t,1000).eq.0) write(*,*) "reading: ", t, "of", Ntimesteps
          velocities(t, :, :) = one_frame_xyz(lun, .false.)
      enddo

      call io_close(lun)

 end subroutine read_quantum_espresso


!-----------------------------------------------------------------------
!----------------- Read one frame of velocities file ------------------
!-----------------------------------------------------------------------
function one_frame(lun)
 implicit none
 real(8), dimension(Natoms, 3) :: one_frame
 integer, intent(in) :: lun
 integer ::  Natoms_file, unused


 read(lun,*) !ITEM: TIMESTEP
 read(lun,*) !timestep
 read(lun,*) !ITEM: NUMBER OF ATOMS
 read(lun,*) Natoms_file
 if (.not.(Natoms .eq. Natoms_file)) then
     write(*,*) "ERROR: N_atoms in velocity file ( ", Natoms_file," ) does", &
                 " not match expected total number of atoms (", &
                 Natoms, ")"
    stop
 endif
 read(lun,*) !ITEM: BOX BOUNDS
 read(lun,*) box(1, :) !bounding box
 read(lun,*) box(2, :)
 read(lun,*) box(3, :)
 read(lun,*) !comment line

 !read velocities
 do ia = 1, Natoms
    read(lun,*) unused,(one_frame(ia, ix), ix=1,3)
 enddo

end function one_frame

!-----------------------------------------------------------------------
!----------------- Read one frame of a standard .xyz file -------------
!-----------------------------------------------------------------------
function one_frame_xyz(lun, has_comment_line)
 implicit none
 real(8), dimension(Natoms, 3) :: one_frame_xyz
 integer, intent(in) :: lun
 logical, optional :: has_comment_line
 integer ::  Natoms_file
 character(2) :: junk

 if (.not. present(has_comment_line)) then
   has_comment_line = .true.
 endif

 read(lun,*) Natoms_file
 if (has_comment_line) read(lun,*) !comment line

 do ia = 1, Natoms
    read(lun,*) junk, (one_frame_xyz(ia, ix), ix=1,3)
 enddo

end function one_frame_xyz




!---------------------------------------------------------------------
!-----------------  Print SED and other outputs ---------------------
!---------------------------------------------------------------------
subroutine print_SED
 implicit none

call system('mkdir /home/gkumar1/Phonon_Lifetimes/PhononSED/'//trim(fileheader)//"_N"//trim(str(Nunitcells))//"_k"//trim(str(Nk)))

 do ik = 1, Nk

     call io_open(lunout, filename="/home/gkumar1/Phonon_Lifetimes/PhononSED/"//trim(fileheader)//"_N"//trim(str(Nunitcells))//&
&"_k"//trim(str(Nk))//"/"//trim(fileheader)//"_N"//trim(str(Nunitcells))//"_k"//trim(str(ik))//"_SED.dat")

     do i = 1, NPointsOut
         write(lunout, '(f12.4,1x)', advance='no') freqs_smoothed(i)

         do ie = 1, Neig-1
             write(lunout, '(e12.6,1x)', advance='no') all_SED_smoothed(ik, ie, i)
         enddo

         write(lunout, '(e12.6)', advance='yes') all_SED_smoothed(ik, Neig, i)
     enddo
     call io_close(lunout)

     call io_open(lunout, filename="/home/gkumar1/Phonon_Lifetimes/PhononSED/"//trim(fileheader)//"_N"//trim(str(Nunitcells))//&
&"_k"//trim(str(Nk))//"/"//trim(fileheader)//"_N"//trim(str(Nunitcells))//"_k"//trim(str(ik))//"_frequencies.dat")
     do ie = 1, Neig
             write(lunout, '(f12.5,1x)') freqs(ik, ie)
     enddo
     call io_close(lunout)


enddo !ie = 1, Neig

 call io_open(lunout, filename="/home/gkumar1/Phonon_Lifetimes/PhononSED/"//trim(fileheader)//"_N"//trim(str(Nunitcells))//&
&"_k"//trim(str(Nk))//"/"//trim(fileheader)//"_N"//trim(str(Nunitcells))//"_k"//trim(str(Nk))//"_kvectors.dat")
 do ik = 1, Nk
    write(lunout,'(f5.3,1x,f5.3,1x,f5.3)') (k_vectors(ik, ix), ix = 1,3)
 enddo
 call io_close(lunout)





end subroutine print_SED


!---------------------------------------------------------------------
!-----------------  Print SED ---------------------------------------
!---------------------------------------------------------------------
subroutine print_corr_fns
 implicit none

 do ik = 1, Nk

     call io_open(lunout, filename=trim(fileheader)//"_N"//trim(str(Nunitcells))//"_k"//trim(str(ik))//"_cor_fun.dat")

     do i = 1, Ncorrptsout
         write(lunout, '(f12.4,1x)', advance='no') i*timestep

         do ie = 1, Neig-1
             write(lunout, '(e12.6,1x)', advance='no') all_corr_fns(ik, ie, i)
         enddo

         write(lunout, '(e12.6)', advance='yes') all_corr_fns(ik, Neig, i)
     enddo
     call io_close(lunout)
enddo !ie = 1, Neig




end subroutine print_corr_fns




!---------------------------------------------------------------------
!-----------------Convert an integer to string ----------------------
!---------------------------------------------------------------------
character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

end module InputOutput
