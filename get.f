      subroutine getsites(system,sites,nsites,nl,ixd,iyd)
c----------------------------------------------------------
c     finds all undercoordinated C atoms and orientation  -
c     of the dangling bond                                -
c                                                         -
c                                                         -
c                      6  1  2                            -
c                       \ | /                             -
c                         C                               -
c                       / | \                             -
c                      5  4  3                            -
c                                                         -
c----------------------------------------------------------

      integer*1 system(800,800,199)
      integer sites(9999,5)
      logical up

      n_s = 0

c-----loop over all sites
      do k=1,nl
         do i=1,ixd
            do j=1,iyd

               if(system(i,j,k).eq.1 .or. system(i,j,k).eq.5
     &              .or. system(i,j,k).eq.6) then

c-----up or down

               if(system(i,md(j+2),k).ne.0) then
                  up = .true.
               else
                  up = .false.
               end if

               if(n_s.gt.9999) then
                  write(*,*)'error: increase size of sites array'
                  stop
               end if

c-----if up, check for nearest neighbour vacancies

               if(up) then

                  if(system(i,md(j+2),k).ne.2  
     &                 .and. system(mb(i-1),md(j-1),k).ne.2  .and. 
     &                 system(mb(i+1),md(j-1),k).ne.2) then
                     go to 100
                  else if(system(i,md(j+2),k).eq.2 .and.
     &                    system(mb(i-1),md(j-1),k).ne.2 .and.
     &                    system(mb(i+1),md(j-1),k).ne.2) then
                     n_s = n_s + 1
                     sites(n_s,1) = i
                     sites(n_s,2) = j
                     sites(n_s,3) = k
                     sites(n_s,4) = 1
                  else if(system(i,md(j+2),k).ne.2 .and.
     &                    system(mb(i-1),md(j-1),k).ne.2 .and.
     &                    system(mb(i+1),md(j-1),k).eq.2) then
                     n_s = n_s + 1
                     sites(n_s,1) = i
                     sites(n_s,2) = j
                     sites(n_s,3) = k
                     sites(n_s,4) = 3
                  else if(system(i,md(j+2),k).ne.2 .and.
     &                    system(mb(i-1),md(j-1),k).eq.2 .and.
     &                    system(mb(i+1),md(j-1),k).ne.2) then
                     n_s = n_s + 1
                     sites(n_s,1) = i
                     sites(n_s,2) = j
                     sites(n_s,3) = k
                     sites(n_s,4) = 5
                  else
                     write(*,*)'Error in site - undercoordinated atom u'
                     write(*,*)system(i,md(j+2),k)
                     write(*,*)system(mb(i+1),md(j-1),k)
                     write(*,*)system(mb(i-1),md(j-1),k)
                     write(*,*)i,j,k
                     stop
                  end if

               else

                  if(system(i,md(j-2),k).ne.2 .and. 
     &                 system(mb(i-1),md(j+1),k).ne.2 .and. 
     &                 system(mb(i+1),md(j+1),k).ne.2) then
                     go to 100
                  else if(system(i,md(j-2),k).eq.2 .and.
     &                    system(mb(i-1),md(j+1),k).ne.2 .and.
     &                    system(mb(i+1),md(j+1),k).ne.2) then
                     n_s = n_s + 1
                     sites(n_s,1) = i
                     sites(n_s,2) = j
                     sites(n_s,3) = k
                     sites(n_s,4) = 4
                  else if(system(i,md(j-2),k).ne.2 .and.
     &                    system(mb(i-1),md(j+1),k).ne.2 .and.
     &                    system(mb(i+1),md(j+1),k).eq.2) then
                     n_s = n_s + 1
                     sites(n_s,1) = i
                     sites(n_s,2) = j
                     sites(n_s,3) = k
                     sites(n_s,4) = 2
                 else if(system(i,md(j-2),k).ne.2 .and.
     &                    system(mb(i-1),md(j+1),k).eq.2 .and.
     &                    system(mb(i+1),md(j+1),k).ne.2) then
                     n_s = n_s + 1
                     sites(n_s,1) = i
                     sites(n_s,2) = j
                     sites(n_s,3) = k
                     sites(n_s,4) = 6
                  else
                     write(*,*)'Error in site - undercoordinated atom d'
                     write(*,*)system(i,md(j-2),k)
                     write(*,*)system(mb(i+1),md(j+1),k)
                     write(*,*)system(mb(i-1),md(j+1),k)
                     write(*,*)i,j,k
                     stop
                  end if

               end if
               end if
 100           continue
            end do
         end do
      end do

      nsites = n_s

      end 

      subroutine getssties(system,ssites,snsites,nl,ixd,iyd)
c----------------------------------------------------------------
c     finds all split interstitial atoms                        -
c----------------------------------------------------------------
      integer*1 system(800,800,199)
      integer ssites(9999,5)
      integer snsites

      n_s = 0

c-----loop over all sites
      do k=1,nl
         do i=1,ixd
            do j=1,iyd
                              
               if(n_s.gt.999) then
                  write(*,*)'error: increase size of ssites array'
                  stop
               end if

               if(system(i,j,k).eq.9) then
                  n_s = n_s + 1
                  ssites(n_s,1) = i
                  ssites(n_s,2) = j
                  ssites(n_s,3) = k
                  if(system(i,md(j+2),k).ne.0) then
                     ssites(n_s,4) = 1
                  else
                     ssites(n_s,4) = 0
                  end if               
               end if

            end do
         end do
      end do

      snsites = n_s

      end

      subroutine getisites(isystem,isites,insites,nl,ixd,iyd)
c----------------------------------------------------------------
c     finds all interstitial atoms in the interlayer plane      -
c     (spiros and grafted)                                      -
c     fourth index is the interstitial type and direction       -
c     spiro interstitial: 1-6                                   -
c     lower grafted interstitial: 11-16                         -
c     upper grafted interstitial: 21-26                         -
c----------------------------------------------------------------
      integer*1 isystem(800,800,199)
      integer isites(9999,5)

      n_s = 0

c-----loop over all sites
      do k=1,nl
         do i=1,ixd
            do j=1,iyd

               if(n_s.gt.999) then
                  write(*,*)'error: increase size of isites array'
                  stop
               end if

               if(isystem(i,j,k).ge.1 .and. isystem(i,j,k).le.6) then
                  n_s = n_s + 1
                  isites(n_s,1) = i
                  isites(n_s,2) = j
                  isites(n_s,3) = k
                  isites(n_s,4) = isystem(i,j,k)
               else if(isystem(i,j,k).ge.11) then
                  n_s = n_s + 1
                  isites(n_s,1) = i
                  isites(n_s,2) = j
                  isites(n_s,3) = k
                  isites(n_s,4) = isystem(i,j,k)
               end if
            end do
         end do
      end do

      insites = n_s

      end

      subroutine getens(system,isystem,sites,nsites,energies,etab,stab,
     &     ntab,nstep,menergies,mtab,emin,enpos)
c----------------------------------------------------------
c     assigns energy barriers for undercoordinated atoms  -
c     to move into vacant lattice positions based on the  -
c     local atomic environment                            -
c                                                         -
c     energies(nsite,1) - in plane barrier                -
c                                                         -     
c     if no interplanar bond:                             -
c     energies(nsite,2) - bond with atom in plane above   -
c     energies(nsite,3) - bond with atom in plane below   -
c                                                         -
c     if interplanar bond:                                -
c     energies(nsite,2) - break with atom in plane above  -
c     energies(nsite,3) - break with atom in plane below  -
c                                                         -
c     energies(nsite,4) - move atom to plane above        -
c     energies(nsite,5) - move atom to plane below        -
c                                                         -
c     we check for spiro and split interstitials in the   -
c     interstitial planes above and below the vacancy     -
c                                                         -
c     whichever match to the library table occurs first,  -
c     above or below the lattie plane, is used for the    -
c     barriers: therefore, place the closest I-V          - 
c     structures at the start of the ens.dat file         -
c----------------------------------------------------------         

      integer*1 system(800,800,199)
      integer*1 isystem(800,800,199)
      integer sites(9999,5),enpos(9999)
      real*8  energies(9999,5),menergies(9999,5)
      character*1 cmat(41,20),c2mat(41,20)
      logical sw1,emin

      integer*2 smat(21,30),umat(21,30),dmat(21,30),amat(21,30)
      integer*2 stab(21,30,12,499),admat(21,30), stab_i,stab_l
      real*8 etab(12,499),mtab(12,499)

      do ns=1,nsites
         

c--------initialise energies list
         do ne=1,5
            energies(ns,ne) = 9.0
            menergies(ns,ne) = 0.0
         end do

c-----in-plane barriers

c.....map local system structure into temporary arrays
         
         ii = sites(ns,1)
         jj = sites(ns,2)
         kk = sites(ns,3)
         ndir = sites(ns,4)

         do i=1,21
            do j=1,30
               smat(i,j) = system(mb(ii-11+i),md(jj-16+j),kk)
c..............interstitial layer above the plane
               umat(i,j) = isystem(mb(ii-11+i),md(jj-16+j),kk)
c..............interstitial layer below the plane
               dmat(i,j) = isystem(mb(ii-11+i),md(jj-16+j),ml(kk-1))

c-----------changing system lables into library format lables

c...........interlayer bound atoms

               if(smat(i,j).eq.7) then
                  smat(i,j) = 2
               end if
               
               if(smat(i,j).eq.5 .or. smat(i,j).eq.6) then
                  smat(i,j) = 4
               end if

               amat(i,j) = 0
               admat(i,j) = 0

c...........interstitials
               if(umat(i,j).eq.8) then
                  umat(i,j) = 0
                  amat(i,j) = -1
               else if(umat(i,j).gt.20) then
                  umat(i,j) = 0
               else if(umat(i,j).ge.11 .and. umat(i,j).le.16) then
                  umat(i,j) = umat(i,j) + 10
               else if(umat(i,j).ge.1 .and. umat(i,j).le.6) then
                  umat(i,j) = umat(i,j) + 10
               end if

c...........interstitials
               if(dmat(i,j).eq.8) then
                  dmat(i,j) = 0
                  admat(i,j) = -1
               else if(dmat(i,j).ge.11 .and. dmat(i,j).le.16) then
                  dmat(i,j) = 0
               else if(dmat(i,j).ge.1 .and. dmat(i,j).le.6) then
                  dmat(i,j) = dmat(i,j) + 10
               end if

            end do
         end do

c         write(*,*)'site: ',ns
c         write(*,*)ii,jj,kk
         
c         do j=24,1,-1
c            write(*,'(17i1)') ( smat(i,j) , i=1,17 )
c         end do

c         do j=24,1,-1
c            write(*,'(17i1)') ( stab(i,j,3,1) , i=1,17 )
c         end do





cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     need to make the program select the closest off-layer defect, not
c     just any layer that matches either umat of dmat
c     
c     at the moment, if it doesn't find a structure matching the adjacent 
c     layer with defect, but then matches the opposite adjacent layer 
c     with nothing in it
c
c     it needs to match the CLOSEST - being higher in the ens.dat file
c     doesn't matter - it merely finds the next closest (or nothing) 
c     further down the table
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC







c........check tables for structure
         
c--------check above

         nchk = 0

         do m=1,ntab
            do n=(ndir*2-1),(ndir*2)
               nchk = 0
               do i=1,21
                  do j=1,30

                     stab_l = stab(i,j,n,m)

c--------------------separate out lattice plane and interstitial plane
c                     in the library record
                     if(stab_l.gt.10) then
                        stab_i = stab_l
                        if(stab_i.eq.17) then
                           stab_l = 1
                        else if(stab_i.eq.18) then
                           stab_l = 1
                        else if(stab_i.ge.31) then
                           stab_l = 2
                           stab_i = stab_i - 20
                        else
                           stab_l = 1
                        end if
                     else if(stab_l.eq.3) then
                        stab_i = 17
                     else
                        stab_i = 0
                     end if

c--------------------if the lattice plane matches and either of the 
c                    up or down interstitial planes

c                    if there is a different interstitial above and 
c                    below the plane, the first library match is used
                     if((smat(i,j).eq.stab_l .or. stab_l.eq.3) .and. 
     &                    (umat(i,j).eq.stab_i .or. stab_i.eq.17 .or.
     &                    (stab_i.eq.18 .and. umat(i,j).ge.11 .and.
     &                    umat(i,j).le.16))) then                           
                        nchk = nchk + 1
                     else
c                        write(*,*)i,j,smat(i,j),stab_l,umat(i,j),stab_i
                        go to 101
                     end if
                  end do
               end do
               if(nchk.eq.630) go to 102
 101           continue
            end do
         end do

c         write(*,*)nchk

c--------check below

         nchk = 0

         do m=1,ntab
            do n=(ndir*2-1),(ndir*2)
               nchk = 0
               do i=1,21
                  do j=1,30

                     stab_l = stab(i,j,n,m)

c--------------------separate out lattice plane and interstitial plane
c                     in the library record
                     if(stab_l.gt.10) then
                        stab_i = stab_l
                        if(stab_i.eq.17) then
                           stab_l = 1
                        else if(stab_i.eq.18) then
                           stab_l = 1
                        else if(stab_i.ge.31) then
                           stab_l = 2
                           stab_i = stab_i - 20
                        else
                           stab_l = 1
                        end if
                     else if(stab_l.eq.3) then
                        stab_i = 17
                     else
                        stab_i = 0
                     end if

c--------------------if the lattice plane matches and either of the 
c                    up or down interstitial planes

c                    if there is a different interstitial above and 
c                    below the plane, the first library match is used
c                     write(*,*)i,j,n,m,dmat(i,j),stab_i
                     if((smat(i,j).eq.stab_l .or. stab_l.eq.3) .and. 
     &                    (dmat(i,j).eq.stab_i .or. stab_i.eq.17 .or.
     &                    (stab_i.eq.18 .and. dmat(i,j).ge.11 .and.
     &                    dmat(i,j).le.16))) then
                        nchk = nchk + 1
                     else
                        go to 113
                     end if
                  end do
               end do
               if(nchk.eq.630) go to 102
 113           continue
            end do
         end do

c         write(*,*)nchk

c........if no match - write unmatched structure and stop

         write(*,*)'unmatched structure'
         write(*,*)ii,jj,kk

         ic = 0
         jc = 0

         sw1 = .false.

c........writing unmatched structure to output

         write(*,*)'lattice plane'

         do j=1,30
            ic = 0
            if(smat(1,j).ge.1 .or. smat(2,j).ge.1 .or. smat(3,j).ge.1
     &           .or. smat(4,j).ge.1 .or. smat(5,j).ge.1) then
               go to 202
            else 
               go to 203
            end if
 202        jc = jc + 1
            do i=1,21
               ic = ic + 1
               if(smat(i,j).eq.0) then
                  cmat(ic,jc) = ' '
               else if(smat(i,j).eq.1) then
                  cmat(ic,jc) = 'o'
               else if(smat(i,j).eq.2) then
                  cmat(ic,jc) = '-'
               else if(smat(i,j).eq.3) then
                  cmat(ic,jc) = '*'
               else if(smat(i,j).eq.4) then
                  cmat(ic,jc) = '+'
               else if(smat(i,j).eq.9) then
                  cmat(ic,jc) = 's'
               end if
               if(ic.eq.41) go to 201
               ic = ic + 1
               cmat(ic,jc) = ' '
            end do
 201        continue
 203        continue
         end do

         do j=20,1,-1
            write(*,'(41a1)') ( cmat(i,j) , i=1,41 )
         end do

         write(*,*)'lower interstitial plane'

         ic = 0
         jc = 0

         do j=1,30
            ic = 0
            if(smat(1,j).ge.1 .or. smat(2,j).ge.1 .or. smat(3,j).ge.1
     &           .or. smat(4,j).ge.1 .or. smat(5,j).ge.1) then
               go to 302
            else 
               go to 303
            end if
 302        jc = jc + 1
            do i=1,21
               ic = ic + 1
               if(dmat(i,j).eq.0) then
                  if(admat(i,j).eq.-1) then
                     c2mat(ic,jc) = '.'
                  else
                     c2mat(ic,jc) = ' '
                  end if
               else if(dmat(i,j).eq.11) then
                  c2mat(ic,jc) = '1'
               else if(dmat(i,j).eq.12) then
                  c2mat(ic,jc) = '2'
               else if(dmat(i,j).eq.13) then
                  c2mat(ic,jc) = '3'
               else if(dmat(i,j).eq.14) then
                  c2mat(ic,jc) = '4'
               else if(dmat(i,j).eq.15) then
                  c2mat(ic,jc) = '5'
               else if(dmat(i,j).eq.16) then
                  c2mat(ic,jc) = '6'
               else if(dmat(i,j).eq.21) then
                  c2mat(ic,jc) = 'a'
               else if(dmat(i,j).eq.22) then
                  c2mat(ic,jc) = 'b'
               else if(dmat(i,j).eq.23) then
                  c2mat(ic,jc) = 'c'
               else if(dmat(i,j).eq.24) then
                  c2mat(ic,jc) = 'd'
               else if(dmat(i,j).eq.25) then
                  c2mat(ic,jc) = 'e'
               else if(dmat(i,j).eq.26) then
                  c2mat(ic,jc) = 'f'
               else if(dmat(i,j).eq.31) then
                  c2mat(ic,jc) = 'A'
               else if(dmat(i,j).eq.32) then
                  c2mat(ic,jc) = 'B'
               else if(dmat(i,j).eq.33) then
                  c2mat(ic,jc) = 'C'
               else if(dmat(i,j).eq.34) then
                  c2mat(ic,jc) = 'D'
               else if(dmat(i,j).eq.35) then
                  c2mat(ic,jc) = 'E'
               else if(dmat(i,j).eq.36) then
                  c2mat(ic,jc) = 'F'
               else
                  write(*,*)'error in inter plane: unrecognised atom'
                  stop
               end if
               if(ic.eq.41) go to 301
               ic = ic + 1
               c2mat(ic,jc) = ' '
            end do
 301        continue
 303        continue
         end do

         do j=20,1,-1
            write(*,'(41a1)') ( c2mat(i,j) , i=1,41 )
         end do

         write(*,*)'upper interstitial plane'

         ic = 0
         jc = 0

         do j=1,30
            ic = 0
            if(smat(1,j).ge.1 .or. smat(2,j).ge.1 .or. smat(3,j).ge.1
     &           .or. smat(4,j).ge.1 .or. smat(5,j).ge.1) then
               go to 402
            else 
               go to 403
            end if
 402        jc = jc + 1
            do i=1,21
               ic = ic + 1
               if(umat(i,j).eq.0) then
                  if(amat(i,j).eq.-1) then
                     c2mat(ic,jc) = '.'
                  else
                     c2mat(ic,jc) = ' '
                  end if
               else if(umat(i,j).eq.11) then
                  c2mat(ic,jc) = '1'
               else if(umat(i,j).eq.12) then
                  c2mat(ic,jc) = '2'
               else if(umat(i,j).eq.13) then
                  c2mat(ic,jc) = '3'
               else if(umat(i,j).eq.14) then
                  c2mat(ic,jc) = '4'
               else if(umat(i,j).eq.15) then
                  c2mat(ic,jc) = '5'
               else if(umat(i,j).eq.16) then
                  c2mat(ic,jc) = '6'
               else if(umat(i,j).eq.21) then
                  c2mat(ic,jc) = 'a'
               else if(umat(i,j).eq.22) then
                  c2mat(ic,jc) = 'b'
               else if(umat(i,j).eq.23) then
                  c2mat(ic,jc) = 'c'
               else if(umat(i,j).eq.24) then
                  c2mat(ic,jc) = 'd'
               else if(umat(i,j).eq.25) then
                  c2mat(ic,jc) = 'e'
               else if(umat(i,j).eq.26) then
                  c2mat(ic,jc) = 'f'
               else if(umat(i,j).eq.31) then
                  c2mat(ic,jc) = 'A'
               else if(umat(i,j).eq.32) then
                  c2mat(ic,jc) = 'B'
               else if(umat(i,j).eq.33) then
                  c2mat(ic,jc) = 'C'
               else if(umat(i,j).eq.34) then
                  c2mat(ic,jc) = 'D'
               else if(umat(i,j).eq.35) then
                  c2mat(ic,jc) = 'E'
               else if(umat(i,j).eq.36) then
                  c2mat(ic,jc) = 'F'

               else
                  write(*,*)'error in inter plane: unrecognised atom'
                  stop
               end if
               if(ic.eq.41) go to 401
               ic = ic + 1
               c2mat(ic,jc) = ' '
            end do
 401        continue
 403        continue
         end do

         do j=20,1,-1
            write(*,'(41a1)') ( c2mat(i,j) , i=1,41 )
         end do

         stop

 102     energies(ns,1) = etab(1,m)
         menergies(ns,1) = mtab(1,m)
         enpos(ns) = m

c        write(*,*)'found structure :',m,n

c--------------------------------------------------------------------         
c-------------check for vacancies above and below -------------------
c--------------------------------------------------------------------
c-------------to fill the interlayer vacancy motion barriers --------
c--------------------------------------------------------------------

c--------loop through all the different possible configurations

c--------first, if the central atom is unbonded to neighbouring layer

         if(system(ii,jj,kk).eq.1) then
            if(ndir.eq.1) then
c-----------check above
c..............DVIL1
               if(system(ii,jj,ml(kk+1)).eq.1 .and. 
     &              system(ii,md(jj-2),ml(kk+1)).eq.2) then

                  iv = ivcoord(system,ii,md(jj-2),ml(kk+1))
                  if(iv.eq.2) then
                     energies(ns,2) = etab(7,m)
                     menergies(ns,2) = mtab(7,m)
                  else
                     energies(ns,2) = etab(4,m)
                     menergies(ns,2) = mtab(4,m)
                  end if

c..............DVIL2
               else if(system(ii,jj,ml(kk+1)).eq.1 .and. 
     &                 (system(mb(ii-1),md(jj+1),ml(kk+1)).eq.2 .or.
     &                 system(mb(ii+1),md(jj+1),ml(kk+1)).eq.2)) then

                  if(system(mb(ii-1),md(jj+1),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii-1),md(jj+1),ml(kk+1))
                  else if(system(mb(ii+1),md(jj+1),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii+1),md(jj+1),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(8,m)
                     menergies(ns,2) = mtab(8,m)
                  else
                     energies(ns,2) = etab(5,m)
                     menergies(ns,2) = mtab(5,m)
                  end if

c..............DVIL3
               else if((system(ii,jj,ml(kk+1)).eq.2 .and. 
     &                 system(mb(ii),md(jj-2),ml(kk+1)).eq.1) .or. 
     &                 (system(ii,md(jj+2),ml(kk+1)).eq.1 .and.
     &                 system(ii,md(jj+4),ml(kk+1)).eq.2)) then

                  if(system(ii,jj,ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk+1))
                  else if(system(ii,md(jj+4),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj+4),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(9,m)
                     menergies(ns,2) = mtab(9,m)
                  else
                     energies(ns,2) = etab(6,m)
                     menergies(ns,2) = mtab(6,m)
                  end if

               end if

c-----------check below
               if(system(ii,jj,ml(kk-1)).eq.1 .and. 
     &              system(ii,md(jj-2),ml(kk-1)).eq.2) then

                  iv = ivcoord(system,ii,md(jj-2),ml(kk-1))
                  if(iv.eq.2) then
                     energies(ns,3) = etab(7,m)
                     menergies(ns,3) = mtab(7,m)
                  else
                     energies(ns,3) = etab(4,m)
                     menergies(ns,3) = mtab(4,m)
                  end if

c..............DVIL2
               else if(system(ii,jj,ml(kk-1)).eq.1 .and. 
     &                 (system(mb(ii-1),md(jj+1),ml(kk-1)).eq.2 .or.
     &                 system(mb(ii+1),md(jj+1),ml(kk-1)).eq.2)) then

                  if(system(mb(ii-1),md(jj+1),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii-1),md(jj+1),ml(kk-1))
                  else if(system(mb(ii+1),md(jj+1),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii+1),md(jj+1),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(8,m)
                     menergies(ns,3) = mtab(8,m)
                  else
                     energies(ns,3) = etab(5,m)
                     menergies(ns,3) = mtab(5,m)
                  end if

c..............DVIL3
               else if((system(ii,jj,ml(kk-1)).eq.2 .and. 
     &                 system(mb(ii),md(jj-2),ml(kk-1)).eq.1) .or. 
     &                 (system(ii,md(jj+2),ml(kk-1)).eq.1 .and.
     &                 system(ii,md(jj+4),ml(kk-1)).eq.2)) then

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk-1))
                  else if(system(ii,md(jj+4),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj+4),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(9,m)
                     menergies(ns,3) = mtab(9,m)
                  else
                     energies(ns,3) = etab(6,m)
                     menergies(ns,3) = mtab(6,m)
                  end if

               end if
               
            else if(ndir.eq.2) then
c-----------check above
c..............DVIL1
               if(system(ii,jj,ml(kk+1)).eq.1 .and. 
     &              system(mb(ii-1),md(jj-1),ml(kk+1)).eq.2) then

                  iv = ivcoord(system,md(ii-1),md(jj-1),ml(kk+1))
                  if(iv.eq.2) then
                     energies(ns,2) = etab(7,m)
                     menergies(ns,2) = mtab(7,m)
                  else
                     energies(ns,2) = etab(4,m)
                     menergies(ns,2) = mtab(4,m)
                  end if

c..............DVIL2
               else if(system(ii,jj,ml(kk+1)).eq.1 .and. 
     &                 (system(mb(ii+1),md(jj-1),ml(kk+1)).eq.2 .or.
     &                 system(ii,md(jj+2),ml(kk+1)).eq.2)) then

                  if(system(mb(ii+1),md(jj-1),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii+1),md(jj-1),ml(kk+1))
                  else if(system(ii,md(jj+2),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj+2),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(8,m)
                     menergies(ns,2) = mtab(8,m)
                  else
                     energies(ns,2) = etab(5,m)
                     menergies(ns,2) = mtab(5,m)
                  end if

c..............DVIL3
               else if((system(ii,jj,ml(kk+1)).eq.2 .and. 
     &                 system(mb(ii-1),md(jj-1),ml(kk+1)).eq.1) .or.
     &                 (system(mb(ii+1),md(jj+1),ml(kk+1)).eq.1 .and.
     &                 system(mb(ii+2),md(jj+2),ml(kk+1)).eq.2)) then

                  if(system(ii,jj,ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk+1))
                  else if(system(mb(ii+2),md(jj+2),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii+2),md(jj+2),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(9,m)
                     menergies(ns,2) = mtab(9,m)
                  else
                     energies(ns,2) = etab(6,m)
                     menergies(ns,2) = mtab(6,m)
                  end if

               end if
c-----------check below
c..............DVIL1
               if(system(ii,jj,ml(kk-1)).eq.1 .and. 
     &              system(mb(ii-1),md(jj-1),ml(kk-1)).eq.2) then

                  iv = ivcoord(system,md(ii-1),md(jj-1),ml(kk-1))
                  if(iv.eq.2) then
                     energies(ns,3) = etab(7,m)
                     menergies(ns,3) = mtab(7,m)
                  else
                     energies(ns,3) = etab(4,m)
                     menergies(ns,3) = mtab(4,m)
                  end if

c..............DVIL2
               else if(system(ii,jj,ml(kk-1)).eq.1 .and. 
     &                 (system(mb(ii+1),md(jj-1),ml(kk-1)).eq.2 .or.
     &                 system(ii,md(jj+2),ml(kk-1)).eq.2)) then

                  if(system(mb(ii+1),md(jj-1),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii+1),md(jj-1),ml(kk-1))
                  else if(system(ii,md(jj+2),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj+2),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(8,m)
                     menergies(ns,3) = mtab(8,m)
                  else
                     energies(ns,3) = etab(5,m)
                     menergies(ns,3) = mtab(5,m)
                  end if

c..............DVIL3
               else if((system(ii,jj,ml(kk-1)).eq.2 .and. 
     &                 system(mb(ii-1),md(jj-1),ml(kk-1)).eq.1) .or.
     &                 (system(mb(ii+1),md(jj+1),ml(kk-1)).eq.1 .and.
     &                 system(mb(ii+2),md(jj+2),ml(kk-1)).eq.2)) then

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk-1))
                  else if(system(mb(ii+2),md(jj+2),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii+2),md(jj+2),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(9,m)
                     menergies(ns,3) = mtab(9,m)
                  else
                     energies(ns,3) = etab(6,m)
                     menergies(ns,3) = mtab(6,m)
                  end if
                  
               end if

            else if(ndir.eq.3) then
c-----------check above
c..............DVIL1
               if(system(ii,jj,ml(kk+1)).eq.1 .and. 
     &              system(mb(ii-1),md(jj+1),ml(kk+1)).eq.2) then

                  iv = ivcoord(system,md(ii-1),md(jj+1),ml(kk+1))
                  if(iv.eq.2) then
                     energies(ns,2) = etab(7,m)
                     menergies(ns,2) = mtab(7,m)
                  else
                     energies(ns,2) = etab(4,m)
                     menergies(ns,2) = mtab(4,m)
                  end if

c..............DVIL2
               else if(system(ii,jj,ml(kk+1)).eq.1 .and. 
     &                 (system(mb(ii+1),md(jj+1),ml(kk+1)).eq.2 .or.
     &                 system(ii,md(jj-2),ml(kk+1)).eq.2)) then

                  if(system(mb(ii+1),md(jj+1),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii+1),md(jj+1),ml(kk+1))
                  else if(system(ii,md(jj-2),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj-2),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(8,m)
                     menergies(ns,2) = mtab(8,m)
                  else
                     energies(ns,2) = etab(5,m)
                     menergies(ns,2) = mtab(5,m)
                  end if

c..............DVIL3
               else if((system(ii,jj,ml(kk+1)).eq.2 .and. 
     &                 system(mb(ii-1),md(jj+1),ml(kk+1)).eq.1) .or.
     &                 (system(mb(ii+1),md(jj-1),ml(kk+1)).eq.1 .and.
     &                 system(mb(ii+2),md(jj-2),ml(kk+1)).eq.2)) then

                  if(system(ii,jj,ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk+1))
                  else if(system(mb(ii+2),md(jj-2),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii+2),md(jj-2),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(9,m)
                     menergies(ns,2) = mtab(9,m)
                  else
                     energies(ns,2) = etab(6,m)
                     menergies(ns,2) = mtab(6,m)
                  end if

               end if
c-----------check below
c..............DVIL1
               if(system(ii,jj,ml(kk-1)).eq.1 .and. 
     &              system(mb(ii-1),md(jj+1),ml(kk-1)).eq.2) then

                  iv = ivcoord(system,md(ii-1),md(jj+1),ml(kk-1))
                  if(iv.eq.2) then
                     energies(ns,3) = etab(7,m)
                     menergies(ns,3) = mtab(7,m)
                  else
                     energies(ns,3) = etab(4,m)
                     menergies(ns,3) = mtab(4,m)
                  end if

c..............DVIL2
               else if(system(ii,jj,ml(kk-1)).eq.1 .and. 
     &                 (system(mb(ii+1),md(jj+1),ml(kk-1)).eq.2 .or.
     &                 system(ii,md(jj-2),ml(kk-1)).eq.2)) then

                  if(system(mb(ii+1),md(jj+1),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii+1),md(jj+1),ml(kk-1))
                  else if(system(ii,md(jj-2),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj-2),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(8,m)
                     menergies(ns,3) = mtab(8,m)
                  else
                     energies(ns,3) = etab(5,m)
                     menergies(ns,3) = mtab(5,m)
                  end if

c..............DVIL3
               else if((system(ii,jj,ml(kk-1)).eq.2 .and. 
     &                 system(mb(ii-1),md(jj+1),ml(kk-1)).eq.1) .or.
     &                 (system(mb(ii+1),md(jj-1),ml(kk-1)).eq.1 .and.
     &                 system(mb(ii+2),md(jj-2),ml(kk-1)).eq.2)) then

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk-1))
                  else if(system(mb(ii+2),md(jj-2),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii+2),md(jj-2),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(9,m)
                     menergies(ns,3) = mtab(9,m)
                  else
                     energies(ns,3) = etab(6,m)
                     menergies(ns,3) = mtab(6,m)
                  end if

               end if

            else if(ndir.eq.4) then
c-----------check above
c..............DVIL1
               if(system(ii,jj,ml(kk+1)).eq.1 .and. 
     &              system(ii,md(jj+2),ml(kk+1)).eq.2) then

                  iv = ivcoord(system,ii,md(jj+2),ml(kk+1))
                  if(iv.eq.2) then
                     energies(ns,2) = etab(7,m)
                     menergies(ns,2) = mtab(7,m)
                  else
                     energies(ns,2) = etab(4,m)
                     menergies(ns,2) = mtab(4,m)
                  end if

c..............DVIL2
               else if(system(ii,jj,ml(kk+1)).eq.1 .and. 
     &                 (system(mb(ii+1),md(jj-1),ml(kk+1)).eq.2 .or.
     &                 system(mb(ii-1),md(jj-1),ml(kk+1)).eq.2)) then

                  if(system(mb(ii+1),md(jj-1),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii+1),md(jj-1),ml(kk+1))
                  else if(system(mb(ii-1),md(jj-1),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii-1),md(jj-1),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(8,m)
                     menergies(ns,2) = mtab(8,m)
                  else
                     energies(ns,2) = etab(5,m)
                     menergies(ns,2) = mtab(5,m)
                  end if

c..............DVIL3
               else if((system(ii,jj,ml(kk+1)).eq.2 .and. 
     &                 system(mb(ii),md(jj+2),ml(kk+1)).eq.1) .or. 
     &                 (system(ii,md(jj-2),ml(kk+1)).eq.1 .and.
     &                 system(ii,md(jj-4),ml(kk+1)).eq.2)) then

                  if(system(ii,jj,ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk+1))
                  else if(system(ii,md(jj-4),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj-4),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(9,m)
                     menergies(ns,2) = mtab(9,m)
                  else
                     energies(ns,2) = etab(6,m)
                     menergies(ns,2) = mtab(6,m)
                  end if

               end if
c-----------check below
               if(system(ii,jj,ml(kk-1)).eq.1 .and. 
     &              system(ii,md(jj+2),ml(kk-1)).eq.2) then

                  iv = ivcoord(system,ii,md(jj+2),ml(kk-1))
                  if(iv.eq.2) then
                     energies(ns,3) = etab(7,m)
                     menergies(ns,3) = mtab(7,m)
                  else
                     energies(ns,3) = etab(4,m)
                     menergies(ns,3) = mtab(4,m)
                  end if

c..............DVIL2
               else if(system(ii,jj,ml(kk-1)).eq.1 .and. 
     &                 (system(mb(ii+1),md(jj-1),ml(kk-1)).eq.2 .or.
     &                 system(mb(ii-1),md(jj-1),ml(kk-1)).eq.2)) then

                  if(system(mb(ii+1),md(jj-1),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii+1),md(jj-1),ml(kk-1))
                  else if(system(mb(ii-1),md(jj-1),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii-1),md(jj-1),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(8,m)
                     menergies(ns,3) = mtab(8,m)
                  else
                     energies(ns,3) = etab(5,m)
                     menergies(ns,3) = mtab(5,m)
                  end if

c..............DVIL3
               else if((system(ii,jj,ml(kk-1)).eq.2 .and. 
     &                 system(mb(ii),md(jj+2),ml(kk-1)).eq.1) .or. 
     &                 (system(ii,md(jj-2),ml(kk-1)).eq.1 .and.
     &                 system(ii,md(jj-4),ml(kk-1)).eq.2)) then

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk-1))
                  else if(system(ii,md(jj-4),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj-4),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(9,m)
                     menergies(ns,3) = mtab(9,m)
                  else
                     energies(ns,3) = etab(6,m)
                     menergies(ns,3) = mtab(6,m)
                  end if

               end if
               
            else if(ndir.eq.5) then
c-----------check above
c..............DVIL1
               if(system(ii,jj,ml(kk+1)).eq.1 .and. 
     &              system(mb(ii+1),md(jj+1),ml(kk+1)).eq.2) then

                  iv = ivcoord(system,md(ii+1),md(jj+1),ml(kk+1))
                  if(iv.eq.2) then
                     energies(ns,2) = etab(7,m)
                     menergies(ns,2) = mtab(7,m)
                  else
                     energies(ns,2) = etab(4,m)
                     menergies(ns,2) = mtab(4,m)
                  end if
                  
c..............DVIL2
               else if(system(ii,jj,ml(kk+1)).eq.1 .and. 
     &                 (system(mb(ii-1),md(jj+1),ml(kk+1)).eq.2 .or.
     &                 system(ii,md(jj-2),ml(kk+1)).eq.2)) then

                  if(system(mb(ii-1),md(jj+1),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii-1),md(jj+1),ml(kk+1))
                  else if(system(ii,md(jj-2),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj-2),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(8,m)
                     menergies(ns,2) = mtab(8,m)
                  else
                     energies(ns,2) = etab(5,m)
                     menergies(ns,2) = mtab(5,m)
                  end if

c..............DVIL3
               else if((system(ii,jj,ml(kk+1)).eq.2 .and. 
     &                 system(mb(ii+1),md(jj+1),ml(kk+1)).eq.1) .or.
     &                 (system(mb(ii-1),md(jj-1),ml(kk+1)).eq.1 .and.
     &                 system(mb(ii-2),md(jj-2),ml(kk+1)).eq.2)) then

                  if(system(ii,jj,ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk+1))
                  else if(system(mb(ii-2),md(jj-2),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii-2),md(jj-2),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(9,m)
                     menergies(ns,2) = mtab(9,m)
                  else
                     energies(ns,2) = etab(6,m)
                     menergies(ns,2) = mtab(6,m)
                  end if
                  
               end if
c-----------check below
c..............DVIL1
               if(system(ii,jj,ml(kk-1)).eq.1 .and. 
     &              system(mb(ii+1),md(jj+1),ml(kk-1)).eq.2) then

                  iv = ivcoord(system,md(ii+1),md(jj+1),ml(kk-1))
                  if(iv.eq.2) then
                     energies(ns,3) = etab(7,m)
                     menergies(ns,3) = mtab(7,m)
                  else
                     energies(ns,3) = etab(4,m)
                     menergies(ns,3) = mtab(4,m)
                  end if

c..............DVIL2
               else if(system(ii,jj,ml(kk-1)).eq.1 .and. 
     &                 (system(mb(ii-1),md(jj+1),ml(kk-1)).eq.2 .or.
     &                 system(ii,md(jj-2),ml(kk-1)).eq.2)) then

                  if(system(mb(ii-1),md(jj+1),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii-1),md(jj+1),ml(kk-1))
                  else if(system(ii,md(jj-2),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj-2),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(8,m)
                     menergies(ns,3) = mtab(8,m)
                  else
                     energies(ns,3) = etab(5,m)
                     menergies(ns,3) = mtab(5,m)
                  end if

c..............DVIL3
               else if((system(ii,jj,ml(kk-1)).eq.2 .and. 
     &                 system(mb(ii+1),md(jj+1),ml(kk-1)).eq.1) .or.
     &                 (system(mb(ii-1),md(jj-1),ml(kk-1)).eq.1 .and.
     &                 system(mb(ii-2),md(jj-2),ml(kk-1)).eq.2)) then

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk-1))
                  else if(system(mb(ii-2),md(jj-2),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii-2),md(jj-2),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(9,m)
                     menergies(ns,3) = mtab(9,m)
                  else
                     energies(ns,3) = etab(6,m)
                     menergies(ns,3) = mtab(6,m)
                  end if

               end if

            else if(ndir.eq.6) then
c-----------check above
c..............DVIL1
               if(system(ii,jj,ml(kk+1)).eq.1 .and. 
     &              system(mb(ii+1),md(jj-1),ml(kk+1)).eq.2) then

                  iv = ivcoord(system,md(ii+1),md(jj-1),ml(kk+1))
                  if(iv.eq.2) then
                     energies(ns,2) = etab(7,m)
                     menergies(ns,2) = mtab(7,m)
                  else
                     energies(ns,2) = etab(4,m)
                     menergies(ns,2) = mtab(4,m)
                  end if

c..............DVIL2
               else if(system(ii,jj,ml(kk+1)).eq.1 .and. 
     &                 (system(mb(ii-1),md(jj-1),ml(kk+1)).eq.2 .or.
     &                 system(ii,md(jj+2),ml(kk+1)).eq.2)) then

                  if(system(mb(ii-1),md(jj-1),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii-1),md(jj-1),ml(kk+1))
                  else if(system(ii,md(jj+2),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj+2),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(8,m)
                     menergies(ns,2) = mtab(8,m)
                  else
                     energies(ns,2) = etab(5,m)
                     menergies(ns,2) = mtab(5,m)
                  end if

c..............DVIL3
               else if((system(ii,jj,ml(kk+1)).eq.2 .and. 
     &                 system(mb(ii+1),md(jj-1),ml(kk+1)).eq.1) .or.
     &                 (system(mb(ii-1),md(jj+1),ml(kk+1)).eq.1 .and.
     &                 system(mb(ii-2),md(jj+2),ml(kk+1)).eq.2)) then

                  if(system(ii,jj,ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk+1))
                  else if(system(mb(ii-2),md(jj+2),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii-2),md(jj+2),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(9,m)
                     menergies(ns,2) = mtab(9,m)
                  else
                     energies(ns,2) = etab(6,m)
                     menergies(ns,2) = mtab(6,m)
                  end if

               end if
c-----------check below
c..............DVIL1
               if(system(ii,jj,ml(kk-1)).eq.1 .and. 
     &              system(mb(ii+1),md(jj-1),ml(kk-1)).eq.2) then

                  iv = ivcoord(system,md(ii+1),md(jj-1),ml(kk-1))
                  if(iv.eq.2) then
                     energies(ns,3) = etab(7,m)
                     menergies(ns,3) = mtab(7,m)
                  else
                     energies(ns,3) = etab(4,m)
                     menergies(ns,3) = mtab(4,m)
                  end if
                  
c..............DVIL2
               else if(system(ii,jj,ml(kk-1)).eq.1 .and. 
     &                 (system(mb(ii-1),md(jj-1),ml(kk-1)).eq.2 .or.
     &                 system(ii,md(jj+2),ml(kk-1)).eq.2)) then

                  if(system(mb(ii-1),md(jj-1),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii-1),md(jj-1),ml(kk-1))
                  else if(system(ii,md(jj+2),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj+2),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(8,m)
                     menergies(ns,3) = mtab(8,m)
                  else
                     energies(ns,3) = etab(5,m)
                     menergies(ns,3) = mtab(5,m)
                  end if

c..............DVIL3
               else if((system(ii,jj,ml(kk-1)).eq.2 .and. 
     &                 system(mb(ii+1),md(jj-1),ml(kk-1)).eq.1) .or.
     &                 (system(mb(ii-1),md(jj+1),ml(kk-1)).eq.1 .and.
     &                 system(mb(ii-2),md(jj+2),ml(kk-1)).eq.2)) then

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk-1))
                  else if(system(mb(ii-2),md(jj+2),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii-2),md(jj+2),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(9,m)
                     menergies(ns,3) = mtab(9,m)
                  else
                     energies(ns,3) = etab(6,m)
                     menergies(ns,3) = mtab(6,m)
                  end if

               end if
            end if



c--------then, if the central atom is bound to the lower layer
c        (this adds the barrier to break the bond - move atom up)
c        (and the barrier to move an atom between the layers)

         else if(system(ii,jj,kk).eq.6) then
c--------find which type of interlayer bond
            
            if(ndir.eq.1) then
c..............DVIL1
               if(system(ii,jj,ml(kk-1)).eq.5 .and. 
     &              system(ii,md(jj-2),ml(kk-1)).eq.2) then

                  iv = ivcoord(system,ii,md(jj-2),ml(kk-1))
                  if(iv.eq.2) then
                     energies(ns,2) = etab(7,m)
                     menergies(ns,2) = mtab(7,m)
                  else
                     energies(ns,2) = etab(4,m)
                     menergies(ns,2) = mtab(4,m)
                  end if

                  energies(ns,5) = etab(10,m)
                  menergies(ns,5) = mtab(10,m)
c..............DVIL2
               else if(system(ii,jj,ml(kk-1)).eq.5 .and. 
     &                 (system(mb(ii-1),md(jj+1),ml(kk-1)).eq.2 .or.
     &                 system(mb(ii+1),md(jj+1),ml(kk-1)).eq.2)) then

                  if(system(mb(ii-1),md(jj+1),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii-1),md(jj+1),ml(kk-1))
                  else if(system(mb(ii+1),md(jj+1),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii+1),md(jj+1),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(8,m)
                     menergies(ns,2) = mtab(8,m)
                  else
                     energies(ns,2) = etab(5,m)
                     menergies(ns,2) = mtab(5,m)
                  end if

                  energies(ns,5) = etab(11,m)
                  menergies(ns,5) = mtab(11,m)
c..............DVIL3
               else

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk-1))
                  else if(system(ii,md(jj+4),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj+4),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(9,m)
                     menergies(ns,2) = mtab(9,m)
                  else
                     energies(ns,2) = etab(6,m)
                     menergies(ns,2) = mtab(6,m)
                  end if

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     energies(ns,5) = etab(12,m)
                     menergies(ns,5) = mtab(12,m)
                  end if
               end if
               
            else if(ndir.eq.2) then
c-----------check above
c..............DVIL1
               if(system(ii,jj,ml(kk-1)).eq.5 .and. 
     &              system(mb(ii-1),md(jj-1),ml(kk-1)).eq.2) then

                  iv = ivcoord(system,md(ii-1),md(jj-1),ml(kk-1))
                  if(iv.eq.2) then
                     energies(ns,2) = etab(7,m)
                     menergies(ns,2) = mtab(7,m)
                  else
                     energies(ns,2) = etab(4,m)
                     menergies(ns,2) = mtab(4,m)
                  end if

                  energies(ns,5) = etab(10,m)
                  menergies(ns,5) = mtab(10,m)
c..............DVIL2
               else if(system(ii,jj,ml(kk-1)).eq.5 .and. 
     &                 (system(mb(ii+1),md(jj-1),ml(kk-1)).eq.2 .or.
     &                 system(ii,md(jj+2),ml(kk-1)).eq.2)) then

                  if(system(mb(ii+1),md(jj-1),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii+1),md(jj-1),ml(kk-1))
                  else if(system(ii,md(jj+2),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj+2),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(8,m)
                     menergies(ns,2) = mtab(8,m)
                  else
                     energies(ns,2) = etab(5,m)
                     menergies(ns,2) = mtab(5,m)
                  end if

                  energies(ns,5) = etab(11,m)
                  menergies(ns,5) = mtab(11,m)
c..............DVIL3
               else 

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk-1))
                  else if(system(mb(ii+2),md(jj+2),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii+2),md(jj+2),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(9,m)
                     menergies(ns,2) = mtab(9,m)
                  else
                     energies(ns,2) = etab(6,m)
                     menergies(ns,2) = mtab(6,m)
                  end if

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     energies(ns,5) = etab(12,m)
                     menergies(ns,5) = mtab(12,m)
                  end if
               end if

            else if(ndir.eq.3) then
c-----------check above
c..............DVIL1
               if(system(ii,jj,ml(kk-1)).eq.5 .and. 
     &              system(mb(ii-1),md(jj+1),ml(kk-1)).eq.2) then

                  iv = ivcoord(system,md(ii-1),md(jj+1),ml(kk-1))
                  if(iv.eq.2) then
                     energies(ns,2) = etab(7,m)
                     menergies(ns,2) = mtab(7,m)
                  else
                     energies(ns,2) = etab(4,m)
                     menergies(ns,2) = mtab(4,m)
                  end if

                 energies(ns,5) = etab(10,m)
c..............DVIL2
               else if(system(ii,jj,ml(kk-1)).eq.5 .and. 
     &                 (system(mb(ii+1),md(jj+1),ml(kk-1)).eq.2 .or.
     &                 system(ii,md(jj-2),ml(kk-1)).eq.2)) then

                  if(system(mb(ii+1),md(jj+1),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii+1),md(jj+1),ml(kk-1))
                  else if(system(ii,md(jj-2),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj-2),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(8,m)
                     menergies(ns,2) = mtab(8,m)
                  else
                     energies(ns,2) = etab(5,m)
                     menergies(ns,2) = mtab(5,m)
                  end if

                  energies(ns,5) = etab(11,m)
                  menergies(ns,5) = mtab(11,m)
c..............DVIL3
               else

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk-1))
                  else if(system(mb(ii+2),md(jj-2),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii+2),md(jj-2),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(9,m)
                     menergies(ns,2) = mtab(9,m)
                  else
                     energies(ns,2) = etab(6,m)
                     menergies(ns,2) = mtab(6,m)
                  end if

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     energies(ns,5) = etab(12,m)
                     menergies(ns,5) = mtab(12,m)
                  end if
               end if

            else if(ndir.eq.4) then
c-----------check above
c..............DVIL1
               if(system(ii,jj,ml(kk-1)).eq.5 .and. 
     &              system(ii,md(jj+2),ml(kk-1)).eq.2) then

                  iv = ivcoord(system,ii,md(jj+2),ml(kk-1))
                  if(iv.eq.2) then
                     energies(ns,2) = etab(7,m)
                     menergies(ns,2) = mtab(7,m)
                  else
                     energies(ns,2) = etab(4,m)
                     menergies(ns,2) = mtab(4,m)
                  end if

                  energies(ns,4) = etab(10,m)
                  menergies(ns,4) = mtab(10,m)
c..............DVIL2
               else if(system(ii,jj,ml(kk-1)).eq.5 .and. 
     &                 (system(mb(ii+1),md(jj-1),ml(kk-1)).eq.2 .or.
     &                 system(mb(ii-1),md(jj-1),ml(kk-1)).eq.2)) then

                  if(system(mb(ii+1),md(jj-1),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii+1),md(jj-1),ml(kk-1))
                  else if(system(mb(ii-1),md(jj-1),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii-1),md(jj-1),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(8,m)
                     menergies(ns,2) = mtab(8,m)
                  else
                     energies(ns,2) = etab(5,m)
                     menergies(ns,2) = mtab(5,m)
                  end if
                  
                  energies(ns,5) = etab(11,m)
                  menergies(ns,5) = mtab(11,m)
c..............DVIL3
               else

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk-1))
                  else if(system(ii,md(jj-4),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj-4),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(9,m)
                     menergies(ns,2) = mtab(9,m)
                  else
                     energies(ns,2) = etab(6,m)
                     menergies(ns,2) = mtab(6,m)
                  end if

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     energies(ns,5) = etab(12,m)
                     menergies(ns,5) = mtab(12,m)
                  end if
               end if
               
            else if(ndir.eq.5) then
c-----------check above
c..............DVIL1
               if(system(ii,jj,ml(kk-1)).eq.5 .and. 
     &              system(mb(ii+1),md(jj+1),ml(kk-1)).eq.2) then

                  iv = ivcoord(system,md(ii+1),md(jj+1),ml(kk-1))
                  if(iv.eq.2) then
                     energies(ns,2) = etab(7,m)
                     menergies(ns,2) = mtab(7,m)
                  else
                     energies(ns,2) = etab(4,m)
                     menergies(ns,2) = mtab(4,m)
                  end if

                  energies(ns,5) = etab(10,m)
                  menergies(ns,5) = mtab(10,m)
c..............DVIL2
               else if(system(ii,jj,ml(kk-1)).eq.5 .and. 
     &                 (system(mb(ii-1),md(jj+1),ml(kk-1)).eq.2 .or.
     &                 system(ii,md(jj-2),ml(kk-1)).eq.2)) then

                  if(system(mb(ii-1),md(jj+1),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii-1),md(jj+1),ml(kk-1))
                  else if(system(ii,md(jj-2),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj-2),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(8,m)
                     menergies(ns,2) = mtab(8,m)
                  else
                     energies(ns,2) = etab(5,m)
                     menergies(ns,2) = mtab(5,m)
                  end if

                  energies(ns,5) = etab(11,m)
                  menergies(ns,5) = mtab(11,m)
c..............DVIL3
               else

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk-1))
                  else if(system(mb(ii-2),md(jj-2),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii-2),md(jj-2),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(9,m)
                     menergies(ns,2) = mtab(9,m)
                  else
                     energies(ns,2) = etab(6,m)
                     menergies(ns,2) = mtab(6,m)
                  end if

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     energies(ns,5) = etab(12,m)
                     menergies(ns,5) = mtab(12,m)
                  end if
               end if

            else if(ndir.eq.6) then
c-----------check above
c..............DVIL1
               if(system(ii,jj,ml(kk-1)).eq.5 .and. 
     &              system(mb(ii+1),md(jj-1),ml(kk-1)).eq.2) then

                  iv = ivcoord(system,md(ii+1),md(jj-1),ml(kk-1))
                  if(iv.eq.2) then
                     energies(ns,2) = etab(7,m)
                     menergies(ns,2) = mtab(7,m)
                  else
                     energies(ns,2) = etab(4,m)
                     menergies(ns,2) = mtab(4,m)
                  end if

                  energies(ns,5) = etab(10,m)
                  menergies(ns,5) = mtab(10,m)
c..............DVIL2
               else if(system(ii,jj,ml(kk-1)).eq.5 .and. 
     &                 (system(mb(ii-1),md(jj-1),ml(kk-1)).eq.2 .or.
     &                 system(ii,md(jj+2),ml(kk-1)).eq.2)) then

                  if(system(mb(ii-1),md(jj-1),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii-1),md(jj-1),ml(kk-1))
                  else if(system(ii,md(jj+2),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj+2),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(8,m)
                     menergies(ns,2) = mtab(8,m)
                  else
                     energies(ns,2) = etab(5,m)
                     menergies(ns,2) = mtab(5,m)
                  end if

                  energies(ns,5) = etab(11,m)
                  menergies(ns,5) = mtab(11,m)
c..............DVIL3
               else

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk-1))
                  else if(system(mb(ii-2),md(jj+2),ml(kk-1)).eq.2) then
                     iv = ivcoord(system,md(ii-2),md(jj+2),ml(kk-1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,2) = etab(9,m)
                     menergies(ns,2) = mtab(9,m)
                  else
                     energies(ns,2) = etab(6,m)
                     menergies(ns,2) = mtab(6,m)
                  end if

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     energies(ns,5) = etab(12,m)
                     menergies(ns,5) = mtab(12,m)
                  end if
               end if
            end if

c--------then, if the central atom is bound to the upper layer
c        (this adds the barrier to break the bond - move atom down)
c        (and the barrier to move an atom between the layers)

         else if(system(ii,jj,kk).eq.5) then
c--------find which type of interlayer bond
            
            if(ndir.eq.1) then
c..............DVIL1
               if(system(ii,jj,ml(kk+1)).eq.6 .and. 
     &              system(ii,md(jj-2),ml(kk+1)).eq.2) then

                  iv = ivcoord(system,ii,md(jj-2),ml(kk+1))
                  if(iv.eq.2) then
                     energies(ns,3) = etab(7,m)
                     menergies(ns,3) = mtab(7,m)
                  else
                     energies(ns,3) = etab(4,m)
                     menergies(ns,3) = mtab(4,m)
                  end if

                  energies(ns,4) = etab(10,m)
                  menergies(ns,4) = mtab(10,m)
c..............DVIL2
               else if(system(ii,jj,ml(kk+1)).eq.6 .and. 
     &                 (system(mb(ii-1),md(jj+1),ml(kk+1)).eq.2 .or.
     &                 system(mb(ii+1),md(jj+1),ml(kk+1)).eq.2)) then

                  if(system(mb(ii-1),md(jj+1),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii-1),md(jj+1),ml(kk+1))
                  else if(system(mb(ii+1),md(jj+1),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii+1),md(jj+1),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(8,m)
                     menergies(ns,3) = mtab(8,m)
                  else
                     energies(ns,3) = etab(5,m)
                     menergies(ns,3) = mtab(5,m)
                  end if

                  energies(ns,4) = etab(11,m)
                  menergies(ns,4) = mtab(11,m)
c..............DVIL3
               else

                  if(system(ii,jj,ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk+1))
                  else if(system(ii,md(jj+4),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj+4),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(9,m)
                     menergies(ns,3) = mtab(9,m)
                  else
                     energies(ns,3) = etab(6,m)
                     menergies(ns,3) = mtab(6,m)
                  end if

                  if(system(ii,jj,ml(kk+1)).eq.2) then
                     energies(ns,4) = etab(12,m)
                     menergies(ns,4) = mtab(12,m)
                  end if
               end if
               
            else if(ndir.eq.2) then
c-----------check above
c..............DVIL1
               if(system(ii,jj,ml(kk+1)).eq.6 .and. 
     &              system(mb(ii-1),md(jj-1),ml(kk+1)).eq.2) then

                  iv = ivcoord(system,md(ii-1),md(jj-1),ml(kk+1))
                  if(iv.eq.2) then
                     energies(ns,3) = etab(7,m)
                     menergies(ns,3) = mtab(7,m)
                  else
                     energies(ns,3) = etab(4,m)
                     menergies(ns,3) = mtab(4,m)
                  end if

                  energies(ns,4) = etab(10,m)
                  menergies(ns,4) = mtab(10,m)
c..............DVIL2
               else if(system(ii,jj,ml(kk+1)).eq.6 .and. 
     &                 (system(mb(ii+1),md(jj-1),ml(kk+1)).eq.2 .or.
     &                 system(ii,md(jj+2),ml(kk+1)).eq.2)) then

                  if(system(mb(ii+1),md(jj-1),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii+1),md(jj-1),ml(kk+1))
                  else if(system(ii,md(jj+2),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj+2),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(8,m)
                     menergies(ns,3) = mtab(8,m)
                  else
                     energies(ns,3) = etab(5,m)
                     menergies(ns,3) = mtab(5,m)
                  end if

                  energies(ns,4) = etab(11,m)
                  menergies(ns,4) = mtab(11,m)
c..............DVIL3
               else 

                  if(system(ii,jj,ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk+1))
                  else if(system(mb(ii+2),md(jj+2),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii+2),md(jj+2),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(9,m)
                     menergies(ns,3) = mtab(9,m)
                  else
                     energies(ns,3) = etab(6,m)
                     menergies(ns,3) = mtab(6,m)
                  end if

                  if(system(ii,jj,ml(kk+1)).eq.2) then
                     energies(ns,4) = etab(12,m)
                     menergies(ns,4) = mtab(12,m)
                  end if
               end if

            else if(ndir.eq.3) then
c-----------check above
c..............DVIL1
               if(system(ii,jj,ml(kk+1)).eq.6 .and. 
     &              system(mb(ii-1),md(jj+1),ml(kk+1)).eq.2) then

                  iv = ivcoord(system,md(ii-1),md(jj+1),ml(kk+1))
                  if(iv.eq.2) then
                     energies(ns,3) = etab(7,m)
                     menergies(ns,3) = mtab(7,m)
                  else
                     energies(ns,3) = etab(4,m)
                     menergies(ns,3) = mtab(4,m)
                  end if

                  energies(ns,4) = etab(10,m)
                  menergies(ns,4) = mtab(10,m)
c..............DVIL2
               else if(system(ii,jj,ml(kk+1)).eq.6 .and. 
     &                 (system(mb(ii+1),md(jj+1),ml(kk+1)).eq.2 .or.
     &                 system(ii,md(jj-2),ml(kk+1)).eq.2)) then

                  if(system(mb(ii+1),md(jj+1),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii+1),md(jj+1),ml(kk+1))
                  else if(system(ii,md(jj-2),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj-2),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(8,m)
                     menergies(ns,3) = mtab(8,m)
                  else
                     energies(ns,3) = etab(5,m)
                     menergies(ns,3) = mtab(5,m)
                  end if

                  energies(ns,4) = etab(11,m)
                  menergies(ns,4) = mtab(11,m)
c..............DVIL3
               else

                  if(system(ii,jj,ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk+1))
                  else if(system(mb(ii+2),md(jj-2),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii+2),md(jj-2),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(9,m)
                     menergies(ns,3) = mtab(9,m)
                  else
                     energies(ns,3) = etab(6,m)
                     menergies(ns,3) = mtab(6,m)
                  end if

                  if(system(ii,jj,ml(kk+1)).eq.2) then
                     energies(ns,4) = etab(12,m)
                     menergies(ns,4) = mtab(12,m)
                  end if
               end if

            else if(ndir.eq.4) then
c-----------check above
c..............DVIL1
               if(system(ii,jj,ml(kk+1)).eq.6 .and. 
     &              system(ii,md(jj+2),ml(kk+1)).eq.2) then

                  iv = ivcoord(system,ii,md(jj+2),ml(kk+1))
                  if(iv.eq.2) then
                     energies(ns,3) = etab(7,m)
                     menergies(ns,3) = mtab(7,m)
                  else
                     energies(ns,3) = etab(4,m)
                     menergies(ns,3) = mtab(4,m)
                  end if

                  energies(ns,4) = etab(10,m)
                  menergies(ns,4) = mtab(10,m)
c..............DVIL2
               else if(system(ii,jj,ml(kk+1)).eq.6 .and. 
     &                 (system(mb(ii+1),md(jj-1),ml(kk+1)).eq.2 .or.
     &                 system(mb(ii-1),md(jj-1),ml(kk+1)).eq.2)) then

                  if(system(mb(ii+1),md(jj-1),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii+1),md(jj-1),ml(kk+1))
                  else if(system(mb(ii-1),md(jj-1),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii-1),md(jj-1),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(8,m)
                     menergies(ns,3) = mtab(8,m)
                  else
                     energies(ns,3) = etab(5,m)
                     menergies(ns,3) = mtab(5,m)
                  end if

                  energies(ns,4) = etab(11,m)
                  menergies(ns,4) = mtab(11,m)
c..............DVIL3
               else
                  
                  if(system(ii,jj,ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk+1))
                  else if(system(ii,md(jj-4),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj-4),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(9,m)
                     menergies(ns,3) = mtab(9,m)
                  else
                     energies(ns,3) = etab(6,m)
                     menergies(ns,3) = mtab(6,m)
                  end if

                  if(system(ii,jj,ml(kk+1)).eq.2) then
                     energies(ns,4) = etab(12,m)
                     menergies(ns,4) = mtab(12,m)
                  end if
               end if
               
            else if(ndir.eq.5) then
c-----------check above
c..............DVIL1
               if(system(ii,jj,ml(kk+1)).eq.6 .and. 
     &              system(mb(ii+1),md(jj+1),ml(kk+1)).eq.2) then

                  iv = ivcoord(system,md(ii+1),md(jj+1),ml(kk+1))
                  if(iv.eq.2) then
                     energies(ns,3) = etab(7,m)
                     menergies(ns,3) = mtab(7,m)
                  else
                     energies(ns,3) = etab(4,m)
                     menergies(ns,3) = mtab(4,m)
                  end if

                  energies(ns,4) = etab(10,m)
                  menergies(ns,4) = mtab(10,m)
c..............DVIL2
               else if(system(ii,jj,ml(kk+1)).eq.6 .and. 
     &                 (system(mb(ii-1),md(jj+1),ml(kk+1)).eq.2 .or.
     &                 system(ii,md(jj-2),ml(kk+1)).eq.2)) then

                  if(system(mb(ii-1),md(jj+1),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii-1),md(jj+1),ml(kk+1))
                  else if(system(ii,md(jj-2),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj-2),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(8,m)
                     menergies(ns,3) = mtab(8,m)
                  else
                     energies(ns,3) = etab(5,m)
                     menergies(ns,3) = mtab(5,m)
                  end if

                  energies(ns,4) = etab(11,m)
                  menergies(ns,4) = mtab(11,m)
c..............DVIL3
               else

                  if(system(ii,jj,ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk+1))
                  else if(system(mb(ii-2),md(jj-2),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii-2),md(jj-2),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(9,m)
                     menergies(ns,3) = mtab(9,m)
                  else
                     energies(ns,3) = etab(6,m)
                     menergies(ns,3) = mtab(6,m)
                  end if

                  if(system(ii,jj,ml(kk+1)).eq.2) then
                     energies(ns,4) = etab(12,m)
                     menergies(ns,4) = mtab(12,m)
                  end if
               end if

            else if(ndir.eq.6) then
c-----------check above
c..............DVIL1
               if(system(ii,jj,ml(kk+1)).eq.6 .and. 
     &              system(mb(ii+1),md(jj-1),ml(kk+1)).eq.2) then

                  iv = ivcoord(system,md(ii+1),md(jj-1),ml(kk+1))
                  if(iv.eq.2) then
                     energies(ns,3) = etab(7,m)
                     menergies(ns,3) = mtab(7,m)
                  else
                     energies(ns,3) = etab(4,m)
                     menergies(ns,3) = mtab(4,m)
                  end if

                  energies(ns,4) = etab(10,m)
                  menergies(ns,4) = mtab(10,m)
c..............DVIL2
               else if(system(ii,jj,ml(kk+1)).eq.6 .and. 
     &                 (system(mb(ii-1),md(jj-1),ml(kk+1)).eq.2 .or.
     &                 system(ii,md(jj+2),ml(kk+1)).eq.2)) then

                  if(system(mb(ii-1),md(jj-1),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii-1),md(jj-1),ml(kk+1))
                  else if(system(ii,md(jj+2),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,ii,md(jj+2),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(8,m)
                     menergies(ns,3) = mtab(8,m)
                  else
                     energies(ns,3) = etab(5,m)
                     menergies(ns,3) = mtab(5,m)
                  end if

                  energies(ns,4) = etab(11,m)
                  menergies(ns,4) = mtab(11,m)
c..............DVIL3
               else

                  if(system(ii,jj,ml(kk-1)).eq.2) then
                     iv = ivcoord(system,ii,jj,ml(kk+1))
                  else if(system(mb(ii-2),md(jj+2),ml(kk+1)).eq.2) then
                     iv = ivcoord(system,md(ii-2),md(jj+2),ml(kk+1))
                  end if
                  if(iv.eq.2) then
                     energies(ns,3) = etab(9,m)
                     menergies(ns,3) = mtab(9,m)
                  else
                     energies(ns,3) = etab(6,m)
                     menergies(ns,3) = mtab(6,m)
                  end if
                  
                  if(system(ii,jj,ml(kk+1)).eq.2) then
                     energies(ns,4) = etab(12,m)
                     menergies(ns,4) = mtab(12,m)
                  end if
               end if
            end if 
         end if

      end do

c      if(n_s.ne.nsites) then
c         write(*,*)'Error in sites'
c         stop
c      end if

c      if(n_s.eq.0) then
c         write(*,*)'Ground state'
c         stop
c      end if
      
      end

      integer function ivcoord(system,i,j,k)
c------------------------------------------------------------
c     coordination number of a vacancy at i,j,k             - 
c------------------------------------------------------------
      integer*1 system(800,800,199)
      
      iv = 0
      
c.....if up

      if(system(i,md(j+2),k).ne.0) then
         
         if(system(i,md(j+2),k).eq.2) iv = iv + 1
         if(system(mb(i+1),md(j-1),k).eq.2) iv = iv + 1
         if(system(mb(i-1),md(j-1),k).eq.2) iv = iv + 1

      else

         if(system(i,md(j-2),k).eq.2) iv = iv + 1
         if(system(mb(i+1),md(j+1),k).eq.2) iv = iv + 1
         if(system(mb(i-1),md(j+1),k).eq.2) iv = iv + 1         

      end if

      if(iv.eq.3) then
         write(*,*)'error in ivcoord'
         stop
      end if

      ivcoord = 3 - iv

      end


      integer function ivcoord2(system,i,j,k)
c------------------------------------------------------------
c     coordination number of a vacancy at i,j,k             - 
c------------------------------------------------------------
      integer*1 system(800,800,199)
      
      iv = 0
      
c.....if up

      if(system(mb(i),md(j+2),k).ne.0) then
         
         if(system(mb(i),md(j+2),k).eq.1) iv = iv + 1
         if(system(mb(i+1),md(j-1),k).eq.1) iv = iv + 1
         if(system(mb(i-1),md(j-1),k).eq.1) iv = iv + 1

      else

         if(system(mb(i),md(j-2),k).eq.1) iv = iv + 1
         if(system(mb(i+1),md(j+1),k).eq.1) iv = iv + 1
         if(system(mb(i-1),md(j+1),k).eq.1) iv = iv + 1         

      end if

      ivcoord2 = iv

      end


      subroutine getiens(system,isystem,isites,insites,
     &     ienergies,ietab,istab,intab,imenergies,imtab,emin,nstep,
     &     ienpos)
c----------------------------------------------------------
c     assigns energy barriers for spiro and grafted       -
c     interstitials to move based on the local atomic     -
c     environment                                         -
c     to move into vacant lattice positions based on the  -
c     local atomic environment                            -
c                                                         -
c     ienergies(nsite,1) - barrier to rotate clockwise    -
c     ienergies(nsite,2) - barrier to rotate anticlockwise-
c     ienergies(nsite,3) - barrier to move to nearest     -
c     neighbour alpha site                                -
c     if spiro:                                           -
c     ienergies(nsite,4) - barrier to move to up grafted  -
c     ienergies(nsite,5) - barrier to move to down grafted-
c     ienergies(nsite,6) - barrier to close intimate      -
c     frenkel pair                                        -
c     if grafted:                                         -
c     ienergies(nsite,4) - barrier to move to spiro       -
c     ienergies(nsite,5) - barrier to move to local       -
c     (alpha) split                                       -
c     ienergies(nsite,6) - barrier to move to non-local   -
c     (beta) split                                        -
c                                                         -
c     interstitial atom implies it sits on an occupied    -
c     lattice site                                        -
c                                                         -
c     we look at the plane of vacancies/splits above and  -
c     below the interstitial plane -                      -
c     whichever matches a library record first is used    -
c     for the energies - therefore place the closest I-V  -
c     interactions first in the iens.dat                  -
c----------------------------------------------------------            

      integer*1 system(800,800,199)
      integer*1 isystem(800,800,199)
      integer isites(9999,4),ienpos(9999)
      real*8  ienergies(9999,6),imenergies(9999,6)
      character*1 cmat(41,20),c2mat(41,20)
      logical sw1,up,emin

      integer*2 smat(21,30),umat(21,30),dmat(21,30),amat(21,30)
      integer*2 sumat(21,30)
      integer*2 istab(21,30,12,499),sdmat(21,30)
      real*8 ietab(6,499),imtab(6,499)
      integer*2 stab_i,stab_l

      do ns=1,insites
         

c--------initialise energies list
         do ne=1,6
            ienergies(ns,ne) = 9.0
            imenergies(ns,ne) = 9.0
         end do

c-----in-plane barriers

c.....map local system structure into temporary arrays
         
         ii = isites(ns,1)
         jj = isites(ns,2)
         kk = isites(ns,3)
         if(isites(ns,4).ge.11 .and. isites(ns,4).le.16) then
            ndir = isites(ns,4) - 10
         else if(isites(ns,4).ge.21 .and. isites(ns,4).le.26) then
            ndir = isites(ns,4) - 20
         else if(isites(ns,4).ge.1 .and. isites(ns,4).le.6) then
            ndir = isites(ns,4)
         else
            write(*,*)'error in getiens with isites(ns,4)'
            write(*,*)ns,ii,jj,kk,isites(ns,4)
            stop
         end if

         do i=1,21
            do j=1,30

               amat(i,j) = 0

c..............upper grafted interstitials or lower grafted
               sumat(i,j) = 0
               sdmat(i,j) = 0
               smat(i,j) = isystem(mb(ii-11+i),md(jj-16+j),kk)
c..............vacancy layer above the plane
               umat(i,j) = system(mb(ii-11+i),md(jj-16+j),ml(kk+1))
c..............vacancy layer below the plane
               dmat(i,j) = system(mb(ii-11+i),md(jj-16+j),kk)

c-----------changing system lables into library format lables

c...........interlayer bound atoms

               if(umat(i,j).eq.5 .or. umat(i,j).eq.6) then
                  umat(i,j) = 4
               end if

               if(dmat(i,j).eq.5 .or. dmat(i,j).eq.6) then
                  dmat(i,j) = 4
               end if  

               if(umat(i,j).eq.7) then
                  umat(i,j) = 1
               end if
               
               if(dmat(i,j).eq.7) then
                  dmat(i,j) = 1
               end if

               if(smat(i,j).ge.1 .and. smat(i,j).le.6) then
                  sdmat(i,j) = smat(i,j) + 10
                  sumat(i,j) = smat(i,j) + 10
                  smat(i,j) = smat(i,j) + 10
               else if(smat(i,j).ge.11 .and. smat(i,j).le.16) then
                  sumat(i,j) = 0
                  sdmat(i,j) = smat(i,j) + 10
                  smat(i,j) = smat(i,j) + 10
               else if(smat(i,j).ge.21 .and. smat(i,j).le.26) then
                  sumat(i,j) = smat(i,j)
                  sdmat(i,j) = 0
                  smat(i,j) = smat(i,j) + 10
               end if

               if(smat(i,j).eq.8) then
                  smat(i,j) = 0
                  amat(i,j) = -1
               end if

            end do
         end do

c         write(*,*)'site: ',ns
c         write(*,*)ii,jj,kk
         
c         do j=24,1,-1
c            write(*,'(17i1)') ( smat(i,j) , i=1,17 )
c         end do

c         do j=24,1,-1
c            write(*,'(17i1)') ( stab(i,j,3,1) , i=1,17 )
c         end do

c........check tables for structure
         
c........check below

         nchk = 0

         do m=1,intab
            do n=(ndir*2-1),(ndir*2)
               nchk = 0
               do i=1,21
                  do j=1,30

                     stab_l = istab(i,j,n,m)

c--------------------separate out lattice plane and interstitial plane
c                     in the library record
                     if(stab_l.gt.10) then
                        stab_i = stab_l
                        if(stab_i.eq.17) then
                           stab_l = 1
                        else if(stab_i.eq.18) then
                           stab_l = 1
                        else if(stab_i.ge.31) then
                           stab_l = 2
                           stab_i = stab_i - 20
                        else
                           stab_l = 1
                        end if
                     else if(stab_l.eq.3) then
                        stab_i = 17
                     else
                        stab_i = 0
                     end if

c                     if(nstep.eq.4801) then
c                        write(*,'(8i5)')m,n,i,j,stab_i,stab_l,
c     &                       sdmat(i,j),dmat(i,j)
c                     end if


c--------------------if the lattice plane matches and either of the 
c                    up or down interstitial planes

c                    if there is a different interstitial above and 
c                    below the plane, the first library match is used
                     if((sdmat(i,j).eq.stab_i .or. stab_i.eq.17  .or.
     &                    (stab_i.eq.18 .and. sdmat(i,j).ge.11 .and.
     &                    sdmat(i,j).le.16)) .and. 
     &                    (dmat(i,j).eq.stab_l .or. stab_l.eq.3)) then
                        nchk = nchk + 1
                     else
                        go to 101
                     end if
                  end do
               end do
               if(nchk.eq.630) then 
                  up = .false.
                  go to 102
               end if
 101           continue
            end do
         end do

c........check above

         nchk = 0

         do m=1,intab
            do n=(ndir*2-1),(ndir*2)
               nchk = 0
               do i=1,21
                  do j=1,30

                     stab_l = istab(i,j,n,m)

c--------------------separate out lattice plane and interstitial plane
c                     in the library record

                     if(stab_l.gt.10) then
                        stab_i = stab_l
                        if(stab_i.eq.17) then
                           stab_l = 1
                        else if(stab_i.eq.18) then
                           stab_l = 1
                        else if(stab_i.ge.31) then
                           stab_l = 2
                           stab_i = stab_i - 20
                        else
                           stab_l = 1
                        end if
                     else if(stab_l.eq.3) then
                        stab_i = 17
                     else
                        stab_i = 0
                     end if
                     
c--------------------if the lattice plane matches and either of the 
c                    up or down interstitial planes

c                    if there is a different interstitial above and 
c                    below the plane, the first library match is used
                     if((sumat(i,j).eq.stab_i .or. stab_i.eq.17  .or.
     &                    (stab_i.eq.18 .and. sumat(i,j).ge.11 .and.
     &                    sumat(i,j).le.16)) .and. 
     &                    (umat(i,j).eq.stab_l .or. stab_l.eq.3)) then
                        nchk = nchk + 1
                     else
                        go to 111
                     end if
                  end do
               end do
               if(nchk.eq.630) then
                  up = .true.
                  go to 102
               end if
 111           continue
            end do
         end do

         write(*,*)'unmatched structure - interstitial'

         write(*,*)'x ',ii,' y ',jj,' z ',kk 
         
         ic = 0
         jc = 0

         sw1 = .false.

c........writing unmatched structure to output

         write(*,*)'interstitial plane'

         ic = 0
         jc = 0

         do j=1,30
            ic = 0
            if(umat(1,j).ge.1 .or. umat(2,j).ge.1 .or. umat(3,j).ge.1
     &           .or. umat(4,j).ge.1 .or. umat(5,j).ge.1) then
               go to 302
            else 
               go to 303
            end if
 302        jc = jc + 1
            do i=1,21
               ic = ic + 1
               if(smat(i,j).eq.0) then
                  if(amat(i,j).eq.-1) then
                     c2mat(ic,jc) = '.'
                  else
                     c2mat(ic,jc) = ' '
                  end if
               else if(smat(i,j).eq.11) then
                  c2mat(ic,jc) = '1'
               else if(smat(i,j).eq.12) then
                  c2mat(ic,jc) = '2'
               else if(smat(i,j).eq.13) then
                  c2mat(ic,jc) = '3'
               else if(smat(i,j).eq.14) then
                  c2mat(ic,jc) = '4'
               else if(smat(i,j).eq.15) then
                  c2mat(ic,jc) = '5'
               else if(smat(i,j).eq.16) then
                  c2mat(ic,jc) = '6'
               else if(smat(i,j).eq.21) then
                  c2mat(ic,jc) = 'a'
               else if(smat(i,j).eq.22) then
                  c2mat(ic,jc) = 'b'
               else if(smat(i,j).eq.23) then
                  c2mat(ic,jc) = 'c'
               else if(smat(i,j).eq.24) then
                  c2mat(ic,jc) = 'd'
               else if(smat(i,j).eq.25) then
                  c2mat(ic,jc) = 'e'
               else if(smat(i,j).eq.26) then
                  c2mat(ic,jc) = 'f'
               else if(smat(i,j).eq.31) then
                  c2mat(ic,jc) = 'A'
               else if(smat(i,j).eq.32) then
                  c2mat(ic,jc) = 'B'
               else if(smat(i,j).eq.33) then
                  c2mat(ic,jc) = 'C'
               else if(smat(i,j).eq.34) then
                  c2mat(ic,jc) = 'D'
               else if(smat(i,j).eq.35) then
                  c2mat(ic,jc) = 'E'
               else if(smat(i,j).eq.36) then
                  c2mat(ic,jc) = 'F'
               else
                  write(*,*)'error in inter plane: unrecognised atom'
                  write(*,*)smat(i,j)
                  stop
               end if
               if(ic.eq.41) go to 301
               ic = ic + 1
               c2mat(ic,jc) = ' '
            end do
 301        continue
 303        continue
         end do

         do j=20,1,-1
            write(*,'(41a1)') ( c2mat(i,j) , i=1,41 )
         end do

         write(*,*)'upper lattice plane'

         ic = 0
         jc = 0

         do j=1,30
            ic = 0
            if(umat(1,j).ge.1 .or. umat(2,j).ge.1 .or. umat(3,j).ge.1
     &           .or. umat(4,j).ge.1 .or. umat(5,j).ge.1) then
               go to 202
            else 
               go to 203
            end if
 202        jc = jc + 1
            do i=1,21
               ic = ic + 1
               if(umat(i,j).eq.0) then
                  cmat(ic,jc) = ' '
               else if(umat(i,j).eq.1) then
                  cmat(ic,jc) = 'o'
               else if(umat(i,j).eq.2) then
                  cmat(ic,jc) = '-'
               else if(umat(i,j).eq.3) then
                  cmat(ic,jc) = '*'
               else if(umat(i,j).eq.4) then
                  cmat(ic,jc) = '+'
               else if(umat(i,j).eq.9) then
                  cmat(ic,jc) = 's'
               end if
               if(ic.eq.41) go to 201
               ic = ic + 1
               cmat(ic,jc) = ' '
            end do
 201        continue
 203        continue
         end do

         do j=20,1,-1
            write(*,'(41a1)') ( cmat(i,j) , i=1,41 )
         end do

         write(*,*)'lower lattice plane'

         ic = 0
         jc = 0

         do j=1,30
            ic = 0
            if(dmat(1,j).ge.1 .or. dmat(2,j).ge.1 .or. dmat(3,j).ge.1
     &           .or. dmat(4,j).ge.1 .or. dmat(5,j).ge.1) then
               go to 402
            else 
               go to 403
            end if
 402        jc = jc + 1
            do i=1,21
               ic = ic + 1
               if(dmat(i,j).eq.0) then
                  cmat(ic,jc) = ' '
               else if(dmat(i,j).eq.1) then
                  cmat(ic,jc) = 'o'
               else if(dmat(i,j).eq.2) then
                  cmat(ic,jc) = '-'
               else if(dmat(i,j).eq.3) then
                  cmat(ic,jc) = '*'
               else if(dmat(i,j).eq.4) then
                  cmat(ic,jc) = '+'
               else if(dmat(i,j).eq.9) then
                  cmat(ic,jc) = 's'
               end if
               if(ic.eq.41) go to 401
               ic = ic + 1
               cmat(ic,jc) = ' '
            end do
 401        continue
 403        continue
         end do

         do j=20,1,-1
            write(*,'(41a1)') ( cmat(i,j) , i=1,41 )
         end do

         stop

c--------if spiro
 102     if(isites(ns,4).le.6) then
            
            ienergies(ns,1) = ietab(1,m)
            ienergies(ns,2) = ietab(2,m)
            ienergies(ns,3) = ietab(3,m)
            ienergies(ns,6) = ietab(6,m)

            imenergies(ns,1) = imtab(1,m)
            imenergies(ns,2) = imtab(2,m)
            imenergies(ns,3) = imtab(3,m)
            imenergies(ns,6) = imtab(6,m)

            if(up) then
               ienergies(ns,4) = ietab(4,m)
               ienergies(ns,5) = ietab(5,m)
               
               imenergies(ns,4) = imtab(4,m)
               imenergies(ns,5) = imtab(5,m)
            else
               ienergies(ns,5) = ietab(4,m)
               ienergies(ns,4) = ietab(5,m)
               
               imenergies(ns,5) = imtab(4,m)
               imenergies(ns,4) = imtab(5,m)
            end if
            
c---------if grafted
         else if(isites(ns,4).ge.11 .and. 
     &           isites(ns,4).le.26) then

            ienergies(ns,1) = ietab(1,m)
            ienergies(ns,2) = ietab(2,m)
            ienergies(ns,3) = ietab(3,m)
            ienergies(ns,4) = ietab(4,m)
            ienergies(ns,5) = ietab(5,m)
            ienergies(ns,6) = ietab(6,m)

            imenergies(ns,1) = imtab(1,m)
            imenergies(ns,2) = imtab(2,m)
            imenergies(ns,3) = imtab(3,m)
            imenergies(ns,4) = imtab(4,m)
            imenergies(ns,5) = imtab(5,m)
            imenergies(ns,6) = imtab(6,m)
         else
            write(*,*)'error in geti - isites'
            stop
         end if

         ienpos(ns) = m

      end do

      end




      subroutine getsens(system,isystem,ssites,snsites,
     &     senergies,ietab,istab,intab,smenergies,imtab,emin,senpos)
c----------------------------------------------------------
c     assigns energy barriers for split interstitial sites-
c     to move to neighbouring grafted interstitial sites  -
c     based on the local atomic environment               -
c                                                         -
c     senergies(nsite,1) - barrier to move up | or /      -
c     senergies(nsite,2) - barrier to move up \ or |      -
c     senergies(nsite,3) - barrier to move up / or \      -
c     senergies(nsite,4) - barrier to move down | or /    -
c     senergies(nsite,5) - barrier to move down \ or |    -
c     senergies(nsite,6) - barrier to move down / or \    -
c----------------------------------------------------------         

      integer*1 system(800,800,199)
      integer*1 isystem(800,800,199)
      integer ssites(9999,5),senpos(9999)
      real*8  senergies(9999,6), smenergies(9999,6)
      character*1 cmat(41,20),c2mat(41,20)
      logical sw1,up,emin
      integer snsites

      integer*2 smat(21,30),umat(21,30),dmat(21,30),amat(21,30)
      integer*2 istab(21,30,12,499),admat(21,30)
      real*8 ietab(6,499),imtab(6,499)
      integer*2 stab_i,stab_l

      do ns=1,snsites

c--------initialise energies list
         do ne=1,6
            senergies(ns,ne) = 9.0
         end do

c-----in-plane barriers

c.....map local system structure into temporary arrays
         
         ii = ssites(ns,1)
         jj = ssites(ns,2)
         kk = ssites(ns,3)

c........ndir = 1 is up:    |
c                          / \
c
c........ndir = 0 is down  \ /
c                           |

         ndir = ssites(ns,4)

         do i=1,21
            do j=1,30
               smat(i,j) = system(mb(ii-11+i),md(jj-16+j),kk)
c..............interstitial layer above the plane
               umat(i,j) = isystem(mb(ii-11+i),md(jj-16+j),kk)
c..............interstitial layer below the plane
               dmat(i,j) = isystem(mb(ii-11+i),md(jj-16+j),ml(kk-1))

c-----------changing system lables into library format lables

c...........interlayer bound atoms

               if(smat(i,j).eq.5 .or. smat(i,j).eq.6) then
                  smat(i,j) = 4
               end if

c               if(smat(i,j).eq.7) then
c                  smat(i,j) = 1
c               end if

               amat(i,j) = 0
               admat(i,j) = 0

c...........interstitials
               if(umat(i,j).eq.8) then
                  umat(i,j) = 0
                  amat(i,j) = -1
               else if(umat(i,j).gt.20) then
                  umat(i,j) = 0
               else if(umat(i,j).ge.11 .and. umat(i,j).le.16) then
                  umat(i,j) = umat(i,j) + 10
               else if(umat(i,j).ge.1 .and. umat(i,j).le.6) then
                  umat(i,j) = umat(i,j) + 10
               end if

c...........interstitials
               if(dmat(i,j).eq.8) then
                  dmat(i,j) = 0
                  admat(i,j) = -1
               else if(dmat(i,j).ge.11 .and. dmat(i,j).le.16) then
                  dmat(i,j) = 0
               else if(dmat(i,j).ge.1 .and. dmat(i,j).le.6) then
                  dmat(i,j) = dmat(i,j) + 10
               end if

            end do
         end do

c         write(*,*)'site: ',ns
c         write(*,*)ii,jj,kk
         
c         do j=24,1,-1
c            write(*,'(17i1)') ( smat(i,j) , i=1,17 )
c         end do

c         do j=24,1,-1
c            write(*,'(17i1)') ( stab(i,j,3,1) , i=1,17 )
c         end do

c........check tables for structure
         
c--------check above

         nchk = 0

         do m=1,intab
            do n=1,12
               nchk = 0
               do i=1,21
                  do j=1,30

                     stab_l = istab(i,j,n,m)

c--------------------separate out lattice plane and interstitial plane
c                     in the library record
                     if(stab_l.gt.10) then
                        stab_i = stab_l
                        if(stab_i.eq.17) then
                           stab_l = 3
                        else if(stab_i.ge.31) then
                           stab_l = 2
                           stab_i = stab_i - 20
                        else
                           stab_l = 1
                        end if
                     else if(stab_l.eq.3) then
                        stab_i = 17
                     else
                        stab_i = 0
                     end if

c--------------------if the lattice plane matches and either of the 
c                    up or down interstitial planes

c                    if there is a different interstitial above and 
c                    below the plane, the first library match is used
                     if((smat(i,j).eq.stab_l .or. stab_l.eq.3) .and. 
     &                    (umat(i,j).eq.stab_i .or. stab_i.eq.17)) then
                        nchk = nchk + 1
                     else
                        go to 101
                     end if
                  end do
               end do
               if(nchk.eq.630) then
                  up = .true.
                  go to 102
               end if
 101           continue
            end do
         end do


c--------check below

         nchk = 0

         do m=1,intab
            do n=1,12
               nchk = 0
               do i=1,21
                  do j=1,30

                     stab_l = istab(i,j,n,m)

c--------------------separate out lattice plane and interstitial plane
c                     in the library record
                     if(stab_l.gt.10) then
                        stab_i = stab_l
                        if(stab_i.eq.17) then
                           stab_l = 3
                        else
                           stab_l = 1
                        end if
                     else
                        stab_i = 0
                     end if

c--------------------if the lattice plane matches and either of the 
c                    up or down interstitial planes

c                    if there is a different interstitial above and 
c                    below the plane, the first library match is used
                     if((smat(i,j).eq.stab_l .or. stab_l.eq.3) .and. 
     &                    (dmat(i,j).eq.stab_i .or. stab_i.eq.17)) then
                        nchk = nchk + 1
                     else
                        go to 111
                     end if
                  end do
               end do
               if(nchk.eq.630) then
                  up = .false.
                  go to 102
               end if
 111           continue
            end do
         end do

c........if no match - write unmatched structure and stop

         write(*,*)'unmatched structure - gets'

         ic = 0
         jc = 0

         sw1 = .false.

c........writing unmatched structure to output

         write(*,*)'lattice plane'

         do j=1,30
            ic = 0
            if(smat(1,j).ge.1 .or. smat(2,j).ge.1 .or. smat(3,j).ge.1
     &           .or. smat(4,j).ge.1 .or. smat(5,j).ge.1) then
               go to 202
            else 
               go to 203
            end if
 202        jc = jc + 1
            do i=1,21
               ic = ic + 1
               if(smat(i,j).eq.0) then
                  cmat(ic,jc) = ' '
               else if(smat(i,j).eq.1) then
                  cmat(ic,jc) = 'o'
               else if(smat(i,j).eq.2) then
                  cmat(ic,jc) = '-'
               else if(smat(i,j).eq.3) then
                  cmat(ic,jc) = '*'
               else if(smat(i,j).eq.4) then
                  cmat(ic,jc) = '+'
               else if(smat(i,j).eq.9) then
                  cmat(ic,jc) = 's'
               end if
               if(ic.eq.41) go to 201
               ic = ic + 1
               cmat(ic,jc) = ' '
            end do
 201        continue
 203        continue
         end do

         do j=20,1,-1
            write(*,'(41a1)') ( cmat(i,j) , i=1,41 )
         end do

         write(*,*)'lower interstitial plane'

         ic = 0
         jc = 0

         do j=1,30
            ic = 0
            if(smat(1,j).ge.1 .or. smat(2,j).ge.1 .or. smat(3,j).ge.1
     &           .or. smat(4,j).ge.1 .or. smat(5,j).ge.1) then
               go to 302
            else 
               go to 303
            end if
 302        jc = jc + 1
            do i=1,21
               ic = ic + 1
               if(dmat(i,j).eq.0) then
                  if(admat(i,j).eq.-1) then
                     c2mat(ic,jc) = '.'
                  else
                     c2mat(ic,jc) = ' '
                  end if
               else if(dmat(i,j).eq.11) then
                  c2mat(ic,jc) = '1'
               else if(dmat(i,j).eq.12) then
                  c2mat(ic,jc) = '2'
               else if(dmat(i,j).eq.13) then
                  c2mat(ic,jc) = '3'
               else if(dmat(i,j).eq.14) then
                  c2mat(ic,jc) = '4'
               else if(dmat(i,j).eq.15) then
                  c2mat(ic,jc) = '5'
               else if(dmat(i,j).eq.16) then
                  c2mat(ic,jc) = '6'
               else if(dmat(i,j).eq.21) then
                  c2mat(ic,jc) = 'a'
               else if(dmat(i,j).eq.22) then
                  c2mat(ic,jc) = 'b'
               else if(dmat(i,j).eq.23) then
                  c2mat(ic,jc) = 'c'
               else if(dmat(i,j).eq.24) then
                  c2mat(ic,jc) = 'd'
               else if(dmat(i,j).eq.25) then
                  c2mat(ic,jc) = 'e'
               else if(dmat(i,j).eq.26) then
                  c2mat(ic,jc) = 'f'
               else
                  write(*,*)'error in sinter plane: unrecognised atom'
                  write(*,*)smat(i,j)
                  stop
               end if
               if(ic.eq.41) go to 301
               ic = ic + 1
               c2mat(ic,jc) = ' '
            end do
 301        continue
 303        continue
         end do

         do j=20,1,-1
            write(*,'(41a1)') ( c2mat(i,j) , i=1,41 )
         end do

         write(*,*)'upper interstitial plane'

         ic = 0
         jc = 0

         do j=1,30
            ic = 0
            if(smat(1,j).ge.1 .or. smat(2,j).ge.1 .or. smat(3,j).ge.1
     &           .or. smat(4,j).ge.1 .or. smat(5,j).ge.1) then
               go to 402
            else 
               go to 403
            end if
 402        jc = jc + 1
            do i=1,21
               ic = ic + 1
               if(umat(i,j).eq.0) then
                  if(amat(i,j).eq.-1) then
                     c2mat(ic,jc) = '.'
                  else
                     c2mat(ic,jc) = ' '
                  end if
               else if(umat(i,j).eq.11) then
                  c2mat(ic,jc) = '1'
               else if(umat(i,j).eq.12) then
                  c2mat(ic,jc) = '2'
               else if(umat(i,j).eq.13) then
                  c2mat(ic,jc) = '3'
               else if(umat(i,j).eq.14) then
                  c2mat(ic,jc) = '4'
               else if(umat(i,j).eq.15) then
                  c2mat(ic,jc) = '5'
               else if(umat(i,j).eq.16) then
                  c2mat(ic,jc) = '6'
               else if(umat(i,j).eq.21) then
                  c2mat(ic,jc) = 'a'
               else if(umat(i,j).eq.22) then
                  c2mat(ic,jc) = 'b'
               else if(umat(i,j).eq.23) then
                  c2mat(ic,jc) = 'c'
               else if(umat(i,j).eq.24) then
                  c2mat(ic,jc) = 'd'
               else if(umat(i,j).eq.25) then
                  c2mat(ic,jc) = 'e'
               else if(umat(i,j).eq.26) then
                  c2mat(ic,jc) = 'f'
               else
                  write(*,*)'error in sinter plane: unrecognised atom'
                  write(*,*)umat(i,j)
                  stop
               end if
               if(ic.eq.41) go to 401
               ic = ic + 1
               c2mat(ic,jc) = ' '
            end do
 401        continue
 403        continue
         end do

         do j=20,1,-1
            write(*,'(41a1)') ( c2mat(i,j) , i=1,41 )
         end do

         stop

 102     if(n.ge.1 .or. n.le.4) then
            if(up) then
               senergies(ns,1) = ietab(1,m)
               senergies(ns,2) = ietab(2,m)
               senergies(ns,3) = ietab(3,m)
               senergies(ns,4) = ietab(4,m)
               senergies(ns,5) = ietab(5,m)
               senergies(ns,6) = ietab(6,m)

               smenergies(ns,1) = imtab(1,m)
               smenergies(ns,2) = imtab(2,m)
               smenergies(ns,3) = imtab(3,m)
               smenergies(ns,4) = imtab(4,m)
               smenergies(ns,5) = imtab(5,m)
               smenergies(ns,6) = imtab(6,m)
            else
               senergies(ns,1) = ietab(4,m)
               senergies(ns,2) = ietab(5,m)
               senergies(ns,3) = ietab(6,m)
               senergies(ns,4) = ietab(1,m)
               senergies(ns,5) = ietab(2,m)
               senergies(ns,6) = ietab(3,m)

               smenergies(ns,1) = imtab(4,m)
               smenergies(ns,2) = imtab(5,m)
               smenergies(ns,3) = imtab(6,m)
               smenergies(ns,4) = imtab(1,m)
               smenergies(ns,5) = imtab(2,m)
               smenergies(ns,6) = imtab(3,m)
            end if
         else if(n.ge.5 .or. n.le.8) then
            if(up) then
               senergies(ns,2) = ietab(1,m)
               senergies(ns,3) = ietab(2,m)
               senergies(ns,1) = ietab(3,m)
               senergies(ns,5) = ietab(4,m)
               senergies(ns,6) = ietab(5,m)
               senergies(ns,4) = ietab(6,m)
               
               smenergies(ns,2) = imtab(1,m)
               smenergies(ns,3) = imtab(2,m)
               smenergies(ns,1) = imtab(3,m)
               smenergies(ns,5) = imtab(4,m)
               smenergies(ns,6) = imtab(5,m)
               smenergies(ns,4) = imtab(6,m)
            else
               senergies(ns,2) = ietab(4,m)
               senergies(ns,3) = ietab(5,m)
               senergies(ns,1) = ietab(6,m)
               senergies(ns,5) = ietab(1,m)
               senergies(ns,6) = ietab(2,m)
               senergies(ns,4) = ietab(3,m)

               smenergies(ns,2) = imtab(4,m)
               smenergies(ns,3) = imtab(5,m)
               smenergies(ns,1) = imtab(6,m)
               smenergies(ns,5) = imtab(1,m)
               smenergies(ns,6) = imtab(2,m)
               smenergies(ns,4) = imtab(3,m)
            end if
         else if(n.ge.9 .or. n.le.12) then
            if(up) then
               senergies(ns,3) = ietab(1,m)
               senergies(ns,1) = ietab(2,m)
               senergies(ns,2) = ietab(3,m)
               senergies(ns,6) = ietab(4,m)
               senergies(ns,4) = ietab(5,m)
               senergies(ns,5) = ietab(6,m)

               smenergies(ns,3) = imtab(1,m)
               smenergies(ns,1) = imtab(2,m)
               smenergies(ns,2) = imtab(3,m)
               smenergies(ns,6) = imtab(4,m)
               smenergies(ns,4) = imtab(5,m)
               smenergies(ns,5) = imtab(6,m)
            else
               senergies(ns,3) = ietab(4,m)
               senergies(ns,1) = ietab(5,m)
               senergies(ns,2) = ietab(6,m)
               senergies(ns,6) = ietab(1,m)
               senergies(ns,4) = ietab(2,m)
               senergies(ns,5) = ietab(3,m)
               
               smenergies(ns,3) = imtab(4,m)
               smenergies(ns,1) = imtab(5,m)
               smenergies(ns,2) = imtab(6,m)
               smenergies(ns,6) = imtab(1,m)
               smenergies(ns,4) = imtab(2,m)
               smenergies(ns,5) = imtab(3,m)
            end if
         end if
      end do

      senpos(ns) = m

      end

      subroutine getsites_lc(system,sites,nsites,nl,ixd,iyd,move)
c----------------------------------------------------------
c     updates the sites array for undercoordinated atoms  -
c     in the vicity of the last transition                -
c                                                         -
c                                                         -
c                      6  1  2                            -
c                       \ | /                             -
c                         C                               -
c                       / | \                             -
c                      5  4  3                            -
c                                                         -
c----------------------------------------------------------

      integer*1 system(800,800,199)
      integer sites(9999,5)
      integer sites_tmp(9999,5)
      integer move(4)
      logical up

c.....first remove all sites in the vicinity of move in the sites array

      nsites_tmp = 0
      nrem = 0

      izd = nl

      do i=1,nsites
         ixm = 0
         iym = 0
         izm = 0
         if(sites(i,1).le.8 .and. move(1).ge.(ixd-7)) ixm = -1
         if(sites(i,1).ge.(ixd-7) .and. move(1).le.8) ixm = 1
         if(sites(i,2).le.8 .and. move(2).ge.(iyd-7)) iym = -1
         if(sites(i,2).ge.(iyd-7) .and. move(2).le.8) iym = 1
         if(sites(i,3).le.2 .and. move(3).ge.(izd-1)) izm = -1
         if(sites(i,3).ge.(izd-1) .and. move(3).le.2) izm = 1
         if(sites(i,1).ge.(move(1)-4+ixm*ixd) .and. 
     &        sites(i,1).le.(move(1)+4+ixm*ixd) .and. 
     &        sites(i,2).ge.(move(2)-4+iym*iyd) .and. 
     &        sites(i,2).le.(move(2)+4+iym*iyd) .and. 
     &        sites(i,3).ge.(move(3)-1+izm*izd) .and. 
     &        sites(i,3).le.(move(3)+1+izm*izd)) then 
            nrem = nrem + 1
         else
            nsites_tmp = nsites_tmp + 1
            sites_tmp(nsites_tmp,1) = sites(i,1)
            sites_tmp(nsites_tmp,2) = sites(i,2)
            sites_tmp(nsites_tmp,3) = sites(i,3)
            sites_tmp(nsites_tmp,4) = sites(i,4)
            sites_tmp(nsites_tmp,5) = sites(i,5)
c            write(*,*)'dir',nsites_tmp,sites(i,4),sites(i,5)
         end if
      end do


      n_s = nsites_tmp

      if(nl.eq.1) then
         klp = 0
      else
         klp = 1
      end if

c-----loop over the vicinity
      do k=move(3)-klp,move(3)+klp
         do i=move(1)-4,move(1)+4
            do j=move(2)-4,move(2)+4

               if(system(mb(i),md(j),ml(k)).eq.1 .or. 
     &              system(mb(i),md(j),ml(k)).eq.5
     &              .or. system(mb(i),md(j),ml(k)).eq.6) then

c-----up or down

                  if(system(mb(i),md(j+2),ml(k)).ne.0) then
                     up = .true.
                  else
                     up = .false.
                  end if
                  
                  if(n_s.gt.9999) then
                     write(*,*)'error: increase size of sites array'
                     stop
                  end if

c-----if up, check for nearest neighbour vacancies

                  if(up) then

                     if(system(mb(i),md(j+2),ml(k)).ne.2  
     &                  .and. system(mb(i-1),md(j-1),ml(k)).ne.2  .and. 
     &                    system(mb(i+1),md(j-1),ml(k)).ne.2) then
                        go to 100
                     else if(system(mb(i),md(j+2),ml(k)).eq.2 .and.
     &                       system(mb(i-1),md(j-1),ml(k)).ne.2 .and.
     &                       system(mb(i+1),md(j-1),ml(k)).ne.2) then
                        n_s = n_s + 1
                        sites_tmp(n_s,1) = mb(i)
                        sites_tmp(n_s,2) = md(j)
                        sites_tmp(n_s,3) = ml(k)
                        sites_tmp(n_s,4) = 1
                  else if(system(mb(i),md(j+2),ml(k)).ne.2 .and.
     &                    system(mb(i-1),md(j-1),ml(k)).ne.2 .and.
     &                    system(mb(i+1),md(j-1),ml(k)).eq.2) then
                     n_s = n_s + 1
                     sites_tmp(n_s,1) = mb(i)
                     sites_tmp(n_s,2) = md(j)
                     sites_tmp(n_s,3) = ml(k)
                     sites_tmp(n_s,4) = 3
                  else if(system(mb(i),md(j+2),ml(k)).ne.2 .and.
     &                    system(mb(i-1),md(j-1),ml(k)).eq.2 .and.
     &                    system(mb(i+1),md(j-1),ml(k)).ne.2) then
                     n_s = n_s + 1
                     sites_tmp(n_s,1) = mb(i)
                     sites_tmp(n_s,2) = md(j)
                     sites_tmp(n_s,3) = ml(k)
                     sites_tmp(n_s,4) = 5
                  else
                     write(*,*)'Error in site - undercoordinated atom u'
                     write(*,*)system(mb(i),md(j+2),ml(k))
                     write(*,*)system(mb(i+1),md(j-1),ml(k))
                     write(*,*)system(mb(i-1),md(j-1),ml(k))
                     write(*,*)i,j,k
                     stop
                  end if

               else

                  if(system(mb(i),md(j-2),ml(k)).ne.2 .and. 
     &                 system(mb(i-1),md(j+1),ml(k)).ne.2 .and. 
     &                 system(mb(i+1),md(j+1),ml(k)).ne.2) then
                     go to 100
                  else if(system(mb(i),md(j-2),ml(k)).eq.2 .and.
     &                    system(mb(i-1),md(j+1),ml(k)).ne.2 .and.
     &                    system(mb(i+1),md(j+1),ml(k)).ne.2) then
                     n_s = n_s + 1
                     sites_tmp(n_s,1) = mb(i)
                     sites_tmp(n_s,2) = md(j)
                     sites_tmp(n_s,3) = ml(k)
                     sites_tmp(n_s,4) = 4
                  else if(system(mb(i),md(j-2),ml(k)).ne.2 .and.
     &                    system(mb(i-1),md(j+1),ml(k)).ne.2 .and.
     &                    system(mb(i+1),md(j+1),ml(k)).eq.2) then
                     n_s = n_s + 1
                     sites_tmp(n_s,1) = mb(i)
                     sites_tmp(n_s,2) = md(j)
                     sites_tmp(n_s,3) = ml(k)
                     sites_tmp(n_s,4) = 2
                 else if(system(mb(i),md(j-2),ml(k)).ne.2 .and.
     &                    system(mb(i-1),md(j+1),ml(k)).eq.2 .and.
     &                    system(mb(i+1),md(j+1),ml(k)).ne.2) then
                     n_s = n_s + 1
                     sites_tmp(n_s,1) = mb(i)
                     sites_tmp(n_s,2) = md(j)
                     sites_tmp(n_s,3) = ml(k)
                     sites_tmp(n_s,4) = 6
                  else
                     write(*,*)'Error in site - undercoordinated atom d'
                     write(*,*)system(mb(i),md(j-2),ml(k))
                     write(*,*)system(mb(i+1),md(j+1),ml(k))
                     write(*,*)system(mb(i-1),md(j+1),ml(k))
                     stop
                  end if

               end if
               end if
 100           continue
            end do
         end do
      end do

      do i=1,n_s
         sites(i,1) = sites_tmp(i,1)
         sites(i,2) = sites_tmp(i,2)
         sites(i,3) = sites_tmp(i,3)
         sites(i,4) = sites_tmp(i,4)
         sites(i,5) = sites_tmp(i,5)
      end do
         
      nsites = n_s

      end 

      subroutine getssties_lc(system,ssites,snsites,nl,ixd,iyd,move)
c----------------------------------------------------------------
c     updates the ssites array in the vicinity of the last 
c     transition
c----------------------------------------------------------------
      integer*1 system(800,800,199)
      integer ssites(9999,5)
      integer ssites_tmp(9999,5)
      integer snsites
      integer move(4)

c.....first remove all sites in the vicinity of move in the ssites array

      nsites_tmp = 0
      nrem = 0

      izd = nl

      do i=1,snsites
         ixm = 0
         iym = 0
         izm = 0
         if(ssites(i,1).le.8 .and. move(1).ge.(ixd-7)) ixm = -1
         if(ssites(i,1).ge.(ixd-7) .and. move(1).le.8) ixm = 1
         if(ssites(i,2).le.8 .and. move(2).ge.(iyd-7)) iym = -1
         if(ssites(i,2).ge.(iyd-7) .and. move(2).le.8) iym = 1
         if(ssites(i,3).le.2 .and. move(3).ge.(izd-1)) izm = -1
         if(ssites(i,3).ge.(izd-1) .and. move(3).le.2) izm = 1
         if(ssites(i,1).ge.(move(1)-4+ixm*ixd) .and. 
     &        ssites(i,1).le.(move(1)+4+ixm*ixd) .and. 
     &        ssites(i,2).ge.(move(2)-4+iym*iyd) .and. 
     &        ssites(i,2).le.(move(2)+4+iym*iyd) .and. 
     &        ssites(i,3).ge.(move(3)-1+izm*izd) .and. 
     &        ssites(i,3).le.(move(3)+1+izm*izd)) then 
            nrem = nrem + 1
         else
            nsites_tmp = nsites_tmp + 1
            ssites_tmp(nsites_tmp,1) = ssites(i,1)
            ssites_tmp(nsites_tmp,2) = ssites(i,2)
            ssites_tmp(nsites_tmp,3) = ssites(i,3)
            ssites_tmp(nsites_tmp,4) = ssites(i,4)
            ssites_tmp(nsites_tmp,5) = ssites(i,5)
         end if
      end do
      
c      write(*,*)nrem,nsites_tmp

      n_s = nsites_tmp

c-----loop over the vicinity
      do k=move(3)-1,move(3)+1
         do i=move(1)-4,move(1)+4
            do j=move(2)-4,move(2)+4
                              
               if(n_s.gt.9999) then
                  write(*,*)'error: increase size of ssites array'
                  stop
               end if

               if(system(mb(i),md(j),ml(k)).eq.9) then
                  n_s = n_s + 1
                  ssites_tmp(n_s,1) = mb(i)
                  ssites_tmp(n_s,2) = md(j)
                  ssites_tmp(n_s,3) = ml(k)
                  if(system(mb(i),md(j+2),ml(k)).ne.0) then
                     ssites_tmp(n_s,4) = 1
                  else
                     ssites_tmp(n_s,4) = 0
                  end if               
               end if

            end do
         end do
      end do

      do i=1,n_s
         ssites(i,1) = ssites_tmp(i,1)
         ssites(i,2) = ssites_tmp(i,2)
         ssites(i,3) = ssites_tmp(i,3)
         ssites(i,4) = ssites_tmp(i,4)
         ssites(i,5) = ssites_tmp(i,5)
      end do

      snsites = n_s

c      write(*,*)snsites

      end

      subroutine getisites_lc(isystem,isites,insites,nl,ixd,iyd,move)
c----------------------------------------------------------------
c     updates the isites array for interstitial sites in the 
c     vicinity of the last transition
c----------------------------------------------------------------
      integer*1 isystem(800,800,199)
      integer isites(9999,5)
      integer isites_tmp(9999,5)
      integer move(4)

c.....first remove all sites in the vicinity of move in the ssites array

      nsites_tmp = 0
      nrem = 0

      izd = nl

      do i=1,insites
         ixm = 0
         iym = 0
         izm = 0
         if(isites(i,1).le.8 .and. move(1).ge.(ixd-7)) ixm = -1
         if(isites(i,1).ge.(ixd-7) .and. move(1).le.8) ixm = 1
         if(isites(i,2).le.8 .and. move(2).ge.(iyd-7)) iym = -1
         if(isites(i,2).ge.(iyd-7) .and. move(2).le.8) iym = 1
         if(isites(i,3).le.2 .and. move(3).ge.(izd-1)) izm = -1
         if(isites(i,3).ge.(izd-1) .and. move(3).le.2) izm = 1
         if(isites(i,1).ge.(move(1)-4+ixm*ixd) .and. 
     &        isites(i,1).le.(move(1)+4+ixm*ixd) .and. 
     &        isites(i,2).ge.(move(2)-4+iym*iyd) .and. 
     &        isites(i,2).le.(move(2)+4+iym*iyd) .and. 
     &        isites(i,3).ge.(move(3)-1+izm*izd) .and. 
     &        isites(i,3).le.(move(3)+1+izm*izd)) then 
            nrem = nrem + 1
         else
            nsites_tmp = nsites_tmp + 1
            isites_tmp(nsites_tmp,1) = isites(i,1)
            isites_tmp(nsites_tmp,2) = isites(i,2)
            isites_tmp(nsites_tmp,3) = isites(i,3)
            isites_tmp(nsites_tmp,4) = isites(i,4)
            isites_tmp(nsites_tmp,5) = isites(i,5)
         end if
      end do
      
c      write(*,*)nrem,nsites_tmp

      n_s = nsites_tmp

c-----loop over the vicinity
      do k=move(3)-1,move(3)+1
         do i=move(1)-4,move(1)+4
            do j=move(2)-4,move(2)+4

               if(n_s.gt.9999) then
                  write(*,*)'error: increase size of isites array'
                  stop
               end if

               if(isystem(mb(i),md(j),ml(k)).ge.1 .and. 
     &              isystem(mb(i),md(j),ml(k)).le.6) then
                  n_s = n_s + 1
                  isites_tmp(n_s,1) = mb(i)
                  isites_tmp(n_s,2) = md(j)
                  isites_tmp(n_s,3) = ml(k)
                  isites_tmp(n_s,4) = isystem(mb(i),md(j),ml(k))
               else if(isystem(mb(i),md(j),ml(k)).ge.11) then
                  n_s = n_s + 1
                  isites_tmp(n_s,1) = mb(i)
                  isites_tmp(n_s,2) = md(j)
                  isites_tmp(n_s,3) = ml(k)
                  isites_tmp(n_s,4) = isystem(mb(i),md(j),ml(k))
               end if
            end do
         end do
      end do

      do i=1,n_s
         isites(i,1) = isites_tmp(i,1)
         isites(i,2) = isites_tmp(i,2)
         isites(i,3) = isites_tmp(i,3)
         isites(i,4) = isites_tmp(i,4)
      end do

      insites = n_s

c      write(*,*)insites

      end
