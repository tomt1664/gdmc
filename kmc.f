      Program main
c-------------------------------------------------------------
c       ____                 _     ____  __  __  ____        -
c      / ___|_ __ __ _ _ __ | |__ |  _ \|  \/  |/ ___|       -
c     | |  _| '__/ _` | '_ \| '_ \| | | | |\/| | |           -
c     | |_| | | | (_| | |_) | | | | |_| | |  | | |___        -
c      \____|_|  \__,_| .__/|_| |_|____/|_|  |_|\____|       -
c                                                            -
c     Kinetic Monte Carlo code for simulating point defect   -
c     dynamics and aggregation in graphite/graphene          -
c                                                            -
c     Rates are determined from explicit ab-initio derived   -
c     activation barriers, supplied in the files ens.dat     -
c     and iens.dat                                           -
c                                                            -
c     Version 1.3.0                                          -
c                                                            -
c     T. Trevethan 2015                                      -
c-------------------------------------------------------------

c     the system: 
c     C atom - 1
c     vacancy - 2
c     interlayer bond up - 5
c     interlayer bond down - 6
c     split interstitial - 9
c
c     the latice system is represented as a 3 dimensional integer
c     grid with i,j,k axes

      integer*1 system(800,800,199),atm,atm2,atm3,atm4

c     the system dimensions in array units
c     one orthorhomic unit cell is 4x6 array units

      integer ixd,iyd,izd

c     periodic boundary conditions - in x,y and z

      logical per_x,per_y,per_z

c     share globally to be available to the periodic functions

      COMMON ixd,iyd,izd,per_x,per_y,per_z

c     spiro and grafted interstitials are represented as an 
c     interweiving grid
c
c     alpha site: 8
c     spiro interstitial: 1-6
c     lower grafted interstitial: 11-16
c     upper grafted interstitial: 21-26

      integer*1 isystem(800,800,199)

c     list of active sites
c     sites(nsite,1-3) - i,j,k indices of the site
c     sites(nsite,4) - dangling bond direction (1-6)
c     sites(nsite,5) - dangling bond saturation (0 - no, 1 - yes)
c
c     interstitials:
c     isites(insites,1-3) - i,j,k indices of the site
c     isites(insites,4) - interstitial type: spiro 1-6
c     lower grafted 11-16
c     upper grafted 21-26
c
c     split:
c     ssites(snsites,1-3) - i,j,k indices of the site
c     ssites(insites,4) - up or down site (1 up, 0 down)

      integer sites(9999,5),isites(9999,4),ssites(9999,5),senpos(9999)
      integer nsites,insites,snsites,enpos(9999),ienpos(9999)

c     transition energies at each active vacancy site
c
c     energies(nsite,1) - barrier to move forward
c     energies(nsite,2) - barrier to move towards upper layer
c     energies(nsite,3) - barrier to move towards lower layer
c     energies(nsite,4) - barrier to move to upper layer
c     energies(nsite,5) - barrier to move to lower layer
c     un-crossable barrier denoted as 9 (eV)

c     total energy change associated with each transition

      real*8 energies(9999,5),menergies(9999,5)

c     transition energies at each active interstitial site
c
c     ienergies(nsite,1) - barrier to rotate clockwise
c     ienergies(nsite,2) - barrier to rotate counter-clockwise
c     ienergies(nsite,3) - barrier to move to nearest neighbour alpha site
c     if spiro:
c     ienergies(nsite,4) - barrier to move to up grafted
c     ienergies(nsite,5) - barrier to move to down grafted
c     ienergies(nsite,6) - barrier to close intimate frenkel pair
c     if grafted:
c     ienergies(nsite,4) - barrier to move to spiro
c     ienergies(nsite,5) - barrier to move to local (alpha) split
c     ienergies(nsite,6) - barrier to move to non-local (beta) split

c     total energy change associated with each transition

      real*8 ienergies(9999,6),imenergies(9999,6)

c     transition energies at each active split interstitial site
c
c     senergies(nsite,1) - barrier to move up | or /
c     senergies(nsite,2) - barrier to move up \ or |
c     senergies(nsite,3) - barrier to move up / or \
c     senergies(nsite,4) - barrier to move down | or /
c     senergies(nsite,5) - barrier to move down \ or |
c     senergies(nsite,6) - barrier to move down / or \

      real*8 senergies(9999,6),smenergies(9999,6)

c     the chosen pathway
c     move(1-3) - i,j,k indices of the site
c     move(4) - the movement direction:
c     
c     1-6 : lattice atom and direction to neighbouring site
c     7: lattice atom move towards upper layer
c     8: lattice atom move towards lower layer
c     9: lattice atom move to upper layer
c     10: lattice atom move to lower layer
c     11: split interstitial to move up | or /
c     12: split interstitial to move up \ or |
c     13: split interstitial to move up / or \
c     14: split interstitial to move down | or /
c     15: split interstitial to move down \ or |
c     16: split interstitial to move down / or \
c
c     21: interstitial to move clockwise
c     22: interstitial to move anti-clockwise
c     23: interstitial to move to neighbouring alpha
c     24: grafted to move to spiro
c         sprio to move up grafted
c     25: grafted to move to local alpha split
c         spiro to move down grafted
c     26: grafted to move to non-local beta split
c         spiro to close intimate pair

      integer move(4)

c     real simulation time - t
c     running total of system energy - en
      real*8 t,en,pfact,temp,tm_int,stopt,en_pos,avac_int,aint_int

      logical emin,do_sites,do_ens,ixyz,ipvacs,do_tm,do_stop
      logical add_vacs,add_ints,iprecomb,iplink,prnt_loss,prnt_recom
      logical recom_evnt,loss_evnt,loop_atoms,healed_recomb

c     tables of structures and energies
      integer*2 smat(21,30)
      integer*2 stab(21,30,12,499),istab(21,30,12,499)
      real*8 etab(12,499),ietab(6,499),mtab(12,499),imtab(6,499)

c     loop arrays and variables
      integer iloop(9999,4),print_threshold
      logical loopl
      real*8 ploop(9999)
      real*8 r_nucl
      real*8 r_init,r_capt,v_capt
      integer ix_rc(9999),iy_rc(9999),iz_rc(9999)
      integer ix_ls(9999),iy_ls(9999),iz_ls(9999)

c     seed for random number generator
      idum = -3

c-----initialising the arrays

      do i=1,800
         do j=1,800
            do k=0,4
               if(k.ge.1) then
                  system(i,j,k) = 0
               end if
               isystem(i,j,k) = 0
            end do
         end do
      end do
      
      do i=1,9999
         do j=1,5
            energies(i,j) = 0.0
            ienergies(i,j) = 0.0
            senergies(i,j) = 0.0
            sites(i,j) = 0
         end do
         enpos(i) = 0
         ienpos(i) = 0
         senpos(i) = 0
         senergies(i,6) = 0.0
         ienergies(i,6) = 0.0
         do j=1,4
            isites(i,j) = 0
            ssites(i,j) = 0
         end do
         do j=1,3
            iloop(i,j) = 0
         end do
         ploop(i) = 0.0
      end do

      t = 0.0
      en = 0.0
      ntab = 0
      intab = 0
      n_trj = 1
      n_ss = 1
      lst_vac = 0
      lst_int = 0
      itime = 1
      ivtime = 1
      iitime = 1
      icmb = 0
      nploop = 0
      recom_evnt = .false.
      loss_evnt = .false.

      write(*,*)'-------------------'
      write(*,*)'|    GraphDMC     |'
      write(*,*)'-------------------'
      write(*,*)
      write(*,*)'   T. Trevethan    '
      write(*,*)'     (c) 2013      '
      write(*,*)'                   '

c.....read input file

      call readinput(emin,nl,nsteps,temp,pfact,do_sites,do_ens,ntraj,
     &     nsys,iseed,ixd,iyd,izd,per_x,per_y,per_z,ixyz,ipvacs,
     &     do_tm,tm_int,do_stop,stopt,add_vacs,add_ints,avac_int,
     &     aint_int,num_avac,num_aint,iprecomb,iplink,loopl,r_nucl,
     &     r_init,r_capt,prnt_loss,loop_atoms,v_capt,healed_recomb,
     &     print_threshold)

      idum = iseed
      izd = nl

      if(per_z .and. mod(nl,2).ne.0) then
         write(*,*)'Cannot have odd number of layers and periodic in z'
         stop
      end if

      open(2,FILE='system.xyz')

      if(ipvacs) then
         open(59,FILE='vacs.xyz')
      end if

      if(iprecomb) then
         open(65,FILE='recomb.dat')
      end if 
      
      if(prnt_loss) then
         open(77,FILE='loss.dat')
      end if 
      
      if(iplink) then
         open(66,FILE='link.xyz')
      end if 

      if(do_sites) then
         open(5,FILE='sites.dat')
      end if
      if(do_ens) then
         open(8,FILE='energies.dat')
      end if
      open(81,FILE='tens.out')

      if(loopl) then
         open(76,FILE='loops.dat')
      end if

c.....reading in the structure database

      call readens(etab,stab,ntab,mtab,emin)
      call readiens(ietab,istab,intab,imtab,emin)

c.....initialising the system and printing structures

c      call makesystem(system,isystem,nat)
c
      write(*,*)' ... reading system file'

      call readsystem(system,isystem,nl,nat,ixd,iyd)

      write(*,'(a15,i10,a6)')'     there are ',nat,' atoms' 

c      call writesystem(system,isystem,nl)

      if(ipvacs) then
         call initpoints(system,isystem,nat,nl,ixd,iyd,n_vac,n_i,n_s)
         call printvacs(system,nl,ixd,iyd,n_vac)
      end if


      if(ixyz) then
         call printsystem(system,isystem,nat,nl,ixd,iyd,per_x,
     &        per_y,per_z)
      else
         call initpoints(system,isystem,nat,nl,ixd,iyd,n_vac,n_i,n_s)
         call printpoints(system,isystem,nat,nl,ixd,iyd,n_vac,n_i,n_s,
     &        loop_atoms,iloop,ploop,nploop,per_x,
     &        per_y,per_z,idum)
      end if
      
      if(loopl) then
         call printloops(iloop,ploop,nploop,nl,ixd,iyd,per_x,
     &        per_y,per_z)
      end if

      nrecom = 0
      nloss = 0

      if(iprecomb) then
         call printrecom(recom_evnt,t,ix_rc,iy_rc,iz_rc,nrecom)
      end if
      if(prnt_loss) then
         call printloss(loss_evnt,t,ix_ls,iy_ls,iz_ls,nloss)
      end if

c      stop

c      write(*,*)r_capt,r_init,r_nucl

c-----start the main KMC loop

      do nstep=1,nsteps

c.....calling getsites

c........at the first timestep, do a full scan of the system
         if(nstep.eq.1) then
            call getsites(system,sites,nsites,nl,ixd,iyd)
            call getisites(isystem,isites,insites,nl,ixd,iyd)
            call getssties(system,ssites,snsites,nl,ixd,iyd)
c........otherwise call local getsites
         else
            call getsites_lc(system,sites,nsites,nl,ixd,iyd,move)
            call getisites_lc(isystem,isites,insites,nl,ixd,iyd,move)
            call getssties_lc(system,ssites,snsites,nl,ixd,iyd,move)
         end if

         if(do_sites) then
            call write_sites(nsites,sites,insites,isites,snsites,
     &           ssites,nstep)
         end if

c--------doing loop coalescence

         if(loopl) then

            call loop_coalescence(iloop,ploop,nploop,r_nucl,r_init,ixd,
     &           iyd,per_x,per_y,per_z,system,isystem,isites,insites,
     &           ssites,snsites,r_capt,idum,en)
            call sortsites(isites,insites,ssites,snsites)

            call vacancy_loop(iloop,ploop,nploop,r_nucl,r_init,ixd,
     &           iyd,per_x,per_y,per_z,system,isystem,isites,insites,
     &           ssites,snsites,r_capt,idum,en,sites,nsites,v_capt)
            call sortnsites(sites,nsites)

         end if
         
         if((nsites+insites+snsites).eq.0) then
            if(loopl) then
               call printloops(iloop,ploop,nploop,nl,ixd,iyd,per_x,
     &              per_y,per_z)
            end if
            write(*,*)'done: no defects'
            stop
         end if
               
c-------labeling healed vacancy lines to direct recombination to ends

         if(healed_recomb) then
            call label_healed(system,sites,nsites)
         end if
         
c.....getting the energies

         call getens(system,isystem,sites,nsites,energies,etab,stab,
     &        ntab,nstep,menergies,mtab,emin,enpos)

         call getiens(system,isystem,isites,insites,ienergies,ietab,
     &        istab,intab,imenergies,imtab,emin,nstep,ienpos)

         call getsens(system,isystem,ssites,snsites,senergies,
     &     ietab,istab,intab,smenergies,imtab,emin,senpos)

          if(do_ens) then
            call write_energies(nsites,energies,insites,ienergies,
     &           snsites,senergies,nstep)
         end if

         if(healed_recomb) then
            call unlabel_healed(system,sites,nsites)
         end if

c........check for crossable barriers

         do ie=1,nsites
            do je=1,5
               if(energies(ie,je).lt.8.9) go to 195
            end do
         end do
         do ie=1,insites
            do je=1,6
               if(ienergies(ie,je).lt.8.9) go to 195
            end do
         end do
         do ie=1,snsites
            do je=1,6
               if(senergies(ie,je).lt.8.9) go to 195
            end do
         end do

         write(*,*)' - no accessible transitions: ground state -'
         stop

 195     continue
         
c.....KMC step

         call kmc(sites,nsites,energies,isites,insites,
     &     ienergies,ssites,snsites,senergies,move,idum,t,nstep,
     &        temp,pfact,en,menergies,imenergies,smenergies,
     &        enpos,ienpos,senpos,n_pos,en_pos)

         write(*,'(i8,5(a,i3),a,f4.2,a,f18.10,a,f12.2)')nstep,' x ',
     &        move(1),
     &        ' y ',move(2),' z ',move(3),' d ', move(4),' p ',n_pos,
     &        ' b ',en_pos,
     &        ' t ',t,' e ',en

c-----move atoms/vacancies

         if(move(4).le.10) then
            call maketrans(system,move)
            recom_evnt = .false.
         else if(move(4).ge.11 .and. move(4).le.16) then
         call strans(system,isystem,move)
            recom_evnt = .false.
         else if(move(4).ge.21 .and. move(4).le.26) then
            recom_evnt = .false.
            if(move(4).eq.26 .and. isystem(move(1),move(2),move(3)).ge.1 
     &           .and. isystem(move(1),move(2),move(3)).le.6) then
               icmb = icmb + 1
               recom_evnt = .true.
               nrecom = nrecom + 1
               ix_rc(nrecom) = move(1)
               iy_rc(nrecom) = move(2)
               iz_rc(nrecom) = move(3)
            end if
            call itrans(system,isystem,move)
         else
            write(*,*)'Error in move index'
            stop
         end if


c-----periodic boundary conditions
c.....if periodic then do nothing
c-----if not periodic in the x direction, remove defects
c     from the edges of the cell

         loss_evnt = .false.

         if(.not.per_x) then

            do i=1,2
               do j=1,iyd
                  do k=1,izd
                     if(system(i,j,k).ne.0) then
                        if(system(i,j,k).ne.1) then
                           system(i,j,k) = 1
                           lst_vac = lst_vac + 1
                           loss_evnt = .true.
                           nloss = nloss + 1
                           ix_ls(nloss) = 1
                           iy_ls(nloss) = j
                           iz_ls(nloss) = k
                           if(system(i,j,k).ne.9) then
                              do nnn=1,nsites
                                 if(sites(nnn,1).eq.i .and. 
     &                                sites(nnn,2).eq.j .and.
     &                                sites(nnn,3).eq.k) then
                                    n_pos = nnn
                                 end if
                              end do
                              sites(n_pos,1) = -999
                           else
                              do nnn=1,snsites
                                 if(ssites(nnn,1).eq.i .and. 
     &                                ssites(nnn,2).eq.j .and.
     &                                ssites(nnn,3).eq.k) then
                                    n_pos = nnn
                                 end if
                              end do
                              ssites(n_pos,1) = -999
                           end if
                        end if
                     end if
                     if(isystem(i,j,k).ne.0) then
                        if(isystem(i,j,k).ne.8) then
                           isystem(i,j,k) = 8
                           lst_int = lst_int + 1
                           loss_evnt = .true.
                           nloss = nloss + 1
                           ix_ls(nloss) = 1
                           iy_ls(nloss) = j
                           iz_ls(nloss) = k
                           do nnn=1,insites
                              if(isites(nnn,1).eq.i .and. 
     &                             isites(nnn,2).eq.j .and.
     &                             isites(nnn,3).eq.k) then
                                 n_pos = nnn
                              end if
                           end do
                           isites(n_pos,1) = -999
                        end if
                     end if
                  end do
               end do
            end do

            do i=ixd-1,ixd
               do j=1,iyd
                  do k=1,izd
                     if(system(i,j,k).ne.0) then
                        if(system(i,j,k).ne.1) then
                           system(i,j,k) = 1
                           lst_vac = lst_vac + 1
                           loss_evnt = .true.
                           nloss = nloss + 1
                           ix_ls(nloss) = ixd
                           iy_ls(nloss) = j
                           iz_ls(nloss) = k
                           if(system(i,j,k).ne.9) then
                              do nnn=1,nsites
                                 if(sites(nnn,1).eq.i .and. 
     &                                sites(nnn,2).eq.j .and.
     &                                sites(nnn,3).eq.k) then
                                    n_pos = nnn
                                 end if
                              end do
                              sites(n_pos,1) = -999
                           else
                              do nnn=1,snsites
                                 if(ssites(nnn,1).eq.i .and. 
     &                                ssites(nnn,2).eq.j .and.
     &                                ssites(nnn,3).eq.k) then
                                    n_pos = nnn
                                 end if
                              end do
                              ssites(n_pos,1) = -999
                           end if
                        end if
                     end if
                     if(isystem(i,j,k).ne.0) then
                        if(isystem(i,j,k).ne.8) then
                           isystem(i,j,k) = 8
                           lst_int = lst_int + 1
                           loss_evnt = .true.
                           nloss = nloss + 1
                           ix_ls(nloss) = ixd
                           iy_ls(nloss) = j
                           iz_ls(nloss) = k
                           do nnn=1,insites
                              if(isites(nnn,1).eq.i .and. 
     &                             isites(nnn,2).eq.j .and.
     &                             isites(nnn,3).eq.k) then
                                 n_pos = nnn
                              end if
                           end do
                           isites(n_pos,1) = -999
                        end if
                     end if
                  end do
               end do
            end do

         end if

c-----if not periodic in the y direction, remove defects
c     from the edges of the cell

         if(.not.per_y) then

            do j=1,2
               do i=1,ixd
                  do k=1,izd
                     if(system(i,j,k).ne.0) then
                       if(system(i,j,k).ne.1) then
                           system(i,j,k) = 1
                           lst_vac = lst_vac + 1
                           loss_evnt = .true.
                           nloss = nloss + 1
                           ix_ls(nloss) = i
                           iy_ls(nloss) = 1
                           iz_ls(nloss) = k
                           if(system(i,j,k).ne.9) then
                              do nnn=1,nsites
                                 if(sites(nnn,1).eq.i .and. 
     &                                sites(nnn,2).eq.j .and.
     &                                sites(nnn,3).eq.k) then
                                    n_pos = nnn
                                 end if
                              end do
                              sites(n_pos,1) = -999
                           else
                              do nnn=1,snsites
                                 if(ssites(nnn,1).eq.i .and. 
     &                                ssites(nnn,2).eq.j .and.
     &                                ssites(nnn,3).eq.k) then
                                    n_pos = nnn
                                 end if
                              end do
                              ssites(n_pos,1) = -999
                           end if
                        end if
                     end if
                     if(isystem(i,j,k).ne.0) then
                        if(isystem(i,j,k).ne.8) then
                           isystem(i,j,k) = 8
                           lst_int = lst_int + 1
                           loss_evnt = .true.
                           nloss = nloss + 1
                           ix_ls(nloss) = i
                           iy_ls(nloss) = 1
                           iz_ls(nloss) = k
                           do nnn=1,insites
                              if(isites(nnn,1).eq.i .and. 
     &                             isites(nnn,2).eq.j .and.
     &                             isites(nnn,3).eq.k) then
                                 n_pos = nnn
                              end if
                           end do
                           isites(n_pos,1) = -999
                        end if
                     end if
                  end do
               end do
            end do

            do j=iyd-1,iyd
               do i=1,ixd
                  do k=1,izd
                     if(system(i,j,k).ne.0) then
                        if(system(i,j,k).ne.1) then
                           system(i,j,k) = 1
                           lst_vac = lst_vac + 1
                           loss_evnt = .true.
                           nloss = nloss + 1
                           ix_ls(nloss) = i
                           iy_ls(nloss) = iyd
                           iz_ls(nloss) = k
                           if(system(i,j,k).ne.9) then
                              do nnn=1,nsites
                                 if(sites(nnn,1).eq.i .and. 
     &                                sites(nnn,2).eq.j .and.
     &                                sites(nnn,3).eq.k) then
                                    n_pos = nnn
                                 end if
                              end do
                              sites(n_pos,1) = -999
                           else
                              do nnn=1,snsites
                                 if(ssites(nnn,1).eq.i .and. 
     &                                ssites(nnn,2).eq.j .and.
     &                                ssites(nnn,3).eq.k) then
                                    n_pos = nnn
                                 end if
                              end do
                              ssites(n_pos,1) = -999
                           end if
                        end if
                     end if
                     if(isystem(i,j,k).ne.0) then
                        if(isystem(i,j,k).ne.8) then
                           isystem(i,j,k) = 8
                           lst_int = lst_int + 1
                           loss_evnt = .true.
                           nloss = nloss + 1
                           ix_ls(nloss) = i
                           iy_ls(nloss) = iyd
                           iz_ls(nloss) = k
                           do nnn=1,insites
                              if(isites(nnn,1).eq.i .and. 
     &                             isites(nnn,2).eq.j .and.
     &                             isites(nnn,3).eq.k) then
                                 n_pos = nnn
                              end if
                           end do
                           isites(n_pos,1) = -999
                        end if
                     end if
                  end do
               end do
            end do

         end if

c-----if not periodic in the z direction, remove interstitial defects
c     from the top layer

        if(.not.per_z) then            
           do i=1,ixd
              do j=1,iyd
                 if(isystem(i,j,izd).ne.0) then
              if(isystem(i,j,izd).ge.1 .and. isystem(i,j,izd).le.6) then
                       isystem(i,j,izd) = 8
                       lst_int = lst_int + 1
                       loss_evnt = .true.
                       nloss = nloss + 1
                       ix_ls(nloss) = i
                       iy_ls(nloss) = j
                       iz_ls(nloss) = izd
                       do nnn=1,insites
                          if(isites(nnn,1).eq.i .and. 
     &                         isites(nnn,2).eq.j .and.
     &                         isites(nnn,3).eq.k) then
                             n_pos = nnn
                          end if
                       end do
                       isites(n_pos,1) = -999
                    end if
                 end if
              end do
           end do
        end if

        call sortsites(isites,insites,ssites,snsites)
        call sortnsites(sites,nsites)

c.....printing structures

        if(do_tm) then
           
           if(t.gt.(itime*tm_int)) then
              if(ixyz) then
                 call printsystem(system,isystem,nat,nl,ixd,iyd,per_x,
     &                per_y,per_z)
                 if(loopl) then
                    call printloops(iloop,ploop,nploop,nl,ixd,iyd,per_x,
     &                   per_y,per_z)
                 end if
                 if(iprecomb) then
                 call printrecom(recom_evnt,t,ix_rc,iy_rc,iz_rc,nrecom)
                    nrecom = 0
                 end if
                 if(prnt_loss) then
                 call printloss(loss_evnt,t,ix_ls,iy_ls,iz_ls,nloss)
                    nloss = 0
                 end if
                 if(ipvacs) then
                    call printvacs(system,nl,ixd,iyd,n_vac)
                 end if
                 if(iplink) then
                    call printlinks(system,nl,ixd,iyd)
                 end if
              else
                 call printpoints(system,isystem,nat,nl,ixd,iyd,n_vac,
     &                n_i,n_s,
     &                loop_atoms,iloop,ploop,nploop,per_x,
     &                per_y,per_z,idum)
                 if(loopl) then
                    call printloops(iloop,ploop,nploop,nl,ixd,iyd,per_x,
     &                   per_y,per_z)
                 end if
                 if(iprecomb) then
                  call printrecom(recom_evnt,t,ix_rc,iy_rc,iz_rc,nrecom)
                    nrecom = 0
                 end if
                 if(prnt_loss) then
                    call printloss(loss_evnt,t,ix_ls,iy_ls,iz_ls,nloss)
                    nloss = 0
                 end if
                 if(iplink) then
                    call printlinks(system,nl,ixd,iyd)
                 end if
              end if
 979          itime = itime + 1
              if((t-((itime-1)*tm_int)).gt.tm_int) then
                 if(ixyz) then
                    call printsystem(system,isystem,nat,nl,ixd,iyd,per_x,
     &                   per_y,per_z)
                    if(loopl) then
                    call printloops(iloop,ploop,nploop,nl,ixd,iyd,per_x,
     &                      per_y,per_z)
                    end if
                    if(iprecomb) then
                  call printrecom(recom_evnt,t,ix_rc,iy_rc,iz_rc,nrecom)
                       nrecom = 0
                    end if
                    if(prnt_loss) then
                    call printloss(loss_evnt,t,ix_ls,iy_ls,iz_ls,nloss)
                       nloss = 0
                    end if
                    if(ipvacs) then
                       call printvacs(system,nl,ixd,iyd,n_vac)
                    end if
                    if(iplink) then
                       call printlinks(system,nl,ixd,iyd)
                    end if
                 else
                   call printpoints(system,isystem,nat,nl,ixd,iyd,n_vac,
     &                   n_i,n_s,
     &                   loop_atoms,iloop,ploop,nploop,per_x,
     &                   per_y,per_z,idum)
                    if(loopl) then
                    call printloops(iloop,ploop,nploop,nl,ixd,iyd,per_x,
     &                      per_y,per_z)
                    end if
                    if(iprecomb) then
                  call printrecom(recom_evnt,t,ix_rc,iy_rc,iz_rc,nrecom)
                       nrecom = 0
                    end if
                    if(prnt_loss) then
                     call printloss(loss_evnt,t,ix_ls,iy_ls,iz_ls,nloss)
                       nloss = 0
                    end if
                    if(iplink) then
                       call printlinks(system,nl,ixd,iyd)
                    end if
                 end if
                 go to 979
              else
                 go to 989
              end if
           end if

        else

           if(nstep.eq.(ntraj*n_trj)) then
              if(ixyz .and. nstep.ge.print_threshold) then
                 call printsystem(system,isystem,nat,nl,ixd,iyd,per_x,
     &                per_y,per_z)
                 if(loopl) then
                  call printloops(iloop,ploop,nploop,nl,ixd,iyd,per_x,
     &                   per_y,per_z)
                 end if
                 if(iprecomb) then
                 call printrecom(recom_evnt,t,ix_rc,iy_rc,iz_rc,nrecom)
                    nrecom = 0
                 end if
                 if(prnt_loss) then
                    call printloss(loss_evnt,t,ix_ls,iy_ls,iz_ls,nloss)
                    nloss = 0
                 end if
                 if(ipvacs) then
                    call printvacs(system,nl,ixd,iyd,n_vac)
                 end if
                 if(iplink) then
                    call printlinks(system,nl,ixd,iyd)
                 end if
              else if(nstep.ge.print_threshold) then
                 call printpoints(system,isystem,nat,nl,ixd,iyd,n_vac,
     &                n_i,n_s,
     &                loop_atoms,iloop,ploop,nploop,per_x,
     &                per_y,per_z,idum)
                 if(loopl) then
                   call printloops(iloop,ploop,nploop,nl,ixd,iyd,per_x,
     &                   per_y,per_z)
                 end if
                 if(iprecomb) then
                  call printrecom(recom_evnt,t,ix_rc,iy_rc,iz_rc,nrecom)
                    nrecom = 0
                 end if
                 if(prnt_loss) then
                    call printloss(loss_evnt,t,ix_ls,iy_ls,iz_ls,nloss)
                    nloss = 0
                 end if
                 if(iplink) then
                    call printlinks(system,nl,ixd,iyd)
                 end if
              end if
              n_trj=n_trj + 1
           end if
           
        end if

        if(nstep.eq.(nsys*n_ss)) then
        if(nstep.ge.print_threshold) then
           call writesystem(system,isystem,nl,ixd,iyd)
        end if
           n_ss=n_ss + 1
        end if         


c.....adding vacancies

        if(add_vacs) then
           if(t.gt.(ivtime*avac_int)) then
              call add_vacancies(system,nl,ixd,iyd,num_avac,idum,
     &             sites,nsites)
              ivtime = ivtime + 1
           end if
        end if


c.....adding interstitials

        if(add_ints) then
           if(t.gt.(iitime*aint_int)) then
              call add_interstitials(isystem,system,nl,ixd,iyd,
     &             num_aint,idum,isites,insites)
              iitime = iitime + 1
           end if
        end if

c--------writing time, energy and lost atoms to file

 989    write(81,'(f30.16,x,f12.2,3(x,i4))')t,en,lst_vac,lst_int,icmb

c--------stoping if stop time reached

        if(do_stop) then
        if(t.gt.stopt) then
           write(*,*)'time limit reached'
           stop
        end if
        end if

      end do

      end




      integer function md(i)
c-----------------------------------
c     y laterial periodic function -
c-----------------------------------
      integer ixd,iyd,izd
      logical per_x,per_y,per_z

      COMMON ixd,iyd,izd,per_x,per_y,per_z

      if(i.gt.iyd) then
         j = i - iyd
      else if(i.lt.1) then
         j = i + iyd
      else
         j = i
      end if

      md = j
      return
      end


      integer function mb(i)
c-----------------------------------
c     x laterial periodic function -
c-----------------------------------
      integer ixd,iyd,izd
      logical per_x,per_y,per_z

      COMMON ixd,iyd,izd,per_x,per_y,per_z

      if(i.gt.ixd) then
         j = i - ixd
      else if(i.lt.1) then
         j = i + ixd
      else
         j = i
      end if

      md = j
      return
      end


      integer function ml(i1)
c------------------------------------
c     c axis periodic function      -
c------------------------------------
      integer ixd,iyd,izd
      logical per_x,per_y,per_z

      COMMON ixd,iyd,izd,per_x,per_y,per_z

      if(i1.gt.izd) then
         i2 = i1 - izd
      else if(i1.lt.1) then
         i2 = i1 + izd
      else
         i2 = i1
      end if

      md = i2

      return
      end

      subroutine kmc(sites,nsites,energies,isites,insites,
     &     ienergies,ssites,snsites,senergies,move,idum,t,nstep,
     &     temp,pfact,en,menergies,imenergies,smenergies,
     &     enpos,ienpos,senpos,n_pos,en_pos)
c----------------------------------------------------------
c     subroutine that selects the transition from the     -
c     energies list using a Kinetic Monte Carlo Algoritm  -
c----------------------------------------------------------
      integer sites(9999,5),enpos(9999)
      integer isites(9999,4),ienpos(9999)
      integer ssites(9999,4),senpos(9999)
      real*8 energies(9999,5),menergies(9999,5)
      real*8 ienergies(9999,6),imenergies(9999,6)
      real*8 senergies(9999,6),smenergies(9999,6)
      real*8 en_list(9999),rate_list(9999)
      integer*2 nen(30)
      integer move(4)
      integer n_en,snsites
      real*8 prefactor,beta,Boltz,Temp
      real*8 rate_total,throw,run_tot
      real*8 dt,t,en,en_ch,pfact,en_pos

      Boltz = 1.380658d-23
      beta = 1.0/(temp*Boltz)
      n_en  = 0
      rate_total = 0.0
      n_en2 = 0




c.....create energy list

      do i=1,nsites

         do j=1,5

            n_en = n_en + 1

            en_list(n_en) = energies(i,j)*1.6e-19

         end do

      end do

      do i=1,snsites

         do j=1,6

            n_en = n_en + 1

            en_list(n_en) = senergies(i,j)*1.6e-19

         end do

      end do

      do i=1,insites

         do j=1,6

            n_en = n_en + 1

            en_list(n_en) = ienergies(i,j)*1.6e-19

         end do

      end do      

c-----calculate rates

      do i=1,n_en
         
         if(en_list(i).gt.1.0e-18) then
            rate_list(i) = 0.0
         else
            rate_list(i) = pfact*exp(-beta*en_list(i))
         end if

         rate_total = rate_total + rate_list(i)

      end do

c-----generate a random number between 0 and 1

      throw = ran1(idum)

c-----pick a transition

      run_tot = 0.0

      do i=1,n_en

         if(throw.gt.run_tot .and. throw.le.(run_tot + 
     &        rate_list(i)/rate_total)) then
            
            go to 12

         end if

         run_tot = run_tot + rate_list(i)/rate_total
         
      end do

      write(*,*)'Error in KMC'
      stop

c-----recover site and transition number

 12   do n=1,nsites

         do j=1,5
            
            n_en2 = n_en2 + 1
            
            if(n_en2.eq.i) then

               move(1) = sites(n,1)
               move(2) = sites(n,2)
               move(3) = sites(n,3)
               n_pos = enpos(n)
               if(j.eq.1) then
                  move(4) = sites(n,4)
               else if(j.eq.2) then
                  move(4) = 7
               else if(j.eq.3) then
                  move(4) = 8
               else if(j.eq.4) then
                  move(4) = 9
               else if(j.eq.5) then
                  move(4) = 10
               else
                  write(*,*)'error in recover imove'
                  stop
               end if
               
               en = en + menergies(n,j)
               en_pos = energies(n,j)

               go to 13

            end if
 
         end do
         
      end do

      do n=1,snsites

         do j=1,6
            
            n_en2 = n_en2 + 1
            
            if(n_en2.eq.i) then

               move(1) = ssites(n,1)
               move(2) = ssites(n,2)
               move(3) = ssites(n,3)
               n_pos = senpos(n)
               if(j.eq.1) then
                  move(4) = 11
               else if(j.eq.2) then
                  move(4) = 12
               else if(j.eq.3) then
                  move(4) = 13
               else if(j.eq.4) then
                  move(4) = 14
               else if(j.eq.5) then
                  move(4) = 15
               else if(j.eq.6) then
                  move(4) = 16
               else
                  write(*,*)'error in recover imove'
                  stop
               end if

               en = en + smenergies(n,j)
               en_pos = senergies(n,j)

               go to 13

            end if
 
         end do
         
      end do

      do n=1,insites

         do j=1,6
            
            n_en2 = n_en2 + 1
            
            if(n_en2.eq.i) then

               move(1) = isites(n,1)
               move(2) = isites(n,2)
               move(3) = isites(n,3)
               n_pos = ienpos(n)
               if(j.eq.1) then
                  move(4) = 21
               else if(j.eq.2) then
                  move(4) = 22
               else if(j.eq.3) then
                  move(4) = 23
               else if(j.eq.4) then
                  move(4) = 24
               else if(j.eq.5) then
                  move(4) = 25
               else if(j.eq.6) then
                  move(4) = 26
               else
                  write(*,*)'error in recover imove'
                  stop
               end if

               en = en + imenergies(n,j)
               en_pos = ienergies(n,j)

               go to 13

            end if
 
         end do
         
      end do

      continue

c.....calculate time

 13   throw = ran1(idum)

      dt = -log(throw)/rate_total
      
      t = t + dt

      end



      function ran1(idum)
c.....................................................................
c     minimal random number generator with Bays-Durham shuffle and    .
c     added safeguards. Returns a uniform deviate between 0.0 and    .
c     1.0 (exclusive of the endpoint values). Call with idum a       .
c     negative integer to initialise; thereafter, do not alter idum  .
c     between successive deviates in a sequence. RNMX should         .
c     approximate the largest floating value that is less than 1.    .
c.....................................................................
      integer idum,ia,im,iq,ir,ntab,ndiv
      real ran1,am,eps,rnmx
      parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836,
     &     ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-eps)
      integer j,k,iv(ntab),iy
      save iv,iy
      data iv /ntab*0/, iy /0/
      if(idum.le.0.or.iy.eq.0) then
         idum = max(-idum,1)
         do j=ntab + 8,1,-1
            k = idum/iq
            idum = ia*(idum-k*iq)-ir*k
            if(idum.lt.0) idum = idum + im
            if(j.le.ntab) iv(j) = idum
         end do
         iy = iv(1)
      end if
      k = idum/iq
      idum = ia*(idum - k*iq)-ir*k
      if(idum.lt.0) idum = idum + im
      j = 1 + iy/ndiv
      iy = iv(j)
      iv(j) = idum
      ran1 = min(am*iy,rnmx)
      return
      end



      subroutine add_vacancies(system,nl,ixd,iyd,num_avac,idum,
     &     sites,nsites)
c---------------------------------------------------------------
c     add num_avac vacancies at random positions in system     -
c---------------------------------------------------------------
      integer*1 system(800,800,199)
      integer sites(9999,5)

      iiv = 0

      do i=1,1000*num_avac
         
c........generating random positions

         dice1 = ran1(idum)
         dice2 = ran1(idum)
         dice3 = ran1(idum)

         iv = int(dice1*ixd)
         jv = int(dice2*iyd)
         kv = int(dice3*nl)

         kv = kv + 1

         if(iv.le.0) iv = 1
         if(jv.le.0) jv = 1
         if(kv.le.0) kv = 1
         if(iv.gt.ixd) iv = ixd
         if(jv.gt.iyd) jv = iyd
         if(kv.gt.nl) kv = nl

c         write(*,*)kv

c........if position has an occupied lattice site then
c...........make sure there are no second nearest neighbour
c           vacancies

         if(system(iv,jv,kv).eq.1) then
            
            n_crd = ivcoord(system,iv,jv,kv)
            
            if(n_crd.eq.3) then
               
               if(system(iv,md(jv+2),kv).ne.0) then

                  n_crd1 = ivcoord(system,iv,jv+2,kv)
                  n_crd2 = ivcoord(system,iv+1,jv-1,kv)
                  n_crd3 = ivcoord(system,iv-1,jv-1,kv)
                  
                  if((n_crd1+n_crd2+n_crd3).eq.9) then
                     iiv = iiv + 1
                     system(iv,jv,kv) = 2

                     sites(nsites+1,1) = iv
                     sites(nsites+1,2) = jv+2
                     sites(nsites+1,3) = kv
                     sites(nsites+1,4) = 4
                     sites(nsites+1,5) = 0

                     sites(nsites+2,1) = iv+1
                     sites(nsites+2,2) = jv-1
                     sites(nsites+2,3) = kv
                     sites(nsites+2,4) = 6
                     sites(nsites+2,5) = 0

                     sites(nsites+3,1) = iv-1
                     sites(nsites+3,2) = jv-1
                     sites(nsites+3,3) = kv
                     sites(nsites+3,4) = 2
                     sites(nsites+3,5) = 0

                     nsites = nsites + 3
                  end if

               else

                  n_crd1 = ivcoord(system,iv,jv-2,kv)
                  n_crd2 = ivcoord(system,iv+1,jv+1,kv)
                  n_crd3 = ivcoord(system,iv-1,jv+1,kv)

                  if((n_crd1+n_crd2+n_crd3).eq.9) then
                     iiv = iiv + 1
                     system(iv,jv,kv) = 2

                     sites(nsites+1,1) = iv
                     sites(nsites+1,2) = jv-2
                     sites(nsites+1,3) = kv
                     sites(nsites+1,4) = 1
                     sites(nsites+1,5) = 0

                     sites(nsites+2,1) = iv+1
                     sites(nsites+2,2) = jv+1
                     sites(nsites+2,3) = kv
                     sites(nsites+2,4) = 5
                     sites(nsites+2,5) = 0

                     sites(nsites+3,1) = iv-1
                     sites(nsites+3,2) = jv+1
                     sites(nsites+3,3) = kv
                     sites(nsites+3,4) = 3
                     sites(nsites+3,5) = 0

                     nsites = nsites + 3

                  end if
                  
               end if

            end if

         end if

         if(iiv.eq.num_avac) goto 10

      end do

      write(*,*)'could not place all vacancies'
      stop

 10   continue

      end


      subroutine add_interstitials(isystem,system,nl,ixd,iyd,
     &     num_aint,idum,isites,insites)
c---------------------------------------------------------------
c     add num_aint vacancies at random positions in isystem    -
c---------------------------------------------------------------
      integer*1 isystem(800,800,199)
      integer*1 system(800,800,199)
      integer isites(9999,4)

      iiv = 0

      do i=1,1000*num_aint
         
c........generating random positions

         dice1 = ran1(idum)
         dice2 = ran1(idum)
         dice3 = ran1(idum)

         iv = int(dice1*ixd)
         jv = int(dice2*iyd)
         kv = int(dice3*nl)

         if(iv.le.0) iv = 1
         if(jv.le.0) jv = 1
         if(kv.le.0) kv = 1
         if(iv.gt.ixd) iv = ixd
         if(jv.gt.iyd) jv = iyd
         if(kv.gt.nl) kv = nl

c........if position has an occupied lattice site then
c...........make sure there are no second nearest neighbour
c           vacancies

         if(isystem(iv,jv,kv).eq.8 .and. system(iv,jv,kv).eq.1 .and. 
     &        system(iv,jv,kv+1).eq.1 .and. 
     &        system(iv+1,jv+1,kv).le.1 .and. 
     &        system(iv-1,jv+1,kv).le.1 .and. 
     &        system(iv,jv+2,kv).le.1 .and. 
     &        system(iv+1,jv-1,kv).le.1 .and. 
     &        system(iv-1,jv-1,kv).le.1 .and. 
     &        system(iv,jv-2,kv).le.1 .and. 
     &        system(iv+1,jv+1,kv+1).le.1 .and. 
     &        system(iv-1,jv+1,kv+1).le.1 .and. 
     &        system(iv,jv+2,kv+1).le.1 .and. 
     &        system(iv+1,jv-1,kv+1).le.1 .and. 
     &        system(iv-1,jv-1,kv+1).le.1 .and. 
     &        system(iv,jv-2,kv+1).le.1 .and. 
     &        isystem(iv,jv,kv-1).eq.8 .and.
     &        isystem(iv,jv,kv+1).eq.8) then
            iiv = iiv + 1
            dice4 = ran1(idum)
            idir = int(dice4*5+1)
            if(idir.ge.6)idir = 6
            if(idir.le.1)idir = 1
            isystem(iv,jv,kv) = idir
            isites(insites+1,1) = iv
            isites(insites+1,2) = jv
            isites(insites+1,3) = kv
            isites(insites+1,4) = idir
            insites = insites + 1
         end if
         
         if(iiv.eq.num_aint) goto 11

      end do
      
      write(*,*)'could not place all interstitials'
      stop
   
 11   continue

      end


      function bond_l(x1,y1,z1,x2,y2,z2)
c...........................................................
c     function to return bond lengths given cartesian      .
c     coordinates and molecule atom index numbers of the   .
c     two atoms forming the bond                           .
c...........................................................
      real*8 x1,y1,z1,x2,y2,z2,xt,yt,zt
c.....calculate bond length
      xt=x1-x2
      yt=y1-y2
      zt=z1-z2
      bond_l=sqrt((xt*xt)+(yt*yt)+(zt*zt))
      end


      subroutine loop_coalescence(iloop,ploop,nploop,r_nucl,r_init,ixd,
     &           iyd,per_x,per_y,per_z,system,isystem,isites,insites,
     &           ssites,snsites,r_capt,idum,en)
c---------------------------------------------------------------------
c     subroutine to create loops for coalescing interstitials        -
c     the loops are stored in the ploop arrays (cartesian            -
c     coordinates and loop radius                                    -
c---------------------------------------------------------------------
      real*8 ploop(9999),r_nucl,r_init,x1,y1,x2,y2,zero,r_capt
      real*8 rloop,rloop2,en
      integer isites(9999,4),iloop(9999,4)
      integer ssites(9999,4),snsites,insites      
      integer*1 system(800,800,199),isystem(800,800,199)
      logical per_x,per_y,per_z

      zero = 0.0

c------------check if any interstitial is within the capture radius of a loop

c...........loop over all loops
      do nli=1,nploop

         ix1 = iloop(nli,1)
         iy1 = iloop(nli,2)
         iz1 = iloop(nli,3)
         rloop = ploop(nli)

c..............check interstitial plane
         do nlj=1,insites

            ix2 = isites(nlj,1)
            iy2 = isites(nlj,2)
            iz2 = isites(nlj,3)                  

            if(iz1.eq.iz2) then
               
               x1 = ix1*1.22
               y1 = iy1*0.71
               x2 = ix2*1.22
               y2 = iy2*0.71

               dist1 = bond_l(x1,y1,zero,x2,y2,zero)
               
               if(dist1.le.(rloop+r_capt)) then
c....................remove interstitial from the isystem array
                  isystem(ix2,iy2,iz2) = 8
                  isites(nlj,1) = -999
c......................increase the loop radius by one atom area
                  ploop(nli) = sqrt((3.142*ploop(nli)*ploop(nli) + 
     &                 2.58844)/3.142)
             en = en - 3.0
                  go to 101
               end if

c--------check periodic images
               if(per_x .and. per_y) then
                  do ip=-1,1
                     do jp=-1,1
                        if(ip.ne.0 .or. jp.ne.0) then

                           x1 = ix1*1.22
                           y1 = iy1*0.71
                           x2 = (ix2+ixd*ip)*1.22
                           y2 = (iy2+iyd*jp)*0.71
                           dist1 = bond_l(x1,y1,zero,x2,y2,zero)

                           if(dist1.le.(rloop+r_capt)) then
c....................remove interstitial from the isystem array
                              isystem(ix2,iy2,iz2) = 8
c......................increase the loop radius by one atom area
                  ploop(nli) = sqrt((3.142*ploop(nli)*ploop(nli) + 
     &                 2.58844)/3.142)
                              isites(nlj,1) = -999
                              en = en - 3.0
                              go to 101
                           end if
                        end if
                     end do
                  end do
               end if

            end if
 101        continue
         end do
            

c--------------check for splits
         do nlj=1,snsites
         
            ix2 = ssites(nlj,1)
            iy2 = ssites(nlj,2)
            iz2 = ssites(nlj,3)                  
            
            if((iz1.eq.iz2).or.(iz1.eq.(iz2+1))) then
               
               x1 = ix1*1.22
               y1 = iy1*0.71
               x2 = ix2*1.22
               y2 = iy2*0.71
               
               dist1 = bond_l(x1,y1,zero,x2,y2,zero)
               
               if(dist1.le.(rloop+r_capt)) then
c....................remove interstitial from the isystem array
                  system(ix2,iy2,iz2) = 1
c......................increase the loop radius by one atom area
                  ploop(nli) = sqrt((3.142*ploop(nli)*ploop(nli) + 
     &                 2.58844)/3.142)
                  ssites(nlj,1) = -999
                  en = en - 3.0
                  go to 102
               end if
      
c--------check periodic images
               if(per_x .and. per_y) then
                  do ip=-1,1
                     do jp=-1,1
                        if(ip.ne.0 .or. jp.ne.0) then

                           x1 = ix1*1.22
                           y1 = iy1*0.71
                           x2 = (ix2+ixd*ip)*1.22
                           y2 = (iy2+iyd*jp)*0.71
                           dist1 = bond_l(x1,y1,zero,x2,y2,zero)
               
                           if(dist1.le.(rloop+r_capt)) then
c....................remove interstitial from the isystem array
                              system(ix2,iy2,iz2) = 1
c......................increase the loop radius by one atom area
                  ploop(nli) = sqrt((3.142*ploop(nli)*ploop(nli) + 
     &                 2.58844)/3.142)
                              ssites(nlj,1) = -999
                              en = en - 3.0
                              go to 102
                           end if
                        end if
                     end do
                  end do
               end if
            end if
 102        continue
         end do
      end do


c-----------THEN check for dimerisation

c............find nucleation

c------------calculating separation between interstitials
      do nli=1,insites
         do nlj=nli+1,insites
c.................only consider if they are in the same layer
            if(isites(nli,3).eq.isites(nlj,3)) then
               ix1 = isites(nli,1)
               iy1 = isites(nli,2)
               iz1 = isites(nli,3)
               
               ix2 = isites(nlj,1)
               iy2 = isites(nlj,2)
               iz2 = isites(nlj,3)                 
               
               x1 = ix1*1.22
               y1 = iy1*0.71
               
               x2 = ix2*1.22
               y2 = iy2*0.71
                     
               dist1 = bond_l(x1,y1,zero,x2,y2,zero)
c--------------if two interstitials are separated by less than r_nucl
c     then nucleate a loop
               if(dist1.le.r_nucl) then

c                  write(*,*)'isite'
c                  write(*,*)isystem(ix1,iy1,iz1)
c                  write(*,*)isystem(ix2,iy2,iz2)
                  nploop = nploop + 1
c....................add an entry to the ploop array: coordinates and radius (r_init)
                  iloop(nploop,1) = ix1
                  iloop(nploop,2) = iy1
                  iloop(nploop,3) = iz1
                  ploop(nploop) = r_init
                  rn1 = ran1(idum)
                  iloop(nploop,4) = int(rn1*5.5+1.01)
                  if(iloop(nploop,4).gt.6)iloop(nploop,4) = 6
c....................remove both interstitials from the isystem array
                  isystem(ix1,iy1,iz1) = 8
                  isystem(ix2,iy2,iz2) = 8
                  isites(nli,1) = -999
                  isites(nlj,1) = -999
                  en = en - 3.0
                  goto 201
               end if
            
c--------check periodic images
               if(per_x .and. per_y) then
                  do ip=-1,1
                     do jp=-1,1
                        if(ip.ne.0 .or. jp.ne.0) then

                           x1 = ix1*1.22
                           y1 = iy1*0.71
                           x2 = (ix2+ixd*ip)*1.22
                           y2 = (iy2+iyd*jp)*0.71
                           dist1 = bond_l(x1,y1,zero,x2,y2,zero)
               
                           if(dist1.le.r_nucl) then
                              nploop = nploop + 1
c....................add an entry to the ploop array: coordinates and radius (r_init)
                              iloop(nploop,1) = ix1
                              iloop(nploop,2) = iy1
                              iloop(nploop,3) = iz1
                              ploop(nploop) = r_init
                              rn1 = ran1(idum)
                              iloop(nploop,4) = int(rn1*5.5+1.01)
                            if(iloop(nploop,4).gt.6)iloop(nploop,4) = 6
c....................remove both interstitials from the isystem array
                              isystem(ix1,iy1,iz1) = 8
                              isystem(ix2,iy2,iz2) = 8
                              isites(nli,1) = -999
                              isites(nlj,1) = -999
                              en = en - 3.0
                              goto 201
                          end if
                        end if
                     end do
                  end do
               end if

            end if
 201        continue
         end do
      end do

c--------------check for splits
      do nli=1,insites
         do nlj=1,snsites
            ix1 = isites(nli,1)
            iy1 = isites(nli,2)
            iz1 = isites(nli,3)
         
            ix2 = ssites(nlj,1)
            iy2 = ssites(nlj,2)
            iz2 = ssites(nlj,3)                  
            if((iz1.eq.iz2).or.(iz1.eq.(iz2-1))) then
               x1 = ix1*0.71
               y1 = iy1*1.22
               x2 = ix2*0.71
               y2 = iy2*1.22
            
               dist1 = bond_l(x1,y1,zero,x2,y2,zero)
               
c               if(ix1.eq.14 .and. iy1.eq.22 .and. iz1.eq.4) then
c               write(*,*)ix2,iy2,iz2,dist1
c               end if

c--------------if two interstitials are separated by less than r_nucl
c              then nucleate a loop
               if(dist1.le.r_nucl) then
                  nploop = nploop + 1
c....................add an entry to the ploop array: coordinates and radius (r_init)
                  iloop(nploop,1) = ix1
                  iloop(nploop,2) = iy1
                  iloop(nploop,3) = iz1
                  ploop(nploop) = r_init
                  rn1 = ran1(idum)
                  iloop(nploop,4) = int(rn1*5.5+1.01)
                  if(iloop(nploop,4).gt.6)iloop(nploop,4) = 6
c....................remove both interstitials from the isystem array
                  isystem(ix1,iy1,iz1) = 8
                  system(ix2,iy2,iz2) = 1
                  isites(nli,1) = -999
                  ssites(nlj,1) = -999
                  en = en - 3.0
                  goto 202
               end if


c--------check periodic images
               if(per_x .and. per_y) then
                  do ip=-1,1
                     do jp=-1,1
                        if(ip.ne.0 .or. jp.ne.0) then

                           x1 = ix1*1.22
                           y1 = iy1*0.71
                           x2 = (ix2+ixd*ip)*1.22
                           y2 = (iy2+iyd*jp)*0.71
                           dist1 = bond_l(x1,y1,zero,x2,y2,zero)
               
                           if(dist1.le.r_nucl) then
                              nploop = nploop + 1
c....................add an entry to the ploop array: coordinates and radius (r_init)
                              iloop(nploop,1) = ix1
                              iloop(nploop,2) = iy1
                              iloop(nploop,3) = iz1
                              ploop(nploop) = r_init
                              rn1 = ran1(idum)
                              iloop(nploop,4) = int(rn1*5.5+1.01)
                             if(iloop(nploop,4).gt.6)iloop(nploop,4) = 6
c....................remove both interstitials from the isystem array
                              isystem(ix1,iy1,iz1) = 8
                              system(ix2,iy2,iz2) = 1
                              isites(nli,1) = -999
                              ssites(nlj,1) = -999
                              en = en - 3.0
                              goto 202
                           end if
                        end if
                     end do
                  end do
               end if

            end if
 202        continue
         end do
      end do
      
c------------calculating separation between split interstitials
      do nli=1,snsites
         do nlj=nli+1,snsites
c.................only consider if they are in the same layer
            if(isites(nli,3).eq.isites(nlj,3)) then
               ix1 = ssites(nli,1)
               iy1 = ssites(nli,2)
               iz1 = ssites(nli,3)
               
               ix2 = ssites(nlj,1)
               iy2 = ssites(nlj,2)
               iz2 = ssites(nlj,3)                 
                     
               x1 = ix1*1.22
               y1 = iy1*0.71
               x2 = ix2*1.22
               y2 = iy2*0.71
               dist1 = bond_l(x1,y1,zero,x2,y2,zero)
c--------------if two interstitials are separated by less than r_nucl
c              then nucleate a loop
               if(dist1.le.r_nucl) then
c                  write(*,*)'split_isite'
                  nploop = nploop + 1
c....................add an entry to the ploop array: coordinates and radius (r_init)
                  iloop(nploop,1) = ix1
                  iloop(nploop,2) = iy1
                  iloop(nploop,3) = iz1
                  ploop(nploop) = r_init
                  rn1 = ran1(idum)
                  iloop(nploop,4) = int(rn1*5.5+1.01)
                  if(iloop(nploop,4).gt.6)iloop(nploop,4) = 6
c....................remove both interstitials from the isystem array
                  system(ix1,iy1,iz1) = 1
                  system(ix2,iy2,iz2) = 1
                  ssites(nli,1) = -999
                  ssites(nlj,1) = -999
                  en = en - 3.0
                  goto 203
               end if

c--------check periodic images
               if(per_x .and. per_y) then
                  do ip=-1,1
                     do jp=-1,1
                        if(ip.ne.0 .or. jp.ne.0) then

                           x1 = ix1*1.22
                           y1 = iy1*0.71
                           x2 = (ix2+ixd*ip)*1.22
                           y2 = (iy2+iyd*jp)*0.71
                           dist1 = bond_l(x1,y1,zero,x2,y2,zero)
               
                           if(dist1.le.r_nucl) then
                              nploop = nploop + 1
c....................add an entry to the ploop array: coordinates and radius (r_init)
                              iloop(nploop,1) = ix1
                              iloop(nploop,2) = iy1
                              iloop(nploop,3) = iz1
                              ploop(nploop) = r_init
                              rn1 = ran1(idum)
                              iloop(nploop,4) = int(rn1*5.5+1.01)
                            if(iloop(nploop,4).gt.6)iloop(nploop,4) = 6
c....................remove both interstitials from the isystem array
                              system(ix1,iy1,iz1) = 1
                              system(ix2,iy2,iz2) = 1
                              ssites(nli,1) = -999
                              ssites(nlj,1) = -999
                              en = en - 3.0
                              goto 203
                           end if
                        end if
                     end do
                  end do
               end if

            end if
 203        continue
         end do
      end do

c------------------What to do when two loops meet/overlap
c     check for touching loops and quit if they are

      do nli=1,nploop

         ix1 = iloop(nli,1)
         iy1 = iloop(nli,2)
         iz1 = iloop(nli,3)
         rloop = ploop(nli)

         do nlj=nli,nploop
            
            if(nli.ne.nlj) then

               ix2 = iloop(nlj,1)
               iy2 = iloop(nlj,2)
               iz2 = iloop(nlj,3)
               rloop2 = ploop(nlj)

               if(iz1.eq.iz2) then

                  x1 = ix1*1.22
                  y1 = iy1*0.71
                  x2 = ix2*1.22
                  y2 = iy2*0.71
               
                  dist1 = bond_l(x1,y1,zero,x2,y2,zero)               
                  
                  if(dist1.lt.(rloop + rloop2)) then

                    call printloops(iloop,ploop,nploop,nl,ixd,iyd,per_x,
     &                    per_y,per_z)

                    write(*,*)'touching loops: ',nli,nlj,x1,y1,x2,y2,
     &                   rloop,rloop2,ix1,iy1,ix2,iy2
                    write(*,*)'here1',nploop
                    stop
                  end if

c--------check periodic images
                  if(per_x .and. per_y) then
                     do ip=-1,1
                        do jp=-1,1
                           if(ip.ne.0 .or. jp.ne.0) then

                              x2 = (ix2+ixd*ip)*1.22
                              y2 = (iy2+iyd*jp)*0.71
                              dist1 = bond_l(x1,y1,zero,x2,y2,zero)

                              if(dist1.lt.(rloop + rloop2)) then

                    call printloops(iloop,ploop,nploop,nl,ixd,iyd,per_x,
     &                                per_y,per_z)

                    write(*,*)'touching loops: ',nli,nlj
                    stop
                              end if
                           end if
                        end do
                     end do
                  end if
               end if
            end if
         end do
      end do

c--------if no PBC then quit if any loop touches the sides of the cell

      if(.not.per_x .or. .not.per_y) then

         do nli=1,nploop

            ix1 = iloop(nli,1)
            iy1 = iloop(nli,2)
            rloop = ploop(nli)         

            x1 = ix1*1.22
            y1 = iy1*0.71
            x2 = ixd*1.22
            y2 = iyd*0.71

            if((ix1 - rloop).le.0 .or. (iy1 - rloop).le.0 .or.
     &           (ix1 + rloop).gt.x2 .or. (iy1 + rloop).gt.y2) then
               call printloops(iloop,ploop,nploop,nl,ixd,iyd,per_x,
     &              per_y,per_z)
               write(*,*)'loop touching boundaries: ',nli
               stop
            end if      
         end do
      end if

      end 

      subroutine sortsites(isites,insites,ssites,snsites)
c--------------------------------------------------------------------
c     subroutine to sort the sites arrays to remove lost defects    -
c--------------------------------------------------------------------
      integer isites(9999,4)
      integer ssites(9999,5),snsites,insites
      integer tmpsites(9999,5)

      itsites = 0

      do i=1,insites
         if(isites(i,1).gt.-900) then
            itsites = itsites + 1
            tmpsites(itsites,1) = isites(i,1)
            tmpsites(itsites,2) = isites(i,2)
            tmpsites(itsites,3) = isites(i,3)
            tmpsites(itsites,4) = isites(i,4)
         end if
      end do

      do i=1,itsites
         isites(i,1) = tmpsites(i,1)
         isites(i,2) = tmpsites(i,2)
         isites(i,3) = tmpsites(i,3)
         isites(i,4) = tmpsites(i,4)
      end do
      
      insites = itsites

      itsites = 0

      do i=1,snsites
         if(ssites(i,1).gt.-900) then
            itsites = itsites + 1
            tmpsites(itsites,1) = ssites(i,1)
            tmpsites(itsites,2) = ssites(i,2)
            tmpsites(itsites,3) = ssites(i,3)
            tmpsites(itsites,4) = ssites(i,4)
         end if
      end do



      do i=1,itsites
         ssites(i,1) = tmpsites(i,1)
         ssites(i,2) = tmpsites(i,2)
         ssites(i,3) = tmpsites(i,3)
         ssites(i,4) = tmpsites(i,4)
      end do

      snsites = itsites

      end 




      subroutine sortnsites(sites,nsites)
c--------------------------------------------------------------------
c     subroutine to sort the sites arrays to remove lost defects    -
c--------------------------------------------------------------------
      integer sites(9999,5),nsites
      integer tmpsites(9999,5)

      itsites = 0

      do i=1,nsites
         if(sites(i,1).gt.-900) then
            itsites = itsites + 1
            tmpsites(itsites,1) = sites(i,1)
            tmpsites(itsites,2) = sites(i,2)
            tmpsites(itsites,3) = sites(i,3)
            tmpsites(itsites,4) = sites(i,4)
            tmpsites(itsites,5) = sites(i,5)
         end if
      end do

      do i=1,itsites
         sites(i,1) = tmpsites(i,1)
         sites(i,2) = tmpsites(i,2)
         sites(i,3) = tmpsites(i,3)
         sites(i,4) = tmpsites(i,4)
         sites(i,5) = tmpsites(i,5)
      end do
      
      nsites = itsites

      end 

      subroutine label_healed(system,sites,nsites)
c---------------------------------------------------------------------
c     subroutine to label vacancy lines with healed status, which    -
c     enables interstitials to diffuse across them and only recombine-
c     with the end vacancies                                         -
c---------------------------------------------------------------------
      integer sites(9999,5),nsites      
      integer*1 system(800,800,199)
      integer vsites(9999,3),vcoord

c-----create vacancy list

      do i=1,nsites
         if(sites(i,4).eq.1) then
            ixadd = 0
            iyadd = 2
         else if(sites(i,4).eq.2) then
            ixadd = 1
            iyadd = 1
         else if(sites(i,4).eq.3) then
            ixadd = 1
            iyadd = -1
         else if(sites(i,4).eq.4) then
            ixadd = 0
            iyadd = -2
         else if(sites(i,4).eq.5) then
            ixadd = -1
            iyadd = -1
         else if(sites(i,4).eq.6) then
            ixadd = -1
            iyadd = 1
         end if

         vsites(i,1) = sites(i,1) + ixadd
         vsites(i,2) = sites(i,2) + iyadd
         vsites(i,3) = sites(i,3)

      end do

      do i=1,nsites
         vcoord = 0
         do j=1,nsites
            if(i.ne.j) then
               idx = vsites(i,1) - vsites(j,1)
               idy = vsites(i,2) - vsites(j,2)
               if(vsites(i,3).eq.vsites(j,3) .and. abs(idx).le.1 .and.
     &              abs(idy).le.2) then
                  vcoord = vcoord + 1
               end if
            end if
         end do
         
         if(vcoord.eq.2) then
            if(system(vsites(i,1),vsites(i,2),vsites(i,3)).eq.2 .or.
     &           system(vsites(i,1),vsites(i,2),vsites(i,3)).eq.7) then
               system(vsites(i,1),vsites(i,2),vsites(i,3)) = 7
            else
               write(*,*)'error in label healed'
               stop
            end if
         else if(vcoord.eq.1) then
            if(system(vsites(i,1),vsites(i,2),vsites(i,3)).eq.2 .or.
     &           system(vsites(i,1),vsites(i,2),vsites(i,3)).eq.7) then
               system(vsites(i,1),vsites(i,2),vsites(i,3)) = 2
            else
               write(*,*)'error in label healed 2'
               stop
            end if
         end if
      end do

      end

      subroutine unlabel_healed(system,sites,nsites)
c---------------------------------------------------------------------               
c     subroutine to label vacancy lines with healed status, which    -   
c     enables interstitials to diffuse across them and only recombine-        
c     with the end vacancies                                         - 
c---------------------------------------------------------------------                         
      integer sites(9999,5),nsites
      integer*1 system(800,800,199)
      integer vsites(9999,3),vcoord

c-----create vacancy list                       

      do i=1,nsites
         if(sites(i,4).eq.1) then
            ixadd = 0
            iyadd = 2
         else if(sites(i,4).eq.2) then
            ixadd = 1
            iyadd = 1
         else if(sites(i,4).eq.3) then
            ixadd = 1
            iyadd = -1
         else if(sites(i,4).eq.4) then
            ixadd = 0
            iyadd = -2
         else if(sites(i,4).eq.5) then
                ixadd = -1
            iyadd = -1
         else if(sites(i,4).eq.6) then
            ixadd = -1
            iyadd = 1
         end if

         vsites(i,1) = sites(i,1) + ixadd
         vsites(i,2) = sites(i,2) + iyadd
         vsites(i,3) = sites(i,3)

      end do

      do i=1,nsites
         if(system(vsites(i,1),vsites(i,2),vsites(i,3)).eq.7) then
            system(vsites(i,1),vsites(i,2),vsites(i,3)) = 2
         end if
      end do
      
      end


      subroutine vacancy_loop(iloop,ploop,nploop,r_nucl,r_init,ixd,
     &           iyd,per_x,per_y,per_z,system,isystem,isites,insites,
     &           ssites,snsites,r_capt,idum,en,sites,nsites,v_capt)
c---------------------------------------------------------------------
c     subroutine to recombine single vacancies with loop edges       -
c---------------------------------------------------------------------
      real*8 ploop(9999),r_nucl,r_init,x1,y1,x2,y2,zero,r_capt
      real*8 rloop,rloop2,en,v_capt
      integer isites(9999,4),iloop(9999,4),nsites
      integer ssites(9999,4),snsites,insites,sites(9999,5)
      integer*1 system(800,800,199),isystem(800,800,199)
      logical per_x,per_y,per_z
      integer itmplp(9999,3),vacloc(9999,3),nvacloc
      real*8 ptmplp(9999)

      zero = 0.0
      nvacloc = 0

c------------check if any SINGLE vacancy is in the capture radius of a loop

c...........loop over all loops
      do nli=1,nploop

         ix1 = iloop(nli,1)
         iy1 = iloop(nli,2)
         iz1 = iloop(nli,3)
         rloop = ploop(nli)

         nvacloc = 0

c..............check lower vacancy plane
         do nlj=1,nsites

            ix2 = sites(nlj,1)
            iy2 = sites(nlj,2)
            iz2 = sites(nlj,3)

            if((iz1.eq.iz2).or.((iz1+1).eq.iz2)) then

               if(sites(nlj,4).eq.1) then
                  ixadd = 0
                  iyadd = 2
               else if(sites(nlj,4).eq.2) then
                  ixadd = 1
                  iyadd = 1
               else if(sites(nlj,4).eq.3) then
                  ixadd = 1
                  iyadd = -1
               else if(sites(nlj,4).eq.4) then
                  ixadd = 0
                  iyadd = -2
               else if(sites(nlj,4).eq.5) then
                  ixadd = -1
                  iyadd = -1
               else if(sites(nlj,4).eq.6) then
                  ixadd = -1
                  iyadd = 1
              end if

              iv2 = ivcoord2(system,ix2+ixadd,iy2+iyadd,iz2)

              if(iv2.eq.3) then
               
c                 do jjj=1,nvacloc
c                    if(vacloc(jjj,1).eq.(ix2+ixadd) .and. 
c     &                   vacloc(jjj,2).eq.(iy2+iyadd) .and. 
c     &                   vacloc(jjj,3).eq.(iz2)) then
c                       sites(nlj,1) = -999
c                       go to 101
c                    end if
c                 end do

c                 nvacloc = nvacloc + 1
c                 vacloc(nvacloc,1) = ix2+ixadd
c                 vacloc(nvacloc,2) = iy2+iyadd
c                 vacloc(nvacloc,3) = iz2
                 
                 x1 = ix1*1.22
                 y1 = iy1*0.71
                 x2 = (ix2+ixadd)*1.22
                 y2 = (iy2+iyadd)*0.71
                 
                 dist1 = bond_l(x1,y1,zero,x2,y2,zero)

c                  write(*,*)x1,y1,x2,y2,dist1

                 if(dist1.le.(rloop+v_capt) .and. dist1.gt.rloop) then
c....................remove vacancy from the system array
                    system(mb(ix2+ixadd),md(iy2+iyadd),iz2) = 1
                    sites(nlj,1) = -999
c......................decrease the loop radius by one atom area
                    ploop(nli) = sqrt((3.142*ploop(nli)*ploop(nli) - 
     &                   2.58844/3.0)/3.142)
                    write(*,*)nlj,ploop(nli)
                    en = en - 11.0/3.0
                    go to 101
                 end if

c--------check periodic images

                  if(per_x .and. per_y) then
                     do ip=-1,1
                        do jp=-1,1
                           if(ip.ne.0 .or. jp.ne.0) then

                              x1 = ix1*1.22
                              y1 = iy1*0.71
                              x2 = (ix2+ixd*ip+ixadd)*1.22
                              y2 = (iy2+iyd*jp+iyadd)*0.71
                              dist1 = bond_l(x1,y1,zero,x2,y2,zero)
                              
                              if(dist1.le.(rloop+v_capt) .and. 
     &                             dist1.gt.rloop) then
c....................remove vacancy from the isystem array
                            system(mb(ix2+ixadd),md(iy2+iyadd),iz2) = 1
c......................decrease the loop radius by one atom area
                        ploop(nli) = sqrt((3.142*ploop(nli)*ploop(nli) - 
     &                           2.58844/3.0)/3.142)
                                 sites(nlj,1) = -999
                                 en = en - 11.0/3.0
                                 go to 101
                              end if
                           end if
                        end do
                     end do
                  end if

               end if
            end if
 101        continue
         end do

c--------release interstitial and remove loop if radius drops below 1

         if(ploop(nli).lt.1.0) then

            if(isystem(ix1,iy1,iz1).eq.8) then
               isystem(ix1,iy1,iz1) = 1
               insites = insites + 1
               isites(insites,1) = ix1
               isites(insites,2) = iy1
               isites(insites,3) = iz1
               isites(insites,4) = 1
            else if(isystem(ix1,md(iy1+2),iz1).eq.8) then
               isystem(ix1,md(iy1+2),iz1) = 1
               insites = insites + 1
               isites(insites,1) = ix1
               isites(insites,2) = md(iy1+2)
               isites(insites,3) = iz1
               isites(insites,4) = 1
            else if(isystem(ix1,md(iy1+2),iz1).eq.8) then
               isystem(ix1,md(iy1+2),iz1) = 1
               insites = insites + 1
               isites(insites,1) = ix1
               isites(insites,2) = md(iy1-2)
               isites(insites,3) = iz1
               isites(insites,4) = 1
            else
               write(*,*)'error in release int'
               stop
            end if
         end if
      
      end do

      itl = 0

      do i=1,nploop
         if(ploop(i).gt.1.0) then
            itl = itl + 1
            ptmplp(itl) = ploop(i)
            itmplp(itl,1) = iloop(i,1) 
            itmplp(itl,2) = iloop(i,2)
            itmplp(itl,3) = iloop(i,3)
         end if
         ploop(i) = 0.0
         iloop(i,1) = 0
         iloop(i,2) = 0
         iloop(i,3) = 0
      end do

      do i=1,itl
         ploop(i) = ptmplp(i)
         iloop(i,1) = itmplp(i,1)
         iloop(i,2) = itmplp(i,2)
         iloop(i,3) = itmplp(i,3)
      end do

      nploop = itl

      end 
