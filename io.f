      subroutine readinput(emin,nl,nsteps,temp,pfact,do_sites,do_ens,
     &     ntraj,nsys,iseed,ixd,iyd,izd,per_x,per_y,per_z,ixyz,ipvacs,
     &     do_tm,tm_int,do_stop,stopt,add_vacs,add_ints,avac_int,
     &     aint_int,num_avac,num_aint,iprecomb,iplink,loopl,r_nucl,
     &     r_init,r_capt,prnt_loss,loop_atoms,v_capt,healed_recomb,
     &     print_threshold)
c----------------------------------------------------
c     read input file
c----------------------------------------------------
      real*8 temp,pfact
      logical emin,do_sites,do_ens,per_x,per_y,per_z,ixyz,ipvacs,do_tm
      logical do_stop,add_vacs,add_ints,iprecomb,iplink,loopl
      logical prnt_loss,prnt_recom,loop_atoms,healed_recomb
      character Line*120
      real*8 tm_int,stopt,avac_int,aint_int,r_nucl,r_init,r_capt,v_capt
      integer LinPos(30),LinEnd(30),print_threshold

c-----default values

      emin = .false.
      nl = 4
      nsteps = 1000
      temp = 273.0
      pfact = 1.0e12
      do_sites = .false.
      do_ens = .false.
      ntraj = 1
      nsys = 100
      iseed = -6
      ixd = 60
      iyd = 60
      izd = 4
      per_x = .true.
      per_y = .true.
      per_z = .true.
      ixyz = .true.
      ipvacs = .false.
      iplink = .false.
      iprecomb = .false.
      do_tm = .false.
      do_stop = .false.
      tm_int = 0.001
      add_vacs = .false.
      add_ints = .false.
      avac_int = 0.01
      aint_int = 0.01
      num_avac = 10
      num_aint = 10
      loopl = .false.
      r_nucl = 5.0
      r_init = 1.0
      r_capt = 6.0
      v_capt = 3.0
      prnt_loss = .false.
      loop_atoms = .true.
      healed_recomb = .false.
      print_threshold = 0

      open(82,FILE='gkin.dat')

      iw=0

      do k=1,40
         call CutStr(Line,NumLin,LinPos,LinEnd,82,iErr)
         if(iErr.ne.0) go to 50
         if(index(line,'end').ne.0) go to 100
         if(index(line(LinPos(1):LinEnd(1)),'en_min').ne.0) then
            if(NumLin.lt.2) go to 50
            if(index(line(LinPos(2):LinEnd(2)),'on').ne.0) then
               emin = .true.
            else if(index(line(LinPos(2):LinEnd(2)),'off').ne.0) then
               emin = .false.
            end if
            iw=iw+1
         else if(index(line(LinPos(1):LinEnd(1)),'num_layers').ne.0)then
            if(NumLin.lt.2) go to 50
            read(Line(LinPos(2):),*,END=50,ERR=50) nl
            iw=iw+1
         else if(index(line(LinPos(1):LinEnd(1)),'num_x').ne.0)then
            if(NumLin.lt.2) go to 50
            read(Line(LinPos(2):),*,END=50,ERR=50) ixd
            iw=iw+1
         else if(index(line(LinPos(1):LinEnd(1)),'num_y').ne.0)then
            if(NumLin.lt.2) go to 50
            read(Line(LinPos(2):),*,END=50,ERR=50) iyd
            iw=iw+1
         else if(index(line(LinPos(1):LinEnd(1)),'num_steps').ne.0) then
            if(NumLin.lt.2) go to 50
            read(Line(LinPos(2):),*,END=50,ERR=50) nsteps
            iw=iw+1
         else if(index(line(LinPos(1):LinEnd(1)),'temp').ne.0) then
            if(NumLin.lt.2) go to 50
            read(Line(LinPos(2):),*,END=50,ERR=50) temp
            iw=iw+1
         else if(index(line(LinPos(1):LinEnd(1)),'pre_factor').ne.0)then
            if(NumLin.lt.2) go to 50
            read(Line(LinPos(2):),*,END=50,ERR=50) pfact
            iw=iw+1
         else if(index(line(LinPos(1):LinEnd(1)),'sites').ne.0) then
            if(NumLin.lt.2) go to 50
            if(index(line(LinPos(2):LinEnd(2)),'on').ne.0) then
               do_sites = .true.
            else if(index(line(LinPos(2):LinEnd(2)),'off').ne.0) then
               do_sites = .false.
            end if
            iw=iw+1
         else if(index(line(LinPos(1):LinEnd(1)),'energies').ne.0) then
            if(NumLin.lt.2) go to 50
            if(index(line(LinPos(2):LinEnd(2)),'on').ne.0) then
               do_ens = .true.
            else if(index(line(LinPos(2):LinEnd(2)),'off').ne.0) then
               do_ens = .false.
            end if
            iw=iw+1
         else if(index(line(LinPos(1):LinEnd(1)),'periodic_x').ne.0)then
            if(NumLin.lt.2) go to 50
            if(index(line(LinPos(2):LinEnd(2)),'on').ne.0) then
               per_x = .true.
            else if(index(line(LinPos(2):LinEnd(2)),'off').ne.0) then
               per_x = .false.
            end if
            iw=iw+1
         else if(index(line(LinPos(1):LinEnd(1)),'periodic_y').ne.0)then
            if(NumLin.lt.2) go to 50
            if(index(line(LinPos(2):LinEnd(2)),'on').ne.0) then
               per_y = .true.
            else if(index(line(LinPos(2):LinEnd(2)),'off').ne.0) then
               per_y = .false.
            end if
            iw=iw+1
         else if(index(line(LinPos(1):LinEnd(1)),'periodic_z').ne.0)then
            if(NumLin.lt.2) go to 50
            if(index(line(LinPos(2):LinEnd(2)),'on').ne.0) then
               per_z = .true.
            else if(index(line(LinPos(2):LinEnd(2)),'off').ne.0) then
               per_z = .false.
            end if
            iw=iw+1
         else if(index(line(LinPos(1):LinEnd(1)),'xyz_full').ne.0) then
            if(NumLin.lt.2) go to 50
            if(index(line(LinPos(2):LinEnd(2)),'on').ne.0) then
               ixyz = .true.
            else if(index(line(LinPos(2):LinEnd(2)),'off').ne.0) then
               ixyz = .false.
            end if
            iw=iw+1
         else if(index(line(LinPos(1):LinEnd(1)),'print_vacs').ne.0)then
            if(NumLin.lt.2) go to 50
            if(index(line(LinPos(2):LinEnd(2)),'on').ne.0) then
               ipvacs = .true.
            else if(index(line(LinPos(2):LinEnd(2)),'off').ne.0) then
               ipvacs = .false.
            end if
            iw=iw+1
         else if(index(line(LinPos(1):LinEnd(1)),'time_out').ne.0)then
            if(NumLin.lt.2) go to 50
            if(index(line(LinPos(2):LinEnd(2)),'on').ne.0) then
               do_tm = .true.
            else if(index(line(LinPos(2):LinEnd(2)),'off').ne.0) then
               do_tm = .false.
            end if
            iw=iw+1
         else if(index(line(LinPos(1):LinEnd(1)),'trajectory').ne.0)then
            if(NumLin.lt.2) go to 50
            read(Line(LinPos(2):),*,END=50,ERR=50) ntraj
            do_tm = .false.
            iw=iw+1
      else if(index(line(LinPos(1):LinEnd(1)),'traj_interval').ne.0)then
            if(NumLin.lt.2) go to 50
            read(Line(LinPos(2):),*,END=50,ERR=50) tm_int
            do_tm = .true.
            iw=iw+1
      else if(index(line(LinPos(1):LinEnd(1)),'stop_time').ne.0)then
            if(NumLin.lt.2) go to 50
            read(Line(LinPos(2):),*,END=50,ERR=50) stopt
            do_stop = .true.
            iw=iw+1
         else if(index(line(LinPos(1):LinEnd(1)),'system').ne.0)then
            if(NumLin.lt.2) go to 50
            read(Line(LinPos(2):),*,END=50,ERR=50) nsys
            iw=iw+1
         else if(index(line(LinPos(1):LinEnd(1)),'seed').ne.0)then
            if(NumLin.lt.2) go to 50
            read(Line(LinPos(2):),*,END=50,ERR=50) iseed
            iw=iw+1
      else if(index(line(LinPos(1):LinEnd(1)),'add_vacs').ne.0)then
            if(NumLin.lt.3) go to 50
            add_vacs = .true.
            read(Line(LinPos(2):),*,END=50,ERR=50) avac_int
            read(Line(LinPos(3):),*,END=50,ERR=50) num_avac
            iw=iw+1
      else if(index(line(LinPos(1):LinEnd(1)),
     &           'print_threshold').ne.0)then
            if(NumLin.lt.2) go to 50
            read(Line(LinPos(2):),*,END=50,ERR=50) print_threshold
            iw=iw+1
      else if(index(line(LinPos(1):LinEnd(1)),'add_ints').ne.0)then
            if(NumLin.lt.3) go to 50
            add_ints = .true.
            read(Line(LinPos(2):),*,END=50,ERR=50) aint_int
            read(Line(LinPos(3):),*,END=50,ERR=50) num_aint
            iw=iw+1
      else if(index(line(LinPos(1):LinEnd(1)),'print_recomb').ne.0)then
            if(NumLin.lt.2) go to 50
            if(index(line(LinPos(2):LinEnd(2)),'on').ne.0) then
               iprecomb = .true.
            else if(index(line(LinPos(2):LinEnd(2)),'off').ne.0) then
               iprecomb = .false.
            end if
            iw=iw+1
      else if(index(line(LinPos(1):LinEnd(1)),'print_link').ne.0)then
            if(NumLin.lt.2) go to 50
            if(index(line(LinPos(2):LinEnd(2)),'on').ne.0) then
               iplink = .true.
            else if(index(line(LinPos(2):LinEnd(2)),'off').ne.0) then
               iplink = .false.
            end if
            iw=iw+1
      else if(index(line(LinPos(1):LinEnd(1)),'print_loss').ne.0)then
            if(NumLin.lt.2) go to 50
            if(index(line(LinPos(2):LinEnd(2)),'on').ne.0) then
               prnt_loss = .true.
            else if(index(line(LinPos(2):LinEnd(2)),'off').ne.0) then
               prnt_loss = .false.
            end if
            iw=iw+1
      else if(index(line(LinPos(1):LinEnd(1)),'healed_recomb').ne.0)then
            if(NumLin.lt.2) go to 50
            if(index(line(LinPos(2):LinEnd(2)),'on').ne.0) then
               healed_recomb = .true.
            else if(index(line(LinPos(2):LinEnd(2)),'off').ne.0) then
               healed_recomb = .false.
            end if
            iw=iw+1
      else if(index(line(LinPos(1):LinEnd(1)),'do_loops').ne.0)then
            if(NumLin.lt.6) go to 50
            loopl = .true.
            read(Line(LinPos(2):),*,END=50,ERR=50) r_nucl
            read(Line(LinPos(3):),*,END=50,ERR=50) r_init
            read(Line(LinPos(4):),*,END=50,ERR=50) r_capt
            read(Line(LinPos(5):),*,END=50,ERR=50) v_capt
            if(index(line(LinPos(6):LinEnd(6)),'atoms').ne.0) then
               loop_atoms = .true.
            else if(index(line(LinPos(6):LinEnd(6)),'discs').ne.0) then
               loop_atoms = .false.
            end if
            iw=iw+1
         end if
      end do
      
 50   write(*,*)'error reading input file'
      stop

 100  write(*,*)' ... reading input file '
      write(*,*)'     layers: ',nl
      write(*,*)'     steps: ',nsteps
      write(*,*)'     temperature : ',temp
      
      end

      subroutine write_sites(nsites,sites,insites,isites,snsites,
     &        ssites,nstep)
c---------------------------------------------------
c     writing sites to file
c---------------------------------------------------
      integer sites(9999,5)
      integer isites(9999,4)
      integer ssites(9999,5)
      integer nsites,snsites,insites

      write(5,'(a6,i4,a8,i4)')'step: ',nstep,' sites: ',nsites

      do i=1,nsites
         write(5,'(4(a6,i3))')'  x = ',sites(i,1),'  y = ',sites(i,2),
     &        '  z = ',sites(i,3),'  dir ',sites(i,4)
      end do      

      write(5,'(a6,i4,a9,i4)')'step: ',nstep,' isites: ',insites

      do i=1,insites
         write(5,'(4(a6,i3))')'  x = ',isites(i,1),'  y = ',isites(i,2),
     &        '  z = ',isites(i,3),'  dir ',isites(i,4)
      end do            
      
      write(5,'(a6,i4,a9,i4)')'step: ',nstep,' ssites: ',snsites

      do i=1,snsites
         write(5,'(4(a6,i3))')'  x = ',ssites(i,1),'  y = ',ssites(i,2),
     &        '  z = ',ssites(i,3),'  dir ',ssites(i,4)
      end do            

      end


      subroutine write_energies(nsites,energies,insites,ienergies,
     &     snsites,senergies,nstep)
c---------------------------------------------------
c     writing sites to file
c---------------------------------------------------
      real*8 energies(9999,5)
      real*8 ienergies(9999,6)
      real*8 senergies(9999,6)
      integer nsites,insites,snsites

      write(8,'(a6,i4,a6)')'step: ',nstep,' sites'

      do i=1,nsites
         write(8,*)'site ',i
         write(8,*)'E1 = ',energies(i,1),' eV'
         write(8,*)'E2 = ',energies(i,2),' eV'
         write(8,*)'E3 = ',energies(i,3),' eV'
         write(8,*)'E4 = ',energies(i,4),' eV'
         write(8,*)'E5 = ',energies(i,5),' eV'
      end do

      write(8,'(a6,i4,a7)')'step: ',nstep,' isites'

      do i=1,insites
         write(8,*)'site ',i
         write(8,*)'E1 = ',ienergies(i,1),' eV'
         write(8,*)'E2 = ',ienergies(i,2),' eV'
         write(8,*)'E3 = ',ienergies(i,3),' eV'
         write(8,*)'E4 = ',ienergies(i,4),' eV'
         write(8,*)'E5 = ',ienergies(i,5),' eV'
         write(8,*)'E6 = ',ienergies(i,6),' eV'
      end do

      write(8,'(a6,i4,a7)')'step: ',nstep,' ssites'

      do i=1,snsites
         write(8,*)'site ',i
         write(8,*)'E1 = ',senergies(i,1),' eV'
         write(8,*)'E2 = ',senergies(i,2),' eV'
         write(8,*)'E3 = ',senergies(i,3),' eV'
         write(8,*)'E4 = ',senergies(i,4),' eV'
         write(8,*)'E5 = ',senergies(i,5),' eV'
         write(8,*)'E6 = ',senergies(i,6),' eV'
      end do      

      end      



      subroutine makesystem(system,isystem,nat)
c-----------------------------------------------------
c     constructs perfect system matrix               -
c                                                    -
c     interstitial system has 8 at alpha sites      -
c-----------------------------------------------------
      integer*1 system(800,800,199)
      integer*1 isystem(800,800,199)

      logical switch1
      
      switch1 = .true.
      n1 = 0
      nat = 0

      do j=1,700
         n1 = n1 + 1
         if(n1.eq.1 .or. n1.eq.2) then
            do i=1,500
               if(switch1) then
                  system(i,j,1) = 1
                  system(i,j,3) = 1
                  if(j.le.(700-2)) then
                     system(i,j+2,2) = 1
                     system(i,j+2,4) = 1
                  else
                     system(i,j-700,2) = 1
                     system(i,j-700,4) = 1
                  end if
                  nat = nat + 4
                  switch1 = .false.
               else
                  switch1 = .true.
               end if
            end do
            if(switch1) then
               switch1 = .false.
            else
               switch1 = .true.
            end if
         else
            n1 = 0
            if(switch1) then
               switch1 = .false.
            else
               switch1 = .true.
            end if
         end if
      end do

      do j=1,700
         do i=1,500
            do k=1,4
               if(system(i,j,ml(k+1)).eq.1 .and.
     &              system(i,j,k).eq.1) then
                  isystem(i,j,k) = 8
               end if
            end do
         end do
      end do



c      system(12,20,1) = 2
c      system(30,40,1) = 2
c      system(14,18,2) = 2
c      system(19,39,2) = 2

c      system(12,20,3) = 2
c      system(30,40,3) = 2
c      system(33,19,4) = 2
c      system(19,39,4) = 2

      end


      subroutine writesystem(system,isystem,nl,ixd,iyd)
c----------------------------------------------------
c     writes system directly
c----------------------------------------------------
      integer*1 system(800,800,199)
      integer*1 isystem(800,800,199)

      character*6 form1
      character*7 form2

      if(ixd.ge.100) then
         write(form2,'(a1,i3,a3)')"(",ixd,"i1)"
      else
         write(form1,'(a1,i2,a3)')"(",ixd,"i1)"
      end if

      open(71,FILE='system.out')
      
      do k=1,nl
         write(71,'(a6,i2)')'layer ',k   
         do j=iyd,1,-1
            if(ixd.ge.100) then
               write(71,form2) ( system(i,j,k) , i=1,ixd )
            else
               write(71,form1) ( system(i,j,k) , i=1,ixd )
            end if
         end do
         write(71,'(a2,i2)')'i ',k
         do j=iyd,1,-1
            if(ixd.ge.100) then
               write(71,form2) ( isystem(i,j,k) , i=1,ixd )
            else
               write(71,form1) ( isystem(i,j,k) , i=1,ixd )
            end if
         end do
      end do
      
      close(71)

      end

      subroutine readsystem(system,isystem,nl,nat,ixd,iyd)
c----------------------------------------------------
c     read system directly
c----------------------------------------------------
      integer*1 system(800,800,199)
      integer*1 isystem(800,800,199)

      character*6 form1
      character*7 form2

      if(ixd.ge.100) then
         write(form2,'(a1,i3,a3)')"(",ixd,"i1)"
      else
         write(form1,'(a1,i2,a3)')"(",ixd,"i1)"
      end if

      nat = 0

      open(72,FILE='system.in')

      do k=1,199
         read(72,*,end=100)   
         do j=iyd,1,-1
            if(ixd.ge.100) then
               read(72,form2,end=100) ( system(i,j,k) , i=1,ixd )
            else
               read(72,form1,end=100) ( system(i,j,k) , i=1,ixd )
            end if
            do i=1,ixd
               if(system(i,j,k).eq.1 .or. system(i,j,k).eq.5 .or.
     &              system(i,j,k).eq.6) then
                  nat = nat + 1
               else if(system(i,j,k).eq.9) then
                  nat = nat + 2
               end if
            end do
         end do
         read(72,*,end=100)
         do j=iyd,1,-1
            if(ixd.ge.100) then
               read(72,form2,end=100) ( isystem(i,j,k) , i=1,ixd )
            else
               read(72,form1,end=100) ( isystem(i,j,k) , i=1,ixd )
            end if
            do i=1,ixd
               if(isystem(i,j,k).ge.1 .and. isystem(i,j,k).le.6) then
                  nat = nat + 1
               end if
            end do
         end do
      end do

 100  nl1 = k-1

      if(nl1.ne.nl) then
         write(*,*)' error in read system: layer number mismatch '
         stop
      end if

c-----cheking input array for likely problems

      do j=3,iyd-2
         do i=2,ixd-1
            do k=1,nl
               if(isystem(i,j,k).ne.0 .and. 
     &              system(i,j,ml(k+1)).eq.0 .and.
     &              system(i,j,k).eq.0) then
                  write(*,*)'error in reading interstitial plane'
                  write(*,*)i,j,k,system(i,j,k)
                  stop
               end if
               if(isystem(i,j,k).ne.0 .and. 
     &              (system(mb(i+1),j,k).ne.0 .or.
     &              system(i,md(j+1),k).ne.0 .or.
     &              system(mb(i-1),j,k).ne.0 .or.
     &              system(i,md(j-1),k).ne.0)) then
                  write(*,*)'error in reading lattice plane'
                  write(*,*)i,j,k
                  write(*,*)system(i,j,k)
                  write(*,*)system(mb(i-1),j,k),
     &              system(mb(i-1),md(j+1),k),system(i,md(j+1),k),
     &              system(mb(i+1),j,k),system(mb(i+1),md(j-1),k),
     &              system(i,md(j-1),k),system(mb(i-1),md(j-1),k)
                  stop
               end if

            end do
         end do
      end do

c      close(72)

      end
      

      subroutine printsystem(system,isystem,nat,nl,ixd,iyd,per_x,
     &     per_y,per_z)
c----------------------------------------------------
c     writes island to an xyz file
c----------------------------------------------------
      integer*1 system(800,800,199)
      integer*1 isystem(800,800,199)
      logical per_x,per_y,per_z

      write(2,'(i8)')nat
      write(2,'(a1)')' '

      nt = 0

      do k=1,nl
         do i=1,ixd
            do j=1,iyd
               if(system(i,j,k).eq.1) then
                  write(2,'(a2,3(x,f8.4))')'C ',i*1.22,
     &                 j*0.71,k*3.34
                  nt = nt + 1
               else if(system(i,j,k).eq.5) then
                  write(2,'(a2,3(x,f8.4))')'C ',i*1.22,
     &                 j*0.71,k*3.34+0.9
                  nt = nt + 1
               else if(system(i,j,k).eq.6) then
                  write(2,'(a2,3(x,f8.4))')'C ',i*1.22,
     &                 j*0.71,k*3.34-0.9
                  nt = nt + 1
               else if(system(i,j,k).eq.9) then
                  write(2,'(a2,3(x,f8.4))')'C ',i*1.22,
     &                 j*0.71,k*3.34-0.69
                  nt = nt + 1
                  write(2,'(a2,3(x,f8.4))')'C ',i*1.22,
     &                 j*0.71,k*3.34+0.69
                  nt = nt + 1
               end if

c-----------interstitial 

c..............spiro
               if(isystem(i,j,k).eq.1) then
                  write(2,'(a2,3(x,f8.4))')'C ',i*1.22+0.355,
     &                 j*0.71+0.615,k*3.34+1.67
                  nt = nt + 1
               else if(isystem(i,j,k).eq.2) then
                  write(2,'(a2,3(x,f8.4))')'C ',i*1.22+0.71,
     &                 j*0.71,k*3.34+1.67
                  nt = nt + 1
               else if(isystem(i,j,k).eq.3) then
                  write(2,'(a2,3(x,f8.4))')'C ',i*1.22+0.355,
     &                 j*0.71-0.615,k*3.34+1.67
                  nt = nt + 1
               else if(isystem(i,j,k).eq.4) then
                  write(2,'(a2,3(x,f8.4))')'C ',i*1.22-0.355,
     &                 j*0.71-0.615,k*3.34+1.67
                  nt = nt + 1
               else if(isystem(i,j,k).eq.5) then
                  write(2,'(a2,3(x,f8.4))')'C ',i*1.22-0.71,
     &                 j*0.71,k*3.34+1.67
                  nt = nt + 1
               else if(isystem(i,j,k).eq.6) then
                  write(2,'(a2,3(x,f8.4))')'C ',i*1.22-0.355,
     &                 j*0.71+0.615,k*3.34+1.67
                  nt = nt + 1
c..............lower grafted
               else if(isystem(i,j,k).eq.11) then
                  write(2,'(a2,3(x,f8.4))')'C ',i*1.22,
     &                 j*0.71+0.71,k*3.34+1.15
                  nt = nt + 1
               else if(isystem(i,j,k).eq.12) then
                  write(2,'(a2,3(x,f8.4))')'C ',i*1.22+0.615,
     &                 j*0.71+0.355,k*3.34+1.15
                  nt = nt + 1
               else if(isystem(i,j,k).eq.13) then
                  write(2,'(a2,3(x,f8.4))')'C ',i*1.22+0.615,
     &                 j*0.71-0.355,k*3.34+1.15
                  nt = nt + 1
               else if(isystem(i,j,k).eq.14) then
                  write(2,'(a2,3(x,f8.4))')'C ',i*1.22,
     &                 j*0.71-0.71,k*3.34+1.15
                  nt = nt + 1
               else if(isystem(i,j,k).eq.15) then
                  write(2,'(a2,3(x,f8.4))')'C ',i*1.22-0.615,
     &                 j*0.71-0.355,k*3.34+1.15
                  nt = nt + 1
               else if(isystem(i,j,k).eq.16) then
                  write(2,'(a2,3(x,f8.4))')'C ',i*1.22-0.615,
     &                 j*0.71+0.355,k*3.34+1.15
                  nt = nt + 1
c..............upper grafted
               else if(isystem(i,j,k).eq.21) then
                  if(k.ne.nl) then
                     write(2,'(a2,3(x,f8.4))')'C ',i*1.22,
     &                    j*0.71+0.71,k*3.34+2.19
                  else
                     write(2,'(a2,3(x,f8.4))')'C ',i*1.22,
     &                    j*0.71+0.71,2.19
                  end if
                  nt = nt + 1
               else if(isystem(i,j,k).eq.22) then
                  if(k.ne.nl) then
                     write(2,'(a2,3(x,f8.4))')'C ',i*1.22+0.615,
     &                    j*0.71+0.355,k*3.34+2.19
                  else
                      write(2,'(a2,3(x,f8.4))')'C ',i*1.22+0.615,
     &                    j*0.71+0.355,2.19
                  end if
                  nt = nt + 1
               else if(isystem(i,j,k).eq.23) then
                  if(k.ne.nl) then
                     write(2,'(a2,3(x,f8.4))')'C ',i*1.22+0.615,
     &                    j*0.71-0.355,k*3.34+2.19
                  else
                     write(2,'(a2,3(x,f8.4))')'C ',i*1.22+0.615,
     &                    j*0.71-0.355,2.19
                  end if
                  nt = nt + 1
               else if(isystem(i,j,k).eq.24) then
                  if(k.ne.nl) then
                     write(2,'(a2,3(x,f8.4))')'C ',i*1.22,
     &                    j*0.71-0.71,k*3.34+2.19
                  else
                     write(2,'(a2,3(x,f8.4))')'C ',i*1.22,
     &                    j*0.71-0.71,2.19
                  end if
                  nt = nt + 1
               else if(isystem(i,j,k).eq.25) then
                  if(k.ne.nl) then
                     write(2,'(a2,3(x,f8.4))')'C ',i*1.22-0.615,
     &                    j*0.71-0.355,k*3.34+2.19
                  else
                     write(2,'(a2,3(x,f8.4))')'C ',i*1.22-0.615,
     &                    j*0.71-0.355,2.19  
                  end if
                  nt = nt + 1
               else if(isystem(i,j,k).eq.26) then
                  if(k.ne.nl) then
                     write(2,'(a2,3(x,f8.4))')'C ',i*1.22-0.615,
     &                    j*0.71+0.355,k*3.34+2.19
                  else
                     write(2,'(a2,3(x,f8.4))')'C ',i*1.22-0.615,
     &                    j*0.71+0.355,2.19
                  end if
                  nt = nt + 1
               end if

            end do
         end do
      end do

      if(nat.ne.nt) then
         write(*,*)'error in print system - lost atom while periodic'
         write(*,*)nat,nt
         stop
      end if
            
      end


      subroutine initpoints(system,isystem,nat,nl,ixd,iyd,n_vac,
     &     n_i,n_s)
c----------------------------------------------------
c     finds initial number of defects
c----------------------------------------------------
      integer*1 system(800,800,199)
      integer*1 isystem(800,800,199)

      nv = 0
      ns = 0
      ni = 0

      do k=1,nl
         do i=1,ixd
            do j=1,iyd
               if(system(i,j,k).eq.2) then
                  nv = nv + 1
               else if(system(i,j,k).eq.9) then
                  ns = ns + 1
               end if


               if(isystem(i,j,k).ge.1 .and. 
     &              isystem(i,j,k).le.6) then
                  ni = ni + 1
               else if(isystem(i,j,k).ge.11) then
                  ni = ni + 1
               end if
               
            end do
         end do
      end do

      n_vac = nv
      n_i = ni
      n_s = ns
      
      write(*,*)'vacs: ',nv
      write(*,*)'splits: ',ns
      write(*,*)'ints: ',ni

      end

      subroutine printvacs(system,nl,ixd,iyd,n_vac)
c----------------------------------------------------
c     writes island to an xyz file
c----------------------------------------------------
      integer*1 system(800,800,199)

      nv = 0

      do k=1,nl
         do i=1,ixd
            do j=1,iyd
               if(system(i,j,k).eq.2 .or. system(i,j,k).eq.7) then
                  nv = nv + 1
               end if
            end do
         end do
      end do

      write(59,'(i4)')nv
      write(59,'(a1)')' '

      do k=1,nl
         do i=1,ixd
            do j=1,iyd
               if(system(i,j,k).eq.2 .or. system(i,j,k).eq.7) then
                  write(59,'(a2,3(x,f8.4))')'V ',i*1.22,
     &                 j*0.71,k*3.34
               end if
            end do
         end do
      end do

      end


      subroutine printlinks(system,nl,ixd,iyd)
c----------------------------------------------------
c     writes island to an xyz file
c----------------------------------------------------
      integer*1 system(800,800,199)

      n1 = 0

      do k=1,nl
         do i=1,ixd
            do j=1,iyd
               if(system(i,j,k).eq.5) then
                  n1 = n1 + 1
               end if
            end do
         end do
      end do

      write(66,'(i4)')n1
      write(66,'(a1)')' '

      do k=1,nl
         do i=1,ixd
            do j=1,iyd
               if(system(i,j,k).eq.5) then
                  write(66,'(a3,3(x,f8.4))')'Li ',i*1.22,
     &                 j*0.71,k*3.34+1.67
               end if
            end do
         end do
      end do

      end




      subroutine printloops(iloop,ploop,nploop,nl,ixd,iyd,per_x,
     &                   per_y,per_z)
c----------------------------------------------------
c     writes island to an xyz file
c----------------------------------------------------
      real*8 ploop(9999)
      integer iloop(9999,4)
      logical per_x,per_y,per_z
      
      write(76,'(i8)')nploop
      write(76,*)' '
      do i=1,nploop
         write(76,'(4(x,f12.4))')iloop(i,1)*1.22,iloop(i,2)*0.71,
     &        iloop(i,3)*3.35,ploop(i)
      end do

      end

      subroutine printrecom(recom_evnt,t,ix_rc,iy_rc,iz_rc,nrecom)
c----------------------------------------------------
c     writes island to an xyz file
c----------------------------------------------------
      logical recom_evnt
      integer ix_rc(9999)
      integer iy_rc(9999)
      integer iz_rc(9999)
      real*8 t

      if(nrecom.ge.1) then
         write(65,'(i4,f16.8)')nrecom,t
         do i=1,nrecom
            write(65,'(i4,3(x,f12.6))')i,ix_rc(i)*1.22,
     &                 iy_rc(i)*0.71,iz_rc(i)*3.34
         end do
      else
         write(65,'(a2,f16.8)')'0 ',t
      end if

      end

      subroutine printloss(loss_evnt,t,ix_ls,iy_ls,iz_ls,nloss)
c----------------------------------------------------
c     writes island to an xyz file
c----------------------------------------------------
      logical loss_evnt
      integer ix_ls(9999)
      integer iy_ls(9999)
      integer iz_ls(9999)
      real*8 t

      if(nloss.ge.1) then
         write(77,'(i4,f16.8)')nloss,t
         do i=1,nloss
         write(77,'(i4,3(x,f12.6))')i,ix_ls(i)*1.22,
     &                 iy_ls(i)*0.71,iz_ls(i)*3.34
         end do
      else
         write(77,'(a2,f16.8)')'0 ',t
      end if

      end


      subroutine printpoints(system,isystem,nat,nl,ixd,iyd,n_vac,
     &     n_i,n_s,
     &     loop_atoms,iloop,ploop,nploop,per_x,
     &     per_y,per_z,idum)
c----------------------------------------------------
c     writes island to an xyz file
c----------------------------------------------------
      integer*1 system(800,800,199)
      integer*1 isystem(800,800,199)
      real*8 ploop(9999),zero,dist1,x1,y1,x2,y2
      integer iloop(9999,4)
      logical per_x,per_y,per_z,loop_atoms

      zero = 0.0

      nv = 0

      do k=1,nl
         do i=1,ixd
            do j=1,iyd
               if(system(i,j,k).eq.2 .or. system(i,j,k).eq.7) then
                  nv = nv + 1
               end if
            end do
         end do
      end do

      ni = 0

      do k=1,nl
         do i=1,ixd
            do j=1,iyd
               if(isystem(i,j,k).ge.1 .and. isystem(i,j,k).le.6) then
                  ni = ni + 1
        else if(isystem(i,j,k).ge.11 .and. isystem(i,j,k).le.16) then
                  ni = ni + 1
        else if(isystem(i,j,k).ge.21 .and. isystem(i,j,k).le.26) then
                  ni = ni + 1
               end if
            end do
         end do
      end do

      ns = 0

      do k=1,nl
         do i=1,ixd
            do j=1,iyd
               if(system(i,j,k).eq.9) then
                  ns = ns + 1
               end if
            end do
         end do
      end do
      num_lp_at = 0
      if(loop_atoms) then
      do il=1,nploop
         do i=-30,30
            do j=-50,50
               ii = iloop(il,1) + i
               jj = iloop(il,2) + j
               kk = iloop(il,3)
               if((system(mb(ii),md(jj),kk).ne.0 .and. 
     &              system(mb(ii),md(jj),ml(kk+1)).eq.0) .or.
     &              (system(mb(ii),md(jj),kk).eq.0 .and. 
     &              system(mb(ii),md(jj),ml(kk+1)).ne.0)) then
                  if(iloop(il,4).eq.1) then
                     add1 = 0.2
                     add2 = 1.0
                  else if(iloop(il,4).eq.2) then
                     add1 = 0.9
                     add2 = 0.4
                  else if(iloop(il,4).eq.3) then
                     add1 = 0.9
                     add2 = -0.6
                  else if(iloop(il,4).eq.4) then
                     add1 = -0.2
                     add2 = -1.0
                  else if(iloop(il,4).eq.5) then
                     add1 = -0.9
                     add2 = -0.4
                  else if(iloop(il,4).eq.6) then
                     add1 = -0.9
                     add2 = 0.6
                  else
                     write(*,*)'error in mix nucleation site'
                  end if
                  x1 = ii*1.22
                  y1 = jj*0.71
                  x2 = iloop(il,1)*1.22+add1
                  y2 = iloop(il,2)*0.71+add2
                  dist1 = bond_l(x1,y1,zero,x2,y2,zero)
c                  write(*,*)x1,y1,x2,y2
                  if(dist1.le.ploop(il)) then
                     num_lp_at = num_lp_at + 1
                  end if
               end if
            end do
         end do
      end do
      end if

      write(991,*)num_lp_at

      write(2,'(i4)')nv+ni+ns+num_lp_at
      write(2,'(a1)')' '

      nv = 0

      do k=1,nl
         do i=1,ixd
            do j=1,iyd
               if(system(i,j,k).eq.2 .or. system(i,j,k).eq.7) then
                  write(2,'(a2,3(x,f12.4))')'V ',i*1.22,
     &                 j*0.71,k*3.34
               end if
            end do
         end do
      end do

      do k=1,nl
         do i=1,ixd
            do j=1,iyd
               if(isystem(i,j,k).eq.1) then
                  write(2,'(a2,3(x,f12.4))')'C ',i*1.22+0.355,
     &                 j*0.71+0.615,k*3.34+1.67
               else if(isystem(i,j,k).eq.2) then
                  write(2,'(a2,3(x,f12.4))')'C ',i*1.22+0.71,
     &                 j*0.71,k*3.34+1.67
               else if(isystem(i,j,k).eq.3) then
                  write(2,'(a2,3(x,f12.4))')'C ',i*1.22+0.355,
     &                 j*0.71-0.615,k*3.34+1.67
               else if(isystem(i,j,k).eq.4) then
                  write(2,'(a2,3(x,f12.4))')'C ',i*1.22-0.355,
     &                 j*0.71-0.615,k*3.34+1.67
               else if(isystem(i,j,k).eq.5) then
                  write(2,'(a2,3(x,f12.4))')'C ',i*1.22-0.71,
     &                 j*0.71,k*3.34+1.67
               else if(isystem(i,j,k).eq.6) then
                  write(2,'(a2,3(x,f12.4))')'C ',i*1.22-0.355,
     &                 j*0.71+0.615,k*3.34+1.67
c..............lower grafted
               else if(isystem(i,j,k).eq.11) then
                  write(2,'(a2,3(x,f12.4))')'C ',i*1.22,
     &                 j*0.71+0.71,k*3.34+1.15
               else if(isystem(i,j,k).eq.12) then
                  write(2,'(a2,3(x,f12.4))')'C ',i*1.22+0.615,
     &                 j*0.71+0.355,k*3.34+1.15
               else if(isystem(i,j,k).eq.13) then
                  write(2,'(a2,3(x,f12.4))')'C ',i*1.22+0.615,
     &                 j*0.71-0.355,k*3.34+1.15
               else if(isystem(i,j,k).eq.14) then
                  write(2,'(a2,3(x,f12.4))')'C ',i*1.22,
     &                 j*0.71-0.71,k*3.34+1.15
               else if(isystem(i,j,k).eq.15) then
                  write(2,'(a2,3(x,f12.4))')'C ',i*1.22-0.615,
     &                 j*0.71-0.355,k*3.34+1.15
               else if(isystem(i,j,k).eq.16) then
                  write(2,'(a2,3(x,f12.4))')'C ',i*1.22-0.615,
     &                 j*0.71+0.355,k*3.34+1.15
c..............upper grafted
               else if(isystem(i,j,k).eq.21) then
                  if(k.ne.nl) then
                     write(2,'(a2,3(x,f12.4))')'C ',i*1.22,
     &                    j*0.71+0.71,k*3.34+2.19
                  else
                     write(2,'(a2,3(x,f12.4))')'C ',i*1.22,
     &                    j*0.71+0.71,2.19
                  end if
               else if(isystem(i,j,k).eq.22) then
                  if(k.ne.nl) then
                     write(2,'(a2,3(x,f12.4))')'C ',i*1.22+0.615,
     &                    j*0.71+0.355,k*3.34+2.19
                  else
                      write(2,'(a2,3(x,f12.4))')'C ',i*1.22+0.615,
     &                    j*0.71+0.355,2.19
                  end if
               else if(isystem(i,j,k).eq.23) then
                  if(k.ne.nl) then
                     write(2,'(a2,3(x,f12.4))')'C ',i*1.22+0.615,
     &                    j*0.71-0.355,k*3.34+2.19
                  else
                     write(2,'(a2,3(x,f12.4))')'C ',i*1.22+0.615,
     &                    j*0.71-0.355,2.19
                  end if
               else if(isystem(i,j,k).eq.24) then
                  if(k.ne.nl) then
                     write(2,'(a2,3(x,f12.4))')'C ',i*1.22,
     &                    j*0.71-0.71,k*3.34+2.19
                  else
                     write(2,'(a2,3(x,f12.4))')'C ',i*1.22,
     &                    j*0.71-0.71,2.19
                  end if
               else if(isystem(i,j,k).eq.25) then
                  if(k.ne.nl) then
                     write(2,'(a2,3(x,f12.4))')'C ',i*1.22-0.615,
     &                    j*0.71-0.355,k*3.34+2.19
                  else
                     write(2,'(a2,3(x,f12.4))')'C ',i*1.22-0.615,
     &                    j*0.71-0.355,2.19  
                  end if
               else if(isystem(i,j,k).eq.26) then
                  if(k.ne.nl) then
                     write(2,'(a2,3(x,f12.4))')'C ',i*1.22-0.615,
     &                    j*0.71+0.355,k*3.34+2.19
                  else
                     write(2,'(a2,3(x,f12.4))')'C ',i*1.22-0.615,
     &                    j*0.71+0.355,2.19
                  end if
               end if

            end do
         end do
      end do

      ns = 0

      do k=1,nl
         do i=1,ixd
            do j=1,iyd
               if(system(i,j,k).eq.9) then
                  write(2,'(a2,3(x,f12.4))')'C ',i*1.22,
     &                 j*0.71,k*3.34
               end if
            end do

         end do
      end do


c-----write interstitial loops
      num_lp_at2 = 0
      if(loop_atoms) then
      do il=1,nploop
         do i=-30,30
            do j=-50,50
               ii = iloop(il,1) + i
               jj = iloop(il,2) + j
               kk = iloop(il,3)
               if((system(mb(ii),md(jj),kk).ne.0 .and. 
     &              system(mb(ii),md(jj),ml(kk+1)).eq.0) .or.
     &              (system(mb(ii),md(jj),kk).eq.0 .and. 
     &              system(mb(ii),md(jj),ml(kk+1)).ne.0)) then
                  if(iloop(il,4).eq.1) then
                     add1 = 0.2
                     add2 = 1.0
                  else if(iloop(il,4).eq.2) then
                     add1 = 0.9
                     add2 = 0.4
                  else if(iloop(il,4).eq.3) then
                     add1 = 0.9
                     add2 = -0.6
                  else if(iloop(il,4).eq.4) then
                     add1 = -0.2
                     add2 = -1.0
                  else if(iloop(il,4).eq.5) then
                     add1 = -0.9
                     add2 = -0.4
                  else if(iloop(il,4).eq.6) then
                     add1 = -0.9
                     add2 = 0.6 
                  else
                     write(*,*)'error in mix nucleation site 2'
                  end if
                 
                  x1 = ii*1.22
                  y1 = jj*0.71
                  x2 = iloop(il,1)*1.22+add1
                  y2 = iloop(il,2)*0.71+add2
                  dist1 = bond_l(x1,y1,zero,x2,y2,zero)
                  if(dist1.le.ploop(il)) then
                     write(2,'(a2,3(x,f12.4))')'I ',mb(ii)*1.22,
     &                    md(jj)*0.71,ml(kk)*3.34+1.67
                     num_lp_at2 = num_lp_at2 + 1
                  end if
               end if
            end do
         end do
      end do

      if(num_lp_at.ne.num_lp_at2) then
         write(*,*)'error with write loop atoms',num_lp_at,num_lp_at2
         stop
      end if

      end if

      end
      

      subroutine readens(etab,stab,ntab,mtab,emin)
c----------------------------------------------------------
c     subroutine to read in the library of structures and -
c     energeis from the ens.dat file                      -
c                                                         -
c     the ens.dat file contains the in plane structure    -
c     schematically in a 33x16 array, where               -
c     o : carbon atom position                            -
c     - : vacancy position                                -
c     + : atom bonded to adjacent layer                   -
c     1-6 : spiro-interstitial                            -
c     a-f : grafted interstitial                          -
c     A-F : spiro directly above vac in target plane      -
c     s : split interstitial                              -
c     * : wildcard                                        -
c     # : interstitial wildcard (no vacancy)              -
c     @ : spiro interstitial of any orientation           -
c     which is converted to a 17x24 array for the code    -
c                                                         -
c     the position of the active site is at 17,9 (9,13)   -
c     (denoted with O)                                    -
c                                                         -
c     the first line is the in-plane barrier              -
c     the second is to either make or break a ILDV1,      -
c     ILDV2 and ILDV3 bond with unsaturated carbon        - 
c     the third is to either make or break a ILDV1,       -
c     ILDV2 and ILDV3 bond with saturated carbon          - 
c     the fourth is to move the atom between the planes   -
c                                                         -
c     interstitial atom implies it sits on an occupied    -
c     lattice site                                        -
c----------------------------------------------------------
      character*1 cmat(41,20)
      integer*2 smat(21,30),fmat(21,30)
      integer*2 stab(21,30,12,499)
      real*8 etab(12,499),mtab(12,499)
      real*8 x1(999),y1(999),x2(999),y2(999)
      real*8 x3(999),y3(999),x4(999),y4(999)
      real*8 x5(999),y5(999),x6(999),y6(999)
      real*8 xf1(999),yf1(999),xf2(999),yf2(999)
      real*8 xf3(999),yf3(999),xf4(999),yf4(999)
      real*8 xf5(999),yf5(999),xf6(999),yf6(999)
      integer*2 a1(499),f1(499)
      logical sw1,sw2,emin

      ntab = 0

      open(21,FILE='ens.dat')
c      open(32,FILE='coords1.xyz')

c-----initialise the library tables

      do k=1,499
         do i=1,21
            do j=1,30
               do n=1,12
                  stab(i,j,n,k) = 3
               end do
            end do
         end do
         do i=1,12
            etab(i,k) = 0.0
            mtab(i,k) = 0.0
         end do
      end do

      write(*,*)'reading ens.dat file'

      do nrec=1,499

         do i=1,21
            do j=1,30
               smat(i,j) = 0
            end do
         end do

c--------read in structure
         
         do j=1,20 
            read(21,'(41a1)',end=100) ( cmat(i,j), i=1,41 ) 
         end do

c         do j=1,16
c            write(*,'(33a1)') ( cmat(i,j), i=1,33 )
c         end do
         
         ii = 1
         sw1 = .true.
         jj = 1
         ji = 2

c--------convert to internal format

         do j=20,1,-1
            do i=1,41
               if(sw1) then 
                  
                  if(cmat(i,j).eq.' ') then
                     smat(ii,jj) = 0
                  else if(cmat(i,j).eq.'o' .or. cmat(i,j).eq.'O') then
                     smat(ii,jj) = 1
                  else if(cmat(i,j).eq.'-') then
                     smat(ii,jj) = 2
                  else if(cmat(i,j).eq.'*') then
                     smat(ii,jj) = 3
                  else if(cmat(i,j).eq.'+') then
                     smat(ii,jj) = 4
                  else if(cmat(i,j).eq.'1') then
                     smat(ii,jj) = 11
                  else if(cmat(i,j).eq.'2') then
                     smat(ii,jj) = 12
                  else if(cmat(i,j).eq.'3') then
                     smat(ii,jj) = 13
                  else if(cmat(i,j).eq.'4') then
                     smat(ii,jj) = 14
                  else if(cmat(i,j).eq.'5') then
                     smat(ii,jj) = 15
                  else if(cmat(i,j).eq.'6') then
                     smat(ii,jj) = 16
                  else if(cmat(i,j).eq.'a') then
                     smat(ii,jj) = 21
                  else if(cmat(i,j).eq.'b') then
                     smat(ii,jj) = 22
                  else if(cmat(i,j).eq.'c') then
                     smat(ii,jj) = 23
                  else if(cmat(i,j).eq.'d') then
                     smat(ii,jj) = 24
                  else if(cmat(i,j).eq.'e') then
                     smat(ii,jj) = 25
                  else if(cmat(i,j).eq.'f') then
                     smat(ii,jj) = 26
                  else if(cmat(i,j).eq.'A') then
                     smat(ii,jj) = 31
                  else if(cmat(i,j).eq.'B') then
                     smat(ii,jj) = 32
                  else if(cmat(i,j).eq.'C') then
                     smat(ii,jj) = 33
                  else if(cmat(i,j).eq.'D') then
                     smat(ii,jj) = 34
                  else if(cmat(i,j).eq.'E') then
                     smat(ii,jj) = 35
                  else if(cmat(i,j).eq.'F') then
                     smat(ii,jj) = 36
                  else if(cmat(i,j).eq.'s') then
                     smat(ii,jj) = 9
                  else if(cmat(i,j).eq.'#') then
                     smat(ii,jj) = 17
                  else if(cmat(i,j).eq.'@') then
                     smat(ii,jj) = 18
                  else
                     write(*,*)'Non-valid character in ens.dat'
                     write(*,*)cmat(i,j),i,j
                     stop
                  end if

                  ii = ii + 1
                  sw1 = .false.
                  
               else
                  sw1 = .true.
               end if
               
            end do
            
            ii = 1
            sw1 = .true.
            
            if(ji.le.1) then
               jj = jj + 1
               ji = ji + 1
            else
               jj = jj + 2
               ji = 1
            end if
            
         end do



c-----store in the structure library table

         do j=1,30
            do i=1,21
               fmat(i,j) = smat(i,j)
            end do
         end do         

         call flip(fmat)

        do j=1,30
            do i=1,21
               if(fmat(i,j).eq.11) then
                  fmat(i,j) = 16
               else if(fmat(i,j).eq.12) then 
                  fmat(i,j) = 15
               else if(fmat(i,j).eq.13) then
                  fmat(i,j) = 14
               else if(fmat(i,j).eq.14) then
                  fmat(i,j) = 13
               else if(fmat(i,j).eq.15) then 
                  fmat(i,j) = 12
               else if(fmat(i,j).eq.16) then
                  fmat(i,j) = 11
               else if(fmat(i,j).eq.21) then 
                  fmat(i,j) = 26
               else if(fmat(i,j).eq.22) then 
                  fmat(i,j) = 25
               else if(fmat(i,j).eq.23) then 
                  fmat(i,j) = 24
               else if(fmat(i,j).eq.24) then
                  fmat(i,j) = 23
               else if(fmat(i,j).eq.25) then
                  fmat(i,j) = 22
               else if(fmat(i,j).eq.26) then
                  fmat(i,j) = 21
               else if(fmat(i,j).eq.31) then 
                  fmat(i,j) = 36
               else if(fmat(i,j).eq.32) then 
                  fmat(i,j) = 35
               else if(fmat(i,j).eq.33) then 
                  fmat(i,j) = 34
               else if(fmat(i,j).eq.34) then
                  fmat(i,j) = 33
               else if(fmat(i,j).eq.35) then
                  fmat(i,j) = 32
               else if(fmat(i,j).eq.36) then
                  fmat(i,j) = 31
               end if
            end do
         end do      

         natom = 0

         do i=1,21
            do j=1,30
               if(smat(i,j).ge.1) then
                  natom = natom + 1
                  x1(natom) = i*1.23 - 13.5300
                  y1(natom) = j*0.71 - 11.36
                  a1(natom) = smat(i,j)
               end if
            end do
         end do

         nfatom = 0

         do i=1,21
            do j=1,30
               if(fmat(i,j).ge.1) then
                  nfatom = nfatom + 1
                  xf1(nfatom) = i*1.23 - 13.5300
                  yf1(nfatom) = j*0.71 - 11.36
                  f1(nfatom) = fmat(i,j)
               end if
            end do
         end do

         call round(smat)
         call round(fmat)

         do j=1,30
            do i=1,21
               stab(i,j,1,nrec) = smat(i,j)
            end do
         end do


         do j=1,30
            do i=1,21
               stab(i,j,2,nrec) = fmat(i,j)
            end do
         end do


c........rotate 60 degrees

         do i=1,natom
            x2(i)  = 0.5*x1(i) + 0.866025*y1(i)
            y2(i) = - 0.866025*x1(i) + 0.5*y1(i)
            xf2(i)  = 0.5*xf1(i) + 0.866025*yf1(i)
            yf2(i) = - 0.866025*xf1(i) + 0.5*yf1(i)
         end do

         do i=1,21
            do j=1,30
               smat(i,j) = 0
               fmat(i,j) = 0
            end do
         end do
         
c........read back into smat array

         do i=1,21
            do j=1,30
               xt = i*1.23 - 13.5300
               yt = j*0.71 - 11.3600
               
               do n=1,natom
                  if(x2(n).gt.(xt-0.2) .and. x2(n).lt.(xt+0.2) .and.
     &                 y2(n).gt.(yt-0.2) .and. y2(n).lt.(yt+0.2))then
                     smat(i,j) = a1(n)
                  end if
               end do

               do n=1,natom
                  if(xf2(n).gt.(xt-0.2) .and. xf2(n).lt.(xt+0.2) .and.
     &                 yf2(n).gt.(yt-0.2) .and. yf2(n).lt.(yt+0.2))then
                     fmat(i,j) = f1(n)
                  end if
               end do
            end do
         end do         
         
c-----store in the structure library table

         call round(smat)
         call round(fmat)

         do j=1,30
            do i=1,21
               if(fmat(i,j).ge.11 .and. fmat(i,j).le.16) then
                  fmat(i,j) = fmat(i,j) + 1
                  if(fmat(i,j).gt.16) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.21 .and. fmat(i,j).le.26) then
                  fmat(i,j) = fmat(i,j) + 1
                  if(fmat(i,j).gt.26) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.31 .and. fmat(i,j).le.36) then
                  fmat(i,j) = fmat(i,j) + 1
                  if(fmat(i,j).gt.36) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if

               if(smat(i,j).ge.11 .and. smat(i,j).le.16) then
                  smat(i,j) = smat(i,j) + 1
                  if(smat(i,j).gt.16) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.21 .and. smat(i,j).le.26) then
                  smat(i,j) = smat(i,j) + 1
                  if(smat(i,j).gt.26) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.31 .and. smat(i,j).le.36) then
                  smat(i,j) = smat(i,j) + 1
                  if(smat(i,j).gt.36) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
            end do
         end do      

         do j=1,30
            do i=1,21
               stab(i,j,3,nrec) = smat(i,j)
            end do
         end do

c........flip and store

         do j=1,30
            do i=1,21
               stab(i,j,4,nrec) = fmat(i,j)
            end do
         end do
        
c........rotate 60 degrees

         do i=1,natom
            x3(i)  = 0.5*x2(i) + 0.866025*y2(i)
            y3(i) = - 0.866025*x2(i) + 0.5*y2(i)
            xf3(i)  = 0.5*xf2(i) + 0.866025*yf2(i)
            yf3(i) = - 0.866025*xf2(i) + 0.5*yf2(i)
         end do

         do i=1,21
            do j=1,30
               smat(i,j) = 0
               fmat(i,j) = 0
            end do
         end do
         
c........read back into smat array

         do i=1,21
            do j=1,30
               xt = i*1.23 - 13.5300
               yt = j*0.71 - 11.3600
               
               do n=1,natom
                  if(x3(n).gt.(xt-0.2) .and. x3(n).lt.(xt+0.2) .and.
     &                 y3(n).gt.(yt-0.2) .and. y3(n).lt.(yt+0.2))then
                     smat(i,j) = a1(n)
                  end if
               end do
               do n=1,natom
                  if(xf3(n).gt.(xt-0.2) .and. xf3(n).lt.(xt+0.2) .and.
     &                 yf3(n).gt.(yt-0.2) .and. yf3(n).lt.(yt+0.2))then
                     fmat(i,j) = f1(n)
                  end if
               end do
            end do
         end do         
         

c-----store in the structure library table

         call round(smat)
         call round(fmat)

         do j=1,30
            do i=1,21
               if(fmat(i,j).ge.11 .and. fmat(i,j).le.16) then
                  fmat(i,j) = fmat(i,j) + 2
                  if(fmat(i,j).gt.16) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.21 .and. fmat(i,j).le.26) then
                  fmat(i,j) = fmat(i,j) + 2
                  if(fmat(i,j).gt.26) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.31 .and. fmat(i,j).le.36) then
                  fmat(i,j) = fmat(i,j) + 2
                  if(fmat(i,j).gt.36) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if

               if(smat(i,j).ge.11 .and. smat(i,j).le.16) then
                  smat(i,j) = smat(i,j) + 2
                  if(smat(i,j).gt.16) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.21 .and. smat(i,j).le.26) then
                  smat(i,j) = smat(i,j) + 2
                  if(smat(i,j).gt.26) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.31 .and. smat(i,j).le.36) then
                  smat(i,j) = smat(i,j) + 2
                  if(smat(i,j).gt.36) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
            end do
         end do      

         do j=1,30
            do i=1,21
               stab(i,j,5,nrec) = smat(i,j)
            end do
         end do

c........flip and store

         do j=1,30
            do i=1,21
               stab(i,j,6,nrec) = fmat(i,j)
            end do
         end do

c........rotate 60 degrees

         do i=1,natom
            x4(i)  = 0.5*x3(i) + 0.866025*y3(i)
            y4(i) = - 0.866025*x3(i) + 0.5*y3(i)
            xf4(i)  = 0.5*xf3(i) + 0.866025*yf3(i)
            yf4(i) = - 0.866025*xf3(i) + 0.5*yf3(i)
         end do

         do i=1,21
            do j=1,30
               smat(i,j) = 0
               fmat(i,j) = 0
            end do
         end do
         
c........read back into smat array

         do i=1,21
            do j=1,30
               xt = i*1.23 - 13.5300
               yt = j*0.71 - 11.3600
               
               do n=1,natom
                  if(x4(n).gt.(xt-0.2) .and. x4(n).lt.(xt+0.2) .and.
     &                 y4(n).gt.(yt-0.2) .and. y4(n).lt.(yt+0.2))then
                     smat(i,j) = a1(n)
                  end if
               end do
               do n=1,natom
                  if(xf4(n).gt.(xt-0.2) .and. xf4(n).lt.(xt+0.2) .and.
     &                 yf4(n).gt.(yt-0.2) .and. yf4(n).lt.(yt+0.2))then
                     fmat(i,j) = f1(n)
                  end if
               end do
            end do
         end do         
         
c-----store in the structure library table

         call round(smat)
         call round(fmat)

         do j=1,30
            do i=1,21
               if(fmat(i,j).ge.11 .and. fmat(i,j).le.16) then
                  fmat(i,j) = fmat(i,j) + 3
                  if(fmat(i,j).gt.16) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.21 .and. fmat(i,j).le.26) then
                  fmat(i,j) = fmat(i,j) + 3
                  if(fmat(i,j).gt.26) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.31 .and. fmat(i,j).le.36) then
                  fmat(i,j) = fmat(i,j) + 3
                  if(fmat(i,j).gt.36) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if

               if(smat(i,j).ge.11 .and. smat(i,j).le.16) then
                  smat(i,j) = smat(i,j) + 3
                  if(smat(i,j).gt.16) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.21 .and. smat(i,j).le.26) then
                  smat(i,j) = smat(i,j) + 3
                  if(smat(i,j).gt.26) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.31 .and. smat(i,j).le.36) then
                  smat(i,j) = smat(i,j) + 3
                  if(smat(i,j).gt.36) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
            end do
         end do      

         do j=1,30
            do i=1,21
               stab(i,j,7,nrec) = smat(i,j)
            end do
         end do

         do j=1,30
            do i=1,21
               stab(i,j,8,nrec) = fmat(i,j)
            end do
         end do

c........rotate 60 degrees

         do i=1,natom
            x5(i)  = 0.5*x4(i) + 0.866025*y4(i)
            y5(i) = - 0.866025*x4(i) + 0.5*y4(i)
            xf5(i)  = 0.5*xf4(i) + 0.866025*yf4(i)
            yf5(i) = - 0.866025*xf4(i) + 0.5*yf4(i)
         end do

         do i=1,21
            do j=1,30
               smat(i,j) = 0
               fmat(i,j) = 0
            end do
         end do
         
c........read back into smat array

         do i=1,21
            do j=1,30
               xt = i*1.23 - 13.5300
               yt = j*0.71 - 11.3600
               
               do n=1,natom
                  if(x5(n).gt.(xt-0.2) .and. x5(n).lt.(xt+0.2) .and.
     &                 y5(n).gt.(yt-0.2) .and. y5(n).lt.(yt+0.2))then
                     smat(i,j) = a1(n)
                  end if
               end do
               do n=1,natom
                  if(xf5(n).gt.(xt-0.2) .and. xf5(n).lt.(xt+0.2) .and.
     &                 yf5(n).gt.(yt-0.2) .and. yf5(n).lt.(yt+0.2))then
                     fmat(i,j) = f1(n)
                  end if
               end do
            end do
         end do         
         
c-----store in the structure library table

         call round(smat)
         call round(fmat)

         do j=1,30
            do i=1,21
               if(fmat(i,j).ge.11 .and. fmat(i,j).le.16) then
                  fmat(i,j) = fmat(i,j) + 4
                  if(fmat(i,j).gt.16) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.21 .and. fmat(i,j).le.26) then
                  fmat(i,j) = fmat(i,j) + 4
                  if(fmat(i,j).gt.26) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.31 .and. fmat(i,j).le.36) then
                  fmat(i,j) = fmat(i,j) + 4
                  if(fmat(i,j).gt.36) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if

               if(smat(i,j).ge.11 .and. smat(i,j).le.16) then
                  smat(i,j) = smat(i,j) + 4
                  if(smat(i,j).gt.16) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.21 .and. smat(i,j).le.26) then
                  smat(i,j) = smat(i,j) + 4
                  if(smat(i,j).gt.26) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.31 .and. smat(i,j).le.36) then
                  smat(i,j) = smat(i,j) + 4
                  if(smat(i,j).gt.36) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
            end do
         end do      

         do j=1,30
            do i=1,21
               stab(i,j,9,nrec) = smat(i,j)
            end do
         end do

         do j=1,30
            do i=1,21
               stab(i,j,10,nrec) = fmat(i,j)
            end do
         end do

c........rotate 60 degrees

         do i=1,natom
            x6(i)  = 0.5*x5(i) + 0.866025*y5(i)
            y6(i) = - 0.866025*x5(i) + 0.5*y5(i)
            xf6(i)  = 0.5*xf5(i) + 0.866025*yf5(i)
            yf6(i) = - 0.866025*xf5(i) + 0.5*yf5(i)
         end do

         do i=1,21
            do j=1,30
               smat(i,j) = 0
               fmat(i,j) = 0
            end do
         end do
         
c........read back into smat array

         do i=1,21
            do j=1,30
               xt = i*1.23 - 13.5300
               yt = j*0.71 - 11.3600
               
               do n=1,natom
                  if(x6(n).gt.(xt-0.2) .and. x6(n).lt.(xt+0.2) .and.
     &                 y6(n).gt.(yt-0.2) .and. y6(n).lt.(yt+0.2))then
                     smat(i,j) = a1(n)
                  end if
               end do
               do n=1,natom
                  if(xf6(n).gt.(xt-0.2) .and. xf6(n).lt.(xt+0.2) .and.
     &                 yf6(n).gt.(yt-0.2) .and. yf6(n).lt.(yt+0.2))then
                     fmat(i,j) = f1(n)
                  end if
               end do
            end do
         end do         
         
c-----store in the structure library table

         call round(smat)
         call round(fmat)

         do j=1,30
            do i=1,21
               if(fmat(i,j).ge.11 .and. fmat(i,j).le.16) then
                  fmat(i,j) = fmat(i,j) + 5
                  if(fmat(i,j).gt.16) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.21 .and. fmat(i,j).le.26) then
                  fmat(i,j) = fmat(i,j) + 5
                  if(fmat(i,j).gt.26) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.31 .and. fmat(i,j).le.36) then
                  fmat(i,j) = fmat(i,j) + 5
                  if(fmat(i,j).gt.36) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if

               if(smat(i,j).ge.11 .and. smat(i,j).le.16) then
                  smat(i,j) = smat(i,j) + 5
                  if(smat(i,j).gt.16) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.21 .and. smat(i,j).le.26) then
                  smat(i,j) = smat(i,j) + 5
                  if(smat(i,j).gt.26) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.31 .and. smat(i,j).le.36) then
                  smat(i,j) = smat(i,j) + 5
                  if(smat(i,j).gt.36) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
            end do
         end do      

         do j=1,30
            do i=1,21
               stab(i,j,11,nrec) = smat(i,j)
            end do
         end do

         do j=1,30
            do i=1,21
               stab(i,j,12,nrec) = fmat(i,j)
            end do
         end do

c-----read in energies

         read(21,*,end=100)etab(1,nrec)
         read(21,*,end=100)etab(4,nrec),etab(5,nrec),etab(6,nrec)
         read(21,*,end=100)etab(7,nrec),etab(8,nrec),etab(9,nrec)
         read(21,*,end=100)etab(10,nrec),etab(11,nrec),etab(12,nrec)

         if(emin) then
            read(21,*,end=100)mtab(1,nrec)
            read(21,*,end=100)mtab(4,nrec),mtab(5,nrec),mtab(6,nrec)
            read(21,*,end=100)mtab(7,nrec),mtab(8,nrec),mtab(9,nrec)
            read(21,*,end=100)mtab(10,nrec),mtab(11,nrec),mtab(12,nrec)
         end if
         
      end do

      write(*,*)'Too many structures - increase tab arrays'
      stop

 100  ntab = nrec - 1

      write(*,*)'found ',ntab,' structures'

c......writing coords

c      do j=24,1,-1
c         write(*,'(17i1)') ( stab(i,j,5,1) , i=1,17 )
c      end do

c      do j=1,24
c         do i=1,17
c            smat(i,j) = stab(i,j,12,1)
c            if(abs((i-9)*(j-13)).ge.50) then
c               smat(i,j) = 3
c            end if
c            
c         end do
c      end do

c      do j=24,1,-1
c         write(*,'(17i1)') ( smat(i,j) , i=1,17 )
c      end do

c      do i=1,17
c         do j=1,24
c            do n=1,11,2
c               write(*,*)n
c            if(stab(i,j,n,1).eq.1) then
c               write(32,'(a2,3(x,f8.4))')'C ',i*1.23 -  11.0700,
c     &              j*0.71 - 9.2300,0.0
c            else if(stab(i,j,n,1).eq.2) then
c               write(32,'(a2,3(x,f8.4))')'B ',i*1.23 -  11.0700,
c     &              j*0.71 - 9.2300,0.0
c            else if(stab(i,j,n,1).eq.3) then
c               write(32,'(a2,3(x,f8.4))')'Si',i*1.23 -  11.0700,
c     &              j*0.71 - 9.2300,0.0
c            else if(stab(i,j,n,1).eq.4) then
c               write(32,'(a2,3(x,f8.4))')'O ',i*1.23 -  11.0700,
c     &              j*0.71 - 9.2300,0.0
c            end if
c            end do
c         end do
c      end do
      
      end 


      subroutine flip(fmat)
c-----------------------------------------------------------
c     flip the pattern matrix left-to-right                -
c-----------------------------------------------------------
      integer*2 fmat(21,30)
      integer*2 tmat(21,30)
  
      do i=1,21
         do j=1,30
            tmat(22-i,j) = fmat(i,j)
         end do
      end do

      do i=1,21
         do j=1,30
            fmat(i,j) =  tmat(i,j)
         end do
      end do
      
      end


      subroutine round(tmat)
c-----------------------------------------------------------
c     round the corners of the matrix with wildcards       -
c-----------------------------------------------------------
      integer*2 tmat(21,30)
  
      do i=1,21
         do j=1,30
            if(j.eq.1 .or. j.eq.2 .or. j.eq.3) then
               tmat(i,j) = 3
            else if(j.le.6 .or. j.ge.28) then
               if(i.le.5 .or. i.ge.17) then
                  tmat(i,j) = 3
               end if
            else if(j.le.9 .or. j.ge.25) then
               if(i.le.4 .or. i.ge.18) then
                  tmat(i,j) = 3
               end if               
            else if(j.le.12 .or. j.ge.22) then
               if(i.le.3 .or. i.ge.19) then
                  tmat(i,j) = 3
               end if
            else if(i.eq.1 .or. i.eq.2 .or. i.eq.20 .or. i.eq.21) then
               tmat(i,j) = 3
            end if
         end do
      end do
      
      end


      subroutine roundi(tmat)
c-----------------------------------------------------------
c     round the corners of the matrix with wildcards       -
c-----------------------------------------------------------
      integer*2 tmat(21,30)
  
      do i=1,21
         do j=1,30
            if(j.eq.1 .or. j.eq.2 .or. j.eq.3) then
               tmat(i,j) = 17
            else if(j.le.6 .or. j.ge.28) then
               if(i.le.5 .or. i.ge.17) then
                  tmat(i,j) = 17
               end if
            else if(j.le.9 .or. j.ge.25) then
               if(i.le.4 .or. i.ge.18) then
                  tmat(i,j) = 17
               end if               
            else if(j.le.12 .or. j.ge.22) then
               if(i.le.3 .or. i.ge.19) then
                  tmat(i,j) = 17
               end if
            else if(i.eq.1 .or. i.eq.2 .or. i.eq.20 .or. i.eq.21) then
               tmat(i,j) = 17
            end if
         end do
      end do
      
      end



      subroutine readiens(ietab,istab,intab,imtab,emin)
c----------------------------------------------------------
c     subroutine to read in the library of structures and -
c     energeis for the interstitail potential energy      -
c     surface from the iens.dat file                      -
c                                                         -
c     the iens.dat file contains the interstitial plane   -
c     structure schematically in a 33x16 array, where     -
c     o : lattice atom position in the adjacent plane     -
c     - : vacancy position                                -
c     + : atom bonded to adjacent layer                   -
c     1-6 : spiro-interstitial                            -
c     A-F: interstitial directly above vac in target plane- 
c     a-f : grafted interstitial                          -
c     s : split interstitial                              -
c     * : wildcard                                        -
c     # : interstitial wildcard (no vacancy, but maybe    - 
c     interstitial)                                       -
c     @ : interstitial of any orientation
c     which is converted to a 17x24 array for the code    -
c                                                         -
c     the position of the active site is at 17,9 (9,13)   -
c                                                         -
c     if the interstitial is a spiro:                     -
c                                                         -
c     the 1st energy is the barrier to rotate the         -
c     interstital clockwise                               -
c     the 2nd energy is the barrier to rotate             -
c     anti-clockwise                                      -
c     the 3rd energy is the barrier to move to the        -
c     nearest alpha spiro                                 -
c     the 4th energy is the barrier to move between spiro -
c     and grafted of the same orientation toward the      -
c     target plane                                        -
c     the 5th energy is the barrier to move between spiro -
c     and grafted of the same orientation away from the   -
c     target plane                                        -
c     the 6th is to close into the nearest vacancy        -
c                                                         -
c     if the interstital is a grafted:                    -
c                                                         -
c     the 1st energy is the barrier to rotate the         -
c     interstital clockwise                               -
c     the 2nd energy is the barrier to rotate             -
c     anti-clockwise                                      -
c     the 3rd energy is the barrier to move to the        -
c     nearest alpha grafted                               -
c     the 4th energy is the barrier to move to spiro of   -
c     the same orientation                                -
c     the 5th energy is the barrier to move to the local  -
c     split                                               -
c     the 6th is to move into the non-local split or      -
c     close                                               -
c                                                         -
c     if the interstitial is split:                       -
c                                                         -
c     the 1st, 2nd and 3rd energies are the barriers to   -
c     the 1, 3 and 5 grafted                              -
c     positions:                                          -
c     same lattice position if alpha, nearest neighbour   -
c     if beta, towards the closest occupied interstitial  -
c     plane                                               -
c                                                         -
c     the 4th, 5th and 6th energies are the barriers to   -
c     the 1, 3 and 5 grafted                              -
c     positions:                                          -
c     same lattice position if alpha, nearest neighbour   -
c     if beta, away from the closest occupied interstitial-
c     plane                                               -
c----------------------------------------------------------
      character*1 cmat(41,20)
      integer*2 smat(21,30),fmat(21,30)
      integer*2 istab(21,30,12,499)
      real*8 ietab(6,499),imtab(6,499)
      real*8 x1(999),y1(999),x2(999),y2(999)
      real*8 x3(999),y3(999),x4(999),y4(999)
      real*8 x5(999),y5(999),x6(999),y6(999)
      real*8 xf1(999),yf1(999),xf2(999),yf2(999)
      real*8 xf3(999),yf3(999),xf4(999),yf4(999)
      real*8 xf5(999),yf5(999),xf6(999),yf6(999)
      integer*2 a1(999),f1(999)
      logical sw1,sw2,emin

      intab = 0

      open(31,FILE='iens.dat')

c-----initialise the library tables

      do k=1,499
         do i=1,21
            do j=1,30
               do n=1,12
                  istab(i,j,n,k) = 3
               end do
            end do
         end do
         ietab(1,k) = 0.0
         ietab(2,k) = 0.0
         ietab(3,k) = 0.0
         ietab(4,k) = 0.0
         ietab(5,k) = 0.0
         ietab(6,k) = 0.0
         imtab(1,k) = 0.0
         imtab(2,k) = 0.0
         imtab(3,k) = 0.0
         imtab(4,k) = 0.0
         imtab(5,k) = 0.0
         imtab(6,k) = 0.0
      end do

      write(*,*)'reading iens.dat file'

      do nrec=1,499

         do i=1,21
            do j=1,30
               smat(i,j) = 0
            end do
         end do

c--------read in structure
         
         do j=1,20 
            read(31,'(41a1)',end=100) ( cmat(i,j), i=1,41 ) 
         end do

c         do j=1,16
c            write(*,'(33a1)') ( cmat(i,j), i=1,33 )
c         end do
         
         ii = 1
         sw1 = .true.
         jj = 1
         ji = 2

c--------convert to internal format

         do j=20,1,-1
            do i=1,41
               if(sw1) then 
                  
                  if(cmat(i,j).eq.' ') then
                     smat(ii,jj) = 0
                  else if(cmat(i,j).eq.'o' .or. cmat(i,j).eq.'O') then
                     smat(ii,jj) = 1
                  else if(cmat(i,j).eq.'-') then
                     smat(ii,jj) = 2
                  else if(cmat(i,j).eq.'*') then
                     smat(ii,jj) = 3
                  else if(cmat(i,j).eq.'+') then
                     smat(ii,jj) = 4
                  else if(cmat(i,j).eq.'1') then
                     smat(ii,jj) = 11
                  else if(cmat(i,j).eq.'2') then
                     smat(ii,jj) = 12
                  else if(cmat(i,j).eq.'3') then
                     smat(ii,jj) = 13
                  else if(cmat(i,j).eq.'4') then
                     smat(ii,jj) = 14
                  else if(cmat(i,j).eq.'5') then
                     smat(ii,jj) = 15
                  else if(cmat(i,j).eq.'6') then
                     smat(ii,jj) = 16
                  else if(cmat(i,j).eq.'a') then
                     smat(ii,jj) = 21
                  else if(cmat(i,j).eq.'b') then
                     smat(ii,jj) = 22
                  else if(cmat(i,j).eq.'c') then
                     smat(ii,jj) = 23
                  else if(cmat(i,j).eq.'d') then
                     smat(ii,jj) = 24
                  else if(cmat(i,j).eq.'e') then
                     smat(ii,jj) = 25
                  else if(cmat(i,j).eq.'f') then
                     smat(ii,jj) = 26
                  else if(cmat(i,j).eq.'A') then
                     smat(ii,jj) = 31
                  else if(cmat(i,j).eq.'B') then
                     smat(ii,jj) = 32
                  else if(cmat(i,j).eq.'C') then
                     smat(ii,jj) = 33
                  else if(cmat(i,j).eq.'D') then
                     smat(ii,jj) = 34
                  else if(cmat(i,j).eq.'E') then
                     smat(ii,jj) = 35
                  else if(cmat(i,j).eq.'F') then
                     smat(ii,jj) = 36
                  else if(cmat(i,j).eq.'s') then
                     smat(ii,jj) = 9
                  else if(cmat(i,j).eq.'#') then
                     smat(ii,jj) = 17
                  else if(cmat(i,j).eq.'@') then
                     smat(ii,jj) = 18
                  else
                     write(*,*)'Non-valid character in iens.dat'
                     write(*,*)cmat(i,j),i,j
                     stop
                  end if
                  
                  ii = ii + 1
                  sw1 = .false.
                  
               else
                  sw1 = .true.
               end if
               
            end do
            
            ii = 1
            sw1 = .true.
            
            if(ji.le.1) then
               jj = jj + 1
               ji = ji + 1
            else
               jj = jj + 2
               ji = 1
            end if
            
         end do



c-----store in the structure library table

         do j=1,30
            do i=1,21
               fmat(i,j) = smat(i,j)
            end do
         end do         

         call flip(fmat)

         do j=1,30
            do i=1,21
               if(fmat(i,j).eq.11) then
                  fmat(i,j) = 16
               else if(fmat(i,j).eq.12) then 
                  fmat(i,j) = 15
               else if(fmat(i,j).eq.13) then
                  fmat(i,j) = 14
               else if(fmat(i,j).eq.14) then
                  fmat(i,j) = 13
               else if(fmat(i,j).eq.15) then 
                  fmat(i,j) = 12
               else if(fmat(i,j).eq.16) then
                  fmat(i,j) = 11
               else if(fmat(i,j).eq.21) then 
                  fmat(i,j) = 26
               else if(fmat(i,j).eq.22) then 
                  fmat(i,j) = 25
               else if(fmat(i,j).eq.23) then 
                  fmat(i,j) = 24
               else if(fmat(i,j).eq.24) then
                  fmat(i,j) = 23
               else if(fmat(i,j).eq.25) then
                  fmat(i,j) = 22
               else if(fmat(i,j).eq.26) then
                  fmat(i,j) = 21
               else if(fmat(i,j).eq.31) then 
                  fmat(i,j) = 36
               else if(fmat(i,j).eq.32) then 
                  fmat(i,j) = 35
               else if(fmat(i,j).eq.33) then 
                  fmat(i,j) = 34
               else if(fmat(i,j).eq.34) then
                  fmat(i,j) = 33
               else if(fmat(i,j).eq.35) then
                  fmat(i,j) = 32
               else if(fmat(i,j).eq.36) then
                  fmat(i,j) = 31
               end if
            end do
         end do      

         natom = 0

         do i=1,21
            do j=1,30
               if(smat(i,j).ge.1) then
                  natom = natom + 1
                  x1(natom) = i*1.23 - 13.5300
                  y1(natom) = j*0.71 - 11.36
                  a1(natom) = smat(i,j)
               end if
            end do
         end do

         nfatom = 0

         do i=1,21
            do j=1,30
               if(fmat(i,j).ge.1) then
                  nfatom = nfatom + 1
                  xf1(nfatom) = i*1.23 - 13.5300
                  yf1(nfatom) = j*0.71 - 11.36
                  f1(nfatom) = fmat(i,j)
               end if
            end do
         end do

         call round(smat)
         call round(fmat)

         do j=1,30
            do i=1,21
               istab(i,j,1,nrec) = smat(i,j)
            end do
         end do


         do j=1,30
            do i=1,21
               istab(i,j,2,nrec) = fmat(i,j)
            end do
         end do


c........rotate 60 degrees

         do i=1,natom
            x2(i)  = 0.5*x1(i) + 0.866025*y1(i)
            y2(i) = - 0.866025*x1(i) + 0.5*y1(i)
            xf2(i)  = 0.5*xf1(i) + 0.866025*yf1(i)
            yf2(i) = - 0.866025*xf1(i) + 0.5*yf1(i)
         end do

         do i=1,21
            do j=1,30
               smat(i,j) = 0
               fmat(i,j) = 0
            end do
         end do
         
c........read back into smat array

         do i=1,21
            do j=1,30
               xt = i*1.23 - 13.5300
               yt = j*0.71 - 11.3600
               
               do n=1,natom
                  if(x2(n).gt.(xt-0.2) .and. x2(n).lt.(xt+0.2) .and.
     &                 y2(n).gt.(yt-0.2) .and. y2(n).lt.(yt+0.2))then
                     smat(i,j) = a1(n)
                  end if
               end do

               do n=1,natom
                  if(xf2(n).gt.(xt-0.2) .and. xf2(n).lt.(xt+0.2) .and.
     &                 yf2(n).gt.(yt-0.2) .and. yf2(n).lt.(yt+0.2))then
                     fmat(i,j) = f1(n)
                  end if
               end do
            end do
         end do         
         
c-----store in the structure library table

         call round(smat)
         call round(fmat)

         do j=1,30
            do i=1,21
               if(fmat(i,j).ge.11 .and. fmat(i,j).le.16) then
                  fmat(i,j) = fmat(i,j) + 1
                  if(fmat(i,j).gt.16) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.21 .and. fmat(i,j).le.26) then
                  fmat(i,j) = fmat(i,j) + 1
                  if(fmat(i,j).gt.26) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.31 .and. fmat(i,j).le.36) then
                  fmat(i,j) = fmat(i,j) + 1
                  if(fmat(i,j).gt.36) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if

               if(smat(i,j).ge.11 .and. smat(i,j).le.16) then
                  smat(i,j) = smat(i,j) + 1
                  if(smat(i,j).gt.16) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.21 .and. smat(i,j).le.26) then
                  smat(i,j) = smat(i,j) + 1
                  if(smat(i,j).gt.26) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.31 .and. smat(i,j).le.36) then
                  smat(i,j) = smat(i,j) + 1
                  if(smat(i,j).gt.36) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
            end do
         end do      

         do j=1,30
            do i=1,21
               istab(i,j,3,nrec) = smat(i,j)
            end do
         end do

c........flip and store

         do j=1,30
            do i=1,21
               istab(i,j,4,nrec) = fmat(i,j)
            end do
         end do
        
c........rotate 60 degrees

         do i=1,natom
            x3(i)  = 0.5*x2(i) + 0.866025*y2(i)
            y3(i) = - 0.866025*x2(i) + 0.5*y2(i)
            xf3(i)  = 0.5*xf2(i) + 0.866025*yf2(i)
            yf3(i) = - 0.866025*xf2(i) + 0.5*yf2(i)
         end do

         do i=1,21
            do j=1,30
               smat(i,j) = 0
               fmat(i,j) = 0
            end do
         end do
         
c........read back into smat array

         do i=1,21
            do j=1,30
               xt = i*1.23 - 13.5300
               yt = j*0.71 - 11.3600
               
               do n=1,natom
                  if(x3(n).gt.(xt-0.2) .and. x3(n).lt.(xt+0.2) .and.
     &                 y3(n).gt.(yt-0.2) .and. y3(n).lt.(yt+0.2))then
                     smat(i,j) = a1(n)
                  end if
               end do
               do n=1,natom
                  if(xf3(n).gt.(xt-0.2) .and. xf3(n).lt.(xt+0.2) .and.
     &                 yf3(n).gt.(yt-0.2) .and. yf3(n).lt.(yt+0.2))then
                     fmat(i,j) = f1(n)
                  end if
               end do
            end do
         end do         
         

c-----store in the structure library table

         call round(smat)
         call round(fmat)

         do j=1,30
            do i=1,21
               if(fmat(i,j).ge.11 .and. fmat(i,j).le.16) then
                  fmat(i,j) = fmat(i,j) + 2
                  if(fmat(i,j).gt.16) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.21 .and. fmat(i,j).le.26) then
                  fmat(i,j) = fmat(i,j) + 2
                  if(fmat(i,j).gt.26) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.31 .and. fmat(i,j).le.36) then
                  fmat(i,j) = fmat(i,j) + 2
                  if(fmat(i,j).gt.36) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if


               if(smat(i,j).ge.11 .and. smat(i,j).le.16) then
                  smat(i,j) = smat(i,j) + 2
                  if(smat(i,j).gt.16) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.21 .and. smat(i,j).le.26) then
                  smat(i,j) = smat(i,j) + 2
                  if(smat(i,j).gt.26) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.31 .and. smat(i,j).le.36) then
                  smat(i,j) = smat(i,j) + 2
                  if(smat(i,j).gt.36) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
            end do
         end do      

         do j=1,30
            do i=1,21
               istab(i,j,5,nrec) = smat(i,j)
            end do
         end do

c........flip and store

         do j=1,30
            do i=1,21
               istab(i,j,6,nrec) = fmat(i,j)
            end do
         end do

c........rotate 60 degrees

         do i=1,natom
            x4(i)  = 0.5*x3(i) + 0.866025*y3(i)
            y4(i) = - 0.866025*x3(i) + 0.5*y3(i)
            xf4(i)  = 0.5*xf3(i) + 0.866025*yf3(i)
            yf4(i) = - 0.866025*xf3(i) + 0.5*yf3(i)
         end do

         do i=1,21
            do j=1,30
               smat(i,j) = 0
               fmat(i,j) = 0
            end do
         end do
         
c........read back into smat array

         do i=1,21
            do j=1,30
               xt = i*1.23 - 13.5300
               yt = j*0.71 - 11.3600
               
               do n=1,natom
                  if(x4(n).gt.(xt-0.2) .and. x4(n).lt.(xt+0.2) .and.
     &                 y4(n).gt.(yt-0.2) .and. y4(n).lt.(yt+0.2))then
                     smat(i,j) = a1(n)
                  end if
               end do
               do n=1,natom
                  if(xf4(n).gt.(xt-0.2) .and. xf4(n).lt.(xt+0.2) .and.
     &                 yf4(n).gt.(yt-0.2) .and. yf4(n).lt.(yt+0.2))then
                     fmat(i,j) = f1(n)
                  end if
               end do
            end do
         end do         
         
c-----store in the structure library table

         call round(smat)
         call round(fmat)

         do j=1,30
            do i=1,21
               if(fmat(i,j).ge.11 .and. fmat(i,j).le.16) then
                  fmat(i,j) = fmat(i,j) + 3
                  if(fmat(i,j).gt.16) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.21 .and. fmat(i,j).le.26) then
                  fmat(i,j) = fmat(i,j) + 3
                  if(fmat(i,j).gt.26) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.31 .and. fmat(i,j).le.36) then
                  fmat(i,j) = fmat(i,j) + 3
                  if(fmat(i,j).gt.36) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if

               if(smat(i,j).ge.11 .and. smat(i,j).le.16) then
                  smat(i,j) = smat(i,j) + 3
                  if(smat(i,j).gt.16) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.21 .and. smat(i,j).le.26) then
                  smat(i,j) = smat(i,j) + 3
                  if(smat(i,j).gt.26) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.31 .and. smat(i,j).le.36) then
                  smat(i,j) = smat(i,j) + 3
                  if(smat(i,j).gt.36) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
            end do
         end do      

         do j=1,30
            do i=1,21
               istab(i,j,7,nrec) = smat(i,j)
            end do
         end do

         do j=1,30
            do i=1,21
               istab(i,j,8,nrec) = fmat(i,j)
            end do
         end do

c........rotate 60 degrees

         do i=1,natom
            x5(i)  = 0.5*x4(i) + 0.866025*y4(i)
            y5(i) = - 0.866025*x4(i) + 0.5*y4(i)
            xf5(i)  = 0.5*xf4(i) + 0.866025*yf4(i)
            yf5(i) = - 0.866025*xf4(i) + 0.5*yf4(i)
         end do

         do i=1,21
            do j=1,30
               smat(i,j) = 0
               fmat(i,j) = 0
            end do
         end do
         
c........read back into smat array

         do i=1,21
            do j=1,30
               xt = i*1.23 - 13.5300
               yt = j*0.71 - 11.3600
               
               do n=1,natom
                  if(x5(n).gt.(xt-0.2) .and. x5(n).lt.(xt+0.2) .and.
     &                 y5(n).gt.(yt-0.2) .and. y5(n).lt.(yt+0.2))then
                     smat(i,j) = a1(n)
                  end if
               end do
               do n=1,natom
                  if(xf5(n).gt.(xt-0.2) .and. xf5(n).lt.(xt+0.2) .and.
     &                 yf5(n).gt.(yt-0.2) .and. yf5(n).lt.(yt+0.2))then
                     fmat(i,j) = f1(n)
                  end if
               end do
            end do
         end do         
         
c-----store in the structure library table

         call round(smat)
         call round(fmat)

         do j=1,30
            do i=1,21
               if(fmat(i,j).ge.11 .and. fmat(i,j).le.16) then
                  fmat(i,j) = fmat(i,j) + 4
                  if(fmat(i,j).gt.16) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.21 .and. fmat(i,j).le.26) then
                  fmat(i,j) = fmat(i,j) + 4
                  if(fmat(i,j).gt.26) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.31 .and. fmat(i,j).le.36) then
                  fmat(i,j) = fmat(i,j) + 4
                  if(fmat(i,j).gt.36) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if

               if(smat(i,j).ge.11 .and. smat(i,j).le.16) then
                  smat(i,j) = smat(i,j) + 4
                  if(smat(i,j).gt.16) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.21 .and. smat(i,j).le.26) then
                  smat(i,j) = smat(i,j) + 4
                  if(smat(i,j).gt.26) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.31 .and. smat(i,j).le.36) then
                  smat(i,j) = smat(i,j) + 4
                  if(smat(i,j).gt.36) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
            end do
         end do      

         do j=1,30
            do i=1,21
               istab(i,j,9,nrec) = smat(i,j)
            end do
         end do

         do j=1,30
            do i=1,21
               istab(i,j,10,nrec) = fmat(i,j)
            end do
         end do

c........rotate 60 degrees

         do i=1,natom
            x6(i)  = 0.5*x5(i) + 0.866025*y5(i)
            y6(i) = - 0.866025*x5(i) + 0.5*y5(i)
            xf6(i)  = 0.5*xf5(i) + 0.866025*yf5(i)
            yf6(i) = - 0.866025*xf5(i) + 0.5*yf5(i)
         end do

         do i=1,21
            do j=1,30
               smat(i,j) = 0
               fmat(i,j) = 0
            end do
         end do
         
c........read back into smat array

         do i=1,21
            do j=1,30
               xt = i*1.23 - 13.5300
               yt = j*0.71 - 11.3600
               
               do n=1,natom
                  if(x6(n).gt.(xt-0.2) .and. x6(n).lt.(xt+0.2) .and.
     &                 y6(n).gt.(yt-0.2) .and. y6(n).lt.(yt+0.2))then
                     smat(i,j) = a1(n)
                  end if
               end do
               do n=1,natom
                  if(xf6(n).gt.(xt-0.2) .and. xf6(n).lt.(xt+0.2) .and.
     &                 yf6(n).gt.(yt-0.2) .and. yf6(n).lt.(yt+0.2))then
                     fmat(i,j) = f1(n)
                  end if
               end do
            end do
         end do         
         
c-----store in the structure library table

         call round(smat)
         call round(fmat)

         do j=1,30
            do i=1,21
               if(fmat(i,j).ge.11 .and. fmat(i,j).le.16) then
                  fmat(i,j) = fmat(i,j) + 5
                  if(fmat(i,j).gt.16) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.21 .and. fmat(i,j).le.26) then
                  fmat(i,j) = fmat(i,j) + 5
                  if(fmat(i,j).gt.26) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if
               if(fmat(i,j).ge.31 .and. fmat(i,j).le.36) then
                  fmat(i,j) = fmat(i,j) + 5
                  if(fmat(i,j).gt.36) then
                     fmat(i,j) = fmat(i,j) - 6
                  end if
               end if


               if(smat(i,j).ge.11 .and. smat(i,j).le.16) then
                  smat(i,j) = smat(i,j) + 5
                  if(smat(i,j).gt.26) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.21 .and. smat(i,j).le.26) then
                  smat(i,j) = smat(i,j) + 5
                  if(smat(i,j).gt.26) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
               if(smat(i,j).ge.31 .and. smat(i,j).le.36) then
                  smat(i,j) = smat(i,j) + 5
                  if(smat(i,j).gt.36) then
                     smat(i,j) = smat(i,j) - 6
                  end if
               end if
            end do
         end do      

         do j=1,30
            do i=1,21
               istab(i,j,11,nrec) = smat(i,j)
            end do
         end do

         do j=1,30
            do i=1,21
               istab(i,j,12,nrec) = fmat(i,j)
            end do
         end do

c-----read in energies

         read(31,*,end=100)ietab(1,nrec),ietab(2,nrec),ietab(3,nrec)
         read(31,*,end=100)ietab(4,nrec),ietab(5,nrec),ietab(6,nrec)

         if(emin) then
            read(31,*,end=100)imtab(1,nrec),imtab(2,nrec),imtab(3,nrec)
            read(31,*,end=100)imtab(4,nrec),imtab(5,nrec),imtab(6,nrec)
         end if
            
      end do

      write(*,*)'Too many structures - increase tab arrays'
      stop

 100  intab = nrec - 1

      write(*,*)'found ',intab,' structures'

c......writing coords

c      do j=24,1,-1
c         write(*,'(17i3)') ( istab(i,j,10,1) , i=1,17 )
c      end do

c      do j=1,24
c         do i=1,17
c            smat(i,j) = stab(i,j,12,1)
c            if(abs((i-9)*(j-13)).ge.50) then
c               smat(i,j) = 3
c            end if
c            
c         end do
c      end do

c      do j=24,1,-1
c         write(*,'(17i1)') ( smat(i,j) , i=1,17 )
c      end do

c      do i=1,17
c         do j=1,24
c            do n=1,11,2
c               write(*,*)n
c            if(stab(i,j,n,1).eq.1) then
c               write(32,'(a2,3(x,f8.4))')'C ',i*1.23 -  11.0700,
c     &              j*0.71 - 9.2300,0.0
c            else if(stab(i,j,n,1).eq.2) then
c               write(32,'(a2,3(x,f8.4))')'B ',i*1.23 -  11.0700,
c     &              j*0.71 - 9.2300,0.0
c            else if(stab(i,j,n,1).eq.3) then
c               write(32,'(a2,3(x,f8.4))')'Si',i*1.23 -  11.0700,
c     &              j*0.71 - 9.2300,0.0
c            else if(stab(i,j,n,1).eq.4) then
c               write(32,'(a2,3(x,f8.4))')'O ',i*1.23 -  11.0700,
c     &              j*0.71 - 9.2300,0.0
c            end if
c            end do
c         end do
c      end do
      
      end 
