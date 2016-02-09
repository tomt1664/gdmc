      subroutine strans(system,isystem,move)
c---------------------------------------------------------
c     move atoms acording to the site and kmc direction  - 
c     move(1-3) : site indices                           -
c     move(4) : direction                                -
c     11: split interstitial to move up | or /           -
c     12: split interstitial to move up \ or |           -
c     13: split interstitial to move up / or \           -
c     14: split interstitial to move down | or /         -
c     15: split interstitial to move down \ or |         -
c     16: split interstitial to move down / or \         -
c---------------------------------------------------------
      integer*1 system(800,800,199),isystem(800,800,199)
      integer move(4)
      integer*2 atm1,atm2,atm3

      i = move(1)
      j = move(2)
      k = move(3)

      if(system(i,j,k).ne.9) then
         write(*,*)'error in strans - not split site'
         stop
      end if

      system(i,j,k) = 1

c*****move to upper layer

      if(move(4).eq.11) then
c........if up
         if(system(i,md(j+2),k).ne.0) then
c-----split at alpha site goes to grafted at same site
            if(isystem(i,j,k).eq.8) then
               isystem(i,j,k) = 11
            else if(isystem(i,j,k).eq.0) then
c-----split at beta site goes to grafted at neighbouring site
               isystem(i,md(j+2),k) = 14
            else
               write(*,*)'error in strans 11 up'
               stop
            end if
c........if down
         else
c-----split at alpha site goes to grafted at same site
            if(isystem(i,j,k).eq.8) then
               isystem(i,j,k) = 12
            else if(isystem(i,j,k).eq.0) then
c-----split at beta site goes to grafted at neighbouring site
               isystem(mb(i+1),md(j+1),k) = 15
            else
               write(*,*)'error in strans 11 down'
               stop
            end if 
         end if


      else if(move(4).eq.12) then
c........if up
         if(system(i,md(j+2),k).ne.0) then
c-----split at alpha site goes to grafted at same site
            if(isystem(i,j,k).eq.8) then
               isystem(i,j,k) = 13
            else if(isystem(i,j,k).eq.0) then
c-----split at beta site goes to grafted at neighbouring site
               isystem(mb(i+1),md(j-1),k) = 16
            else
               write(*,*)'error in strans 12 up'
               stop
            end if
c........if down
         else
c-----split at alpha site goes to grafted at same site
            if(isystem(i,j,k).eq.8) then
               isystem(i,j,k) = 14
            else if(isystem(i,j,k).eq.0) then
c-----split at beta site goes to grafted at neighbouring site
               isystem(i,md(j-2),k) = 11
            else
               write(*,*)'error in strans 12 down'
               stop
            end if 
         end if


      else if(move(4).eq.13) then
c........if up
         if(system(i,md(j+2),k).ne.0) then
c-----split at alpha site goes to grafted at same site
            if(isystem(i,j,k).eq.8) then
               isystem(i,j,k) = 15
            else if(isystem(i,j,k).eq.0) then
c-----split at beta site goes to grafted at neighbouring site
               isystem(mb(i-1),md(j-1),k) = 12
            else
               write(*,*)'error in strans 13 up'
               stop
            end if
c........if down
         else
c-----split at alpha site goes to grafted at same site
            if(isystem(i,j,k).eq.8) then
               isystem(i,j,k) = 16
            else if(isystem(i,j,k).eq.0) then
c-----split at beta site goes to grafted at neighbouring site
               isystem(mb(i-1),md(j+1),k) = 13
            else
               write(*,*)'error in strans 13 down'
               stop
            end if 
         end if

c*******move to lower layer

      else if(move(4).eq.14) then
c........if up
         if(system(i,md(j+2),k).ne.0) then
c-----split at alpha site goes to grafted at same site
            if(isystem(i,j,ml(k-1)).eq.8) then
               isystem(i,j,ml(k-1)) = 21
            else if(isystem(i,j,ml(k-1)).eq.0) then
c-----split at beta site goes to grafted at neighbouring site
               isystem(i,md(j+2),ml(k-1)) = 24
            else
               write(*,*)'error in strans 14 up'
               stop
            end if

c........if down
         else
c-----split at alpha site goes to grafted at same site
            if(isystem(i,j,ml(k-1)).eq.8) then
               isystem(i,j,ml(k-1)) = 22
            else if(isystem(i,j,ml(k-1)).eq.0) then
c-----split at beta site goes to grafted at neighbouring site
               isystem(mb(i+1),md(j+1),ml(k-1)) = 25
            else
               write(*,*)'error in strans 14 down'
               stop
            end if 
         end if


      else if(move(4).eq.15) then
c........if up
         if(system(i,md(j+2),k).ne.0) then
c-----split at alpha site goes to grafted at same site
            if(isystem(i,j,ml(k-1)).eq.8) then
               isystem(i,j,ml(k-1)) = 23
            else if(isystem(i,j,ml(k-1)).eq.0) then
c-----split at beta site goes to grafted at neighbouring site
               isystem(mb(i+1),md(j-1),ml(k-1)) = 26
            else
               write(*,*)'error in strans 15 up'
               stop
            end if
c........if down
         else
c-----split at alpha site goes to grafted at same site
            if(isystem(i,j,ml(k-1)).eq.8) then
               isystem(i,j,ml(k-1)) = 24
            else if(isystem(i,j,ml(k-1)).eq.0) then
c-----split at beta site goes to grafted at neighbouring site
               isystem(i,md(j-2),ml(k-1)) = 21
            else
               write(*,*)'error in strans 15 down'
               stop
            end if 
         end if


      else if(move(4).eq.16) then
c........if up
         if(system(i,md(j+2),k).ne.0) then
c-----split at alpha site goes to grafted at same site
            if(isystem(i,j,ml(k-1)).eq.8) then
               isystem(i,j,ml(k-1)) = 25
            else if(isystem(i,j,ml(k-1)).eq.0) then
c-----split at beta site goes to grafted at neighbouring site
               isystem(mb(i-1),md(j-1),ml(k-1)) = 22
            else
               write(*,*)'error in strans 16 up'
               stop
            end if
c........if down
         else
c-----split at alpha site goes to grafted at same site
            if(isystem(i,j,ml(k-1)).eq.8) then
               isystem(i,j,ml(k-1)) = 26
            else if(isystem(i,j,ml(k-1)).eq.0) then
c-----split at beta site goes to grafted at neighbouring site
               isystem(mb(i-1),md(j+1),ml(k-1)) = 23
            else
               write(*,*)'error in strans 16 down'
               stop
            end if 
         end if
         
      end if
      
      end


      subroutine itrans(system,isystem,move)
c---------------------------------------------------------
c     move atoms acording to the site and kmc direction  - 
c     move(1-3) : site indices                           -
c     move(4) : direction                                -
c     21: interstitial to move clockwise                 -
c     22: interstitial to move anti-clockwise            -
c     23: interstitial to move to neighbouring alpha     -
c     24: grafted to move to spiro                       -
c         sprio to move up grafted                       -
c     25: grafted to move to local alpha split           -
c         spiro to move down grafted                     -
c     26: grafted to move to non-local beta split        -
c         spiro to close intimate pair                   -
c---------------------------------------------------------
      integer*1 system(800,800,199),isystem(800,800,199)
      integer move(4)
      logical up

      i = move(1)
      j = move(2)
      k = move(3)

c*****spiro

      if(isystem(i,j,k).ge.1 .and.isystem(i,j,k).le.6) then
         
         if(move(4).eq.21) then
            isystem(i,j,k) = isystem(i,j,k) + 1
            if(isystem(i,j,k).gt.6) then
               isystem(i,j,k) = isystem(i,j,k) - 6
            end if
         else if(move(4).eq.22) then
            isystem(i,j,k) = isystem(i,j,k) - 1
            if(isystem(i,j,k).lt.1) then
               isystem(i,j,k) = isystem(i,j,k) + 6
            end if
         else if(move(4).eq.23) then
            if(isystem(i,j,k).eq.1) then
               if(isystem(mb(i+1),md(j+3),k).eq.8) then
                  isystem(mb(i+1),md(j+3),k) = 4
               else
                  write(*,*)'error in itrans spiro 23 1'
                  stop
               end if
            else if(isystem(i,j,k).eq.2) then
               if(isystem(mb(i+2),j,k).eq.8) then
                  isystem(mb(i+2),j,k) = 5
               else
                  write(*,*)'error in itrans spiro 23 2'
                  stop
               end if
            else if(isystem(i,j,k).eq.3) then
               if(isystem(mb(i+1),md(j-3),k).eq.8) then
                  isystem(mb(i+1),md(j-3),k) = 6
               else
                  write(*,*)'error in itrans spiro 23 3'
                  stop
               end if
            elseif(isystem(i,j,k).eq.4) then
               if(isystem(mb(i-1),md(j-3),k).eq.8) then
                  isystem(mb(i-1),md(j-3),k) = 1
               else
                  write(*,*)'error in itrans spiro 23 4'
                  stop
               end if
            else if(isystem(i,j,k).eq.5) then
               if(isystem(mb(i-2),j,k).eq.8) then
                  isystem(mb(i-2),j,k) = 2
               else
                  write(*,*)'error in itrans spiro 23 5'
                  stop
               end if
            else if(isystem(i,j,k).eq.6) then
               if(isystem(mb(i-1),md(j+3),k).eq.8) then
                  isystem(mb(i-1),md(j+3),k) = 3
               else
                  write(*,*)'error in itrans spiro 23 6'
                  stop
               end if
            end if
            isystem(i,j,k) = 8
         else if(move(4).eq.24) then
            if(system(i,md(j+2),ml(k+1)).ne.0) then
               if(isystem(i,j,k).eq.1 .or. isystem(i,j,k).eq.6) then
                  isystem(i,j,k) = 21
               else if(isystem(i,j,k).eq.2 .or. isystem(i,j,k).eq.3)then
                  isystem(i,j,k) = 23
               else if(isystem(i,j,k).eq.4 .or. isystem(i,j,k).eq.5)then
                  isystem(i,j,k) = 25
               end if
            else
               if(isystem(i,j,k).eq.1 .or. isystem(i,j,k).eq.2) then
                  isystem(i,j,k) = 22
               else if(isystem(i,j,k).eq.3 .or. isystem(i,j,k).eq.4)then
                  isystem(i,j,k) = 24
               else if(isystem(i,j,k).eq.5 .or. isystem(i,j,k).eq.6)then
                  isystem(i,j,k) = 26
               end if 
            end if
         else if(move(4).eq.25) then
            if(system(i,md(j+2),k).ne.0) then
               if(isystem(i,j,k).eq.1 .or. isystem(i,j,k).eq.6) then
                  isystem(i,j,k) = 11
               else if(isystem(i,j,k).eq.2 .or. isystem(i,j,k).eq.3)then
                  isystem(i,j,k) = 13
               else if(isystem(i,j,k).eq.4 .or. isystem(i,j,k).eq.5)then
                  isystem(i,j,k) = 15
               end if
            else
               if(isystem(i,j,k).eq.1 .or. isystem(i,j,k).eq.2) then
                  isystem(i,j,k) = 12
               else if(isystem(i,j,k).eq.3 .or. isystem(i,j,k).eq.4)then
                  isystem(i,j,k) = 14
               else if(isystem(i,j,k).eq.5 .or. isystem(i,j,k).eq.6)then
                  isystem(i,j,k) = 16
               end if 
            end if




         else if(move(4).eq.26) then
c-----------close frenkel close
            if(system(i,j,k).eq.2) then
               system(i,j,k) = 1
c - - - - - - -close interlayer divac

               if(system(i,md(j+2),k).eq.5) then
                  system(i,md(j+2),k) = 1
                  if(system(i,md(j+2),ml(k+1)).eq.6) then
                     system(i,md(j+2),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j+1),k).eq.5) then
                  system(mb(i+1),md(j+1),k) = 1
                  if(system(mb(i+1),md(j+1),ml(k+1)).eq.6) then
                     system(mb(i+1),md(j+1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1),md(j-1),k).eq.5) then
                  system(mb(i+1),md(j-1),k) = 1
                  if(system(mb(i+1),md(j-1),ml(k+1)).eq.6) then
                     system(mb(i+1),md(j-1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(i,md(j-2),k).eq.5) then
                  system(i,md(j-2),k) = 1
                  if(system(i,md(j-2),ml(k+1)).eq.6) then
                     system(i,md(j-2),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1),md(j-1),k).eq.5) then
                  system(mb(i-1),md(j-1),k) = 1
                  if(system(mb(i-1),md(j-1),ml(k+1)).eq.6) then
                     system(mb(i-1),md(j-1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j+1),k).eq.5) then
                  system(mb(i-1),md(j+1),k) = 1
                  if(system(mb(i-1),md(j+1),ml(k+1)).eq.6) then
                     system(mb(i-1),md(j+1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 501
               end if

               if(up) then
                  if(system(i,md(j+4),ml(k+1)).eq.6) then
                     system(i,md(j+4),ml(k+1)) = 1
                  else if(system(mb(i+2),md(j-2),ml(k+1)).eq.6) then
                     system(mb(i+2),md(j-2),ml(k+1)) = 1
                  else if(system(mb(i-2),md(j-2),ml(k+1)).eq.6) then
                     system(mb(i-2),md(j-2),ml(k+1)) = 1
                  else
                     write(*,*)'error in move 26 - 1 - 1'
                  end if
                  go to 333
               else
                  if(system(i,md(j-4),ml(k+1)).eq.6) then
                     system(i,md(j-4),ml(k+1)) = 1
                  else if(system(mb(i-2),md(j+2),ml(k+1)).eq.6) then
                     system(mb(i-2),md(j+2),ml(k+1)) = 1
                  else if(system(mb(i+2),md(j+2),ml(k+1)).eq.6) then
                     system(mb(i+2),md(j+2),ml(k+1)) = 1
                  else
                     write(*,*)'error in move 26 - 1 - 2'
                  end if
                  go to 333
               end if


 501           if(system(i,md(j+2),k).eq.6) then
                  system(i,md(j+2),k) = 1
                  if(system(i,md(j+2),ml(k-1)).eq.5) then
                     system(i,md(j+2),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j+1),k).eq.6) then
                  system(mb(i+1),md(j+1),k) = 1
                  if(system(mb(i+1),md(j+1),ml(k-1)).eq.5) then
                     system(mb(i+1),md(j+1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1),md(j-1),k).eq.6) then
                  system(mb(i+1),md(j-1),k) = 1
                  if(system(mb(i+1),md(j-1),ml(k-1)).eq.5) then
                     system(mb(i+1),md(j-1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(i,md(j-2),k).eq.6) then
                  system(i,md(j-2),k) = 1
                  if(system(i,md(j-2),ml(k-1)).eq.5) then
                     system(i,md(j-2),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1),md(j-1),k).eq.6) then
                  system(mb(i-1),md(j-1),k) = 1
                  if(system(mb(i-1),md(j-1),ml(k-1)).eq.5) then
                     system(mb(i-1),md(j-1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j+1),k).eq.6) then
                  system(mb(i-1),md(j+1),k) = 1
                  if(system(mb(i-1),md(j+1),ml(k-1)).eq.5) then
                     system(mb(i-1),md(j+1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 333
               end if

               if(up) then
                  if(system(i,md(j+4),ml(k-1)).eq.5) then
                     system(i,md(j+4),ml(k-1)) = 1
                  else if(system(mb(i+2),md(j-2),ml(k-1)).eq.5) then
                     system(mb(i+2),md(j-2),ml(k-1)) = 1
                  else if(system(mb(i-2),md(j-2),ml(k-1)).eq.5) then
                     system(mb(i-2),md(j-2),ml(k-1)) = 1
                  else
                     write(*,*)'error in move 26 - 2 - 1'
                  end if
                  go to 333
               else
                  if(system(i,md(j-4),ml(k-1)).eq.5) then
                     system(i,md(j-4),ml(k-1)) = 1
                  else if(system(mb(i-2),md(j+2),ml(k-1)).eq.5) then
                     system(mb(i-2),md(j+2),ml(k-1)) = 1
                  else if(system(mb(i+2),md(j+2),ml(k-1)).eq.5) then
                     system(mb(i+2),md(j+2),ml(k-1)) = 1
                  else
                     write(*,*)'error in move 26 - 2 - 2'
                  end if
                  go to 333
               end if               

            else if(system(i,j,ml(k+1)).eq.2) then
               system(i,j,ml(k+1)) = 1
c - - - - - - -close interlayer divac

               if(system(i,md(j+2),ml(k+1)).eq.5) then
                  system(i,md(j+2),ml(k+1)) = 1
                  if(system(i,md(j+2),ml(k+2)).eq.6) then
                     system(i,md(j+2),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j+1),ml(k+1)).eq.5) then
                  system(mb(i+1),md(j+1),ml(k+1)) = 1
                  if(system(mb(i+1),md(j+1),ml(k+2)).eq.6) then
                     system(mb(i+1),md(j+1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1),md(j-1),ml(k+1)).eq.5) then
                  system(mb(i+1),md(j-1),ml(k+1)) = 1
                  if(system(mb(i+1),md(j-1),ml(k+2)).eq.6) then
                     system(mb(i+1),md(j-1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(i,md(j-2),ml(k+1)).eq.5) then
                  system(i,md(j-2),ml(k+1)) = 1
                  if(system(i,md(j-2),ml(k+2)).eq.6) then
                     system(i,md(j-2),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1),md(j-1),ml(k+1)).eq.5) then
                  system(mb(i-1),md(j-1),ml(k+1)) = 1
                  if(system(mb(i-1),md(j-1),ml(k+2)).eq.6) then
                     system(mb(i-1),md(j-1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j+1),ml(k+1)).eq.5) then
                  system(mb(i-1),md(j+1),ml(k+1)) = 1
                  if(system(mb(i-1),md(j+1),ml(k+2)).eq.6) then
                     system(mb(i-1),md(j+1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 502
               end if

               if(up) then
                  if(system(i,md(j+4),ml(k+2)).eq.6) then
                     system(i,md(j+4),ml(k+2)) = 1
                  else if(system(mb(i+2),md(j-2),ml(k+2)).eq.6) then
                     system(mb(i+2),md(j-2),ml(k+2)) = 1
                  else if(system(mb(i-2),md(j-2),ml(k+2)).eq.6) then
                     system(mb(i-2),md(j-2),ml(k+2)) = 1
                  else
                     write(*,*)'error in move 26 - 3 - 1'
                  end if
                  go to 333
               else
                  if(system(i,md(j-4),ml(k+2)).eq.6) then
                     system(i,md(j-4),ml(k+2)) = 1
                  else if(system(mb(i-2),md(j+2),ml(k+2)).eq.6) then
                     system(mb(i-2),md(j+2),ml(k+2)) = 1
                  else if(system(mb(i+2),md(j+2),ml(k+2)).eq.6) then
                     system(mb(i+2),md(j+2),ml(k+2)) = 1
                  else
                     write(*,*)'error in move 26 - 3 - 2'
                  end if
                  go to 333
               end if



 502           if(system(i,md(j+2),ml(k+1)).eq.6) then
                  system(i,md(j+2),ml(k+1)) = 1
                  if(system(i,md(j+2),ml(k)).eq.5) then
                     system(i,md(j+2),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j+1),ml(k+1)).eq.6) then
                  system(mb(i+1),md(j+1),ml(k+1)) = 1
                  if(system(mb(i+1),md(j+1),ml(k)).eq.5) then
                     system(mb(i+1),md(j+1),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1),md(j-1),ml(k+1)).eq.6) then
                  system(mb(i+1),md(j-1),ml(k+1)) = 1
                  if(system(mb(i+1),md(j-1),ml(k)).eq.5) then
                     system(mb(i+1),md(j-1),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(i,md(j-2),ml(k+1)).eq.6) then
                  system(i,md(j-2),ml(k+1)) = 1
                  if(system(i,md(j-2),ml(k)).eq.5) then
                     system(i,md(j-2),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1),md(j-1),ml(k+1)).eq.6) then
                  system(mb(i-1),md(j-1),ml(k+1)) = 1
                  if(system(mb(i-1),md(j-1),ml(k)).eq.5) then
                     system(mb(i-1),md(j-1),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j+1),ml(k+1)).eq.6) then
                  system(mb(i-1),md(j+1),ml(k+1)) = 1
                  if(system(mb(i-1),md(j+1),ml(k)).eq.5) then
                     system(mb(i-1),md(j+1),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 333
               end if

               if(up) then
                  if(system(i,md(j+4),ml(k)).eq.5) then
                     system(i,md(j+4),ml(k)) = 1
                  else if(system(mb(i+2),md(j-2),ml(k)).eq.5) then
                     system(mb(i+2),md(j-2),ml(k)) = 1
                  else if(system(mb(i-2),md(j-2),ml(k)).eq.5) then
                     system(mb(i-2),md(j-2),ml(k)) = 1
                  else
                     write(*,*)'error in move 26 - 4 - 1'
                  end if
                  go to 333
               else
                  if(system(i,md(j-4),ml(k)).eq.5) then
                     system(i,md(j-4),ml(k)) = 1
                  else if(system(mb(i-2),md(j+2),ml(k)).eq.5) then
                     system(mb(i-2),md(j+2),ml(k)) = 1
                  else if(system(mb(i+2),md(j+2),ml(k)).eq.5) then
                     system(mb(i+2),md(j+2),ml(k)) = 1
                  else
                     write(*,*)'error in move 26 - 4 - 2'
                  end if
                  go to 333
               end if               


c..............intimate dir 1
            else if(system(i,md(j+2),k).eq.2) then
               system(i,md(j+2),k) = 1

c - - - - - - -close interlayer divac

               if(system(i,md(j+2+2),k).eq.5) then
                  system(i,md(j+2+2),k) = 1
                  if(system(i,md(j+2+2),ml(k+1)).eq.6) then
                     system(i,md(j+2+2),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j+1+2),k).eq.5) then
                  system(mb(i+1),md(j+1+2),k) = 1
                  if(system(mb(i+1),md(j+1+2),ml(k+1)).eq.6) then
                     system(mb(i+1),md(j+1+2),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1),md(j-1+2),k).eq.5) then
                  system(mb(i+1),md(j-1+2),k) = 1
                  if(system(mb(i+1),md(j-1+2),ml(k+1)).eq.6) then
                     system(mb(i+1),md(j-1+2),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(i,md(j-2+2),k).eq.5) then
                  system(i,md(j-2+2),k) = 1
                  if(system(i,md(j-2+2),ml(k+1)).eq.6) then
                     system(i,md(j-2+2),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1),md(j-1+2),k).eq.5) then
                  system(mb(i-1),md(j-1+2),k) = 1
                  if(system(mb(i-1),md(j-1+2),ml(k+1)).eq.6) then
                     system(mb(i-1),md(j-1+2),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j+1+2),k).eq.5) then
                  system(mb(i-1),md(j+1+2),k) = 1
                  if(system(mb(i-1),md(j+1+2),ml(k+1)).eq.6) then
                     system(mb(i-1),md(j+1+2),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 503
               end if

               if(up) then
                  if(system(i,md(j+4+2),ml(k+1)).eq.6) then
                     system(i,md(j+4+2),ml(k+1)) = 1
                  else if(system(mb(i+2),md(j-2+2),ml(k+1)).eq.6) then
                     system(mb(i+2),md(j-2+2),ml(k+1)) = 1
                  else if(system(mb(i-2),md(j-2+2),ml(k+1)).eq.6) then
                     system(mb(i-2),md(j-2+2),ml(k+1)) = 1
                  else
                     write(*,*)'error in move 26 - 5 - 1'
                  end if
                  go to 333
               else
                  if(system(i,md(j-4+2),ml(k+1)).eq.6) then
                     system(i,md(j-4+2),ml(k+1)) = 1
                  else if(system(mb(i-2),md(j+2+2),ml(k+1)).eq.6) then
                     system(mb(i-2),md(j+2+2),ml(k+1)) = 1
                  else if(system(mb(i+2),md(j+2+2),ml(k+1)).eq.6) then
                     system(mb(i+2),md(j+2+2),ml(k+1)) = 1
                  else
                     write(*,*)'error in move 26 - 5 - 2'
                  end if
                  go to 333
               end if



 503           if(system(i,md(j+2+2),k).eq.6) then
                  system(i,md(j+2+2),k) = 1
                  if(system(i,md(j+2+2),ml(k-1)).eq.5) then
                     system(i,md(j+2+2),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j+1+2),k).eq.6) then
                  system(mb(i+1),md(j+1+2),k) = 1
                  if(system(mb(i+1),md(j+1+2),ml(k-1)).eq.5) then
                     system(mb(i+1),md(j+1+2),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1),md(j-1+2),k).eq.6) then
                  system(mb(i+1),md(j-1+2),k) = 1
                  if(system(mb(i+1),md(j-1+2),ml(k-1)).eq.5) then
                     system(mb(i+1),md(j-1+2),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(i,md(j-2+2),k).eq.6) then
                  system(i,md(j-2+2),k) = 1
                  if(system(i,md(j-2+2),ml(k-1)).eq.5) then
                     system(i,md(j-2+2),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1),md(j-1+2),k).eq.6) then
                  system(mb(i-1),md(j-1+2),k) = 1
                  if(system(mb(i-1),md(j-1+2),ml(k-1)).eq.5) then
                     system(mb(i-1),md(j-1+2),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j+1+2),k).eq.6) then
                  system(mb(i-1),md(j+1+2),k) = 1
                  if(system(mb(i-1),md(j+1+2),ml(k-1)).eq.5) then
                     system(mb(i-1),md(j+1+2),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 333   
               end if

               if(up) then
                  if(system(i,md(j+4+2),ml(k-1)).eq.5) then
                     system(i,md(j+4+2),ml(k-1)) = 1
                  else if(system(mb(i+2),md(j-2+2),ml(k-1)).eq.5) then
                     system(mb(i+2),md(j-2+2),ml(k-1)) = 1
                  else if(system(mb(i-2),md(j-2+2),ml(k-1)).eq.5) then
                     system(mb(i-2),md(j-2+2),ml(k-1)) = 1
                  else
                     write(*,*)'error in move 26 - 6 - 1'
                  end if
                  go to 333
               else
                  if(system(i,md(j-4+2),ml(k-1)).eq.5) then
                     system(i,md(j-4+2),ml(k-1)) = 1
                  else if(system(mb(i-2),md(j+2+2),ml(k-1)).eq.5) then
                     system(mb(i-2),md(j+2+2),ml(k-1)) = 1
                  else if(system(mb(i+2),md(j+2+2),ml(k-1)).eq.5) then
                     system(mb(i+2),md(j+2+2),ml(k-1)) = 1
                  else
                     write(*,*)'error in move 26 - 6 - 2'
                  end if
                  go to 333
               end if               



            else if(system(i,md(j+2),ml(k+1)).eq.2) then
               system(i,md(j+2),ml(k+1)) = 1

c - - - - - - -close interlayer divac
               if(system(i,md(j+2+2),ml(k+1)).eq.5) then
                  system(i,md(j+2+2),ml(k+1)) = 1
                  if(system(i,md(j+2+2),ml(k+2)).eq.6) then
                     system(i,md(j+2+2),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j+1+2),ml(k+1)).eq.5) then
                  system(mb(i+1),md(j+1+2),ml(k+1)) = 1
                  if(system(mb(i+1),md(j+1+2),ml(k+2)).eq.6) then
                     system(mb(i+1),md(j+1+2),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1),md(j-1+2),ml(k+1)).eq.5) then
                  system(mb(i+1),md(j-1+2),ml(k+1)) = 1
                  if(system(mb(i+1),md(j-1+2),ml(k+2)).eq.6) then
                     system(mb(i+1),md(j-1+2),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(i,md(j-2+2),ml(k+1)).eq.5) then
                  system(i,md(j-2+2),ml(k+1)) = 1
                  if(system(i,md(j-2+2),ml(k+2)).eq.6) then
                     system(i,md(j-2+2),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1),md(j-1+2),ml(k+1)).eq.5) then
                  system(mb(i-1),md(j-1+2),ml(k+1)) = 1
                  if(system(mb(i-1),md(j-1+2),ml(k+2)).eq.6) then
                     system(mb(i-1),md(j-1+2),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j+1+2),ml(k+1)).eq.5) then
                  system(mb(i-1),md(j+1+2),ml(k+1)) = 1
                  if(system(mb(i-1),md(j+1+2),ml(k+2)).eq.6) then
                     system(mb(i-1),md(j+1+2),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 504
               end if

               if(up) then
                  if(system(i,md(j+4+2),ml(k+2)).eq.6) then
                     system(i,md(j+4+2),ml(k+2)) = 1
                  else if(system(mb(i+2),md(j-2+2),ml(k+2)).eq.6) then
                     system(mb(i+2),md(j-2+2),ml(k+2)) = 1
                  else if(system(mb(i-2),md(j-2+2),ml(k+2)).eq.6) then
                     system(mb(i-2),md(j-2+2),ml(k+2)) = 1
                  else
                     write(*,*)'error in move 26 - 3 - 1'
                  end if
                  go to 333
               else
                  if(system(i,md(j-4+2),ml(k+2)).eq.6) then
                     system(i,md(j-4+2),ml(k+2)) = 1
                  else if(system(mb(i-2),md(j+2+2),ml(k+2)).eq.6) then
                     system(mb(i-2),md(j+2+2),ml(k+2)) = 1
                  else if(system(mb(i+2),md(j+2+2),ml(k+2)).eq.6) then
                     system(mb(i+2),md(j+2+2),ml(k+2)) = 1
                  else
                     write(*,*)'error in move 26 - 3 - 2'
                  end if
                  go to 333
               end if



 504           if(system(i,md(j+2+2),ml(k+1)).eq.6) then
                  system(i,md(j+2+2),ml(k+1)) = 1
                  if(system(i,md(j+2+2),ml(k)).eq.5) then
                     system(i,md(j+2+2),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j+1+2),ml(k+1)).eq.6) then
                  system(mb(i+1),md(j+1+2),ml(k+1)) = 1
                  if(system(mb(i+1),md(j+1+2),ml(k)).eq.5) then
                     system(mb(i+1),md(j+1+2),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1),md(j-1+2),ml(k+1)).eq.6) then
                  system(mb(i+1),md(j-1+2),ml(k+1)) = 1
                  if(system(mb(i+1),md(j-1+2),ml(k)).eq.5) then
                     system(mb(i+1),md(j-1+2),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(i,md(j-2+2),ml(k+1)).eq.6) then
                  system(i,md(j-2+2),ml(k+1)) = 1
                  if(system(i,md(j-2+2),ml(k)).eq.5) then
                     system(i,md(j-2+2),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1),md(j-1+2),ml(k+1)).eq.6) then
                  system(mb(i-1),md(j-1+2),ml(k+1)) = 1
                  if(system(mb(i-1),md(j-1+2),ml(k)).eq.5) then
                     system(mb(i-1),md(j-1+2),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j+1+2),ml(k+1)).eq.6) then
                  system(mb(i-1),md(j+1+2),ml(k+1)) = 1
                  if(system(mb(i-1),md(j+1+2),ml(k)).eq.5) then
                     system(mb(i-1),md(j+1+2),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 333
               end if

               if(up) then
                  if(system(i,md(j+4+2),ml(k)).eq.5) then
                     system(i,md(j+4+2),ml(k)) = 1
                  else if(system(mb(i+2),md(j-2+2),ml(k)).eq.5) then
                     system(mb(i+2),md(j-2+2),ml(k)) = 1
                  else if(system(mb(i-2),md(j-2+2),ml(k)).eq.5) then
                     system(mb(i-2),md(j-2+2),ml(k)) = 1
                  else
                     write(*,*)'error in move 26 - 4 - 1'
                  end if
                  go to 333
               else
                  if(system(i,md(j-4+2),ml(k)).eq.5) then
                     system(i,md(j-4+2),ml(k)) = 1
                  else if(system(mb(i-2),md(j+2+2),ml(k)).eq.5) then
                     system(mb(i-2),md(j+2+2),ml(k)) = 1
                  else if(system(mb(i+2),md(j+2+2),ml(k)).eq.5) then
                     system(mb(i+2),md(j+2+2),ml(k)) = 1
                  else
                     write(*,*)'error in move 26 - 4 - 2'
                  end if
                  go to 333
               end if               


c..............intimate dir 2
            else if(system(mb(i+1),md(j+1),k).eq.2) then
               system(mb(i+1),md(j+1),k) = 1

c - - - - - - -close interlayer divac

               if(system(mb(i+1),md(j+2+1),k).eq.5) then
                  system(mb(i+1),md(j+2+1),k) = 1
                  if(system(mb(i+1),md(j+2+1),ml(k+1)).eq.6) then
                     system(mb(i+1),md(j+2+1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1+1),md(j+1+1),k).eq.5) then
                  system(mb(i+1+1),md(j+1+1),k) = 1
                  if(system(mb(i+1+1),md(j+1+1),ml(k+1)).eq.6) then
                     system(mb(i+1+1),md(j+1+1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1+1),md(j-1+1),k).eq.5) then
                  system(mb(i+1+1),md(j-1+1),k) = 1
                  if(system(mb(i+1+1),md(j-1+1),ml(k+1)).eq.6) then
                     system(mb(i+1+1),md(j-1+1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j-2+1),k).eq.5) then
                  system(mb(i+1),md(j-2+1),k) = 1
                  if(system(mb(i+1),md(j-2+1),ml(k+1)).eq.6) then
                     system(mb(i+1),md(j-2+1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1+1),md(j-1+1),k).eq.5) then
                  system(mb(i-1+1),md(j-1+1),k) = 1
                  if(system(mb(i-1+1),md(j-1+1),ml(k+1)).eq.6) then
                     system(mb(i-1+1),md(j-1+1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1+1),md(j+1+1),k).eq.5) then
                  system(mb(i-1+1),md(j+1+1),k) = 1
                  if(system(mb(i-1+1),md(j+1+1),ml(k+1)).eq.6) then
                     system(mb(i-1+1),md(j+1+1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 505
               end if

               if(up) then
                  if(system(mb(i+1),md(j+4+1),ml(k+1)).eq.6) then
                     system(mb(i+1),md(j+4+1),ml(k+1)) = 1
                  else if(system(mb(i+2+1),md(j-2+1),ml(k+1)).eq.6) then
                     system(mb(i+2+1),md(j-2+1),ml(k+1)) = 1
                  else if(system(mb(i-2+1),md(j-2+1),ml(k+1)).eq.6) then
                     system(mb(i-2+1),md(j-2+1),ml(k+1)) = 1
                  else
                     write(*,*)'error in move 26 - 1 - 1'
                  end if
                  go to 333
               else
                  if(system(mb(i+1),md(j-4+1),ml(k+1)).eq.6) then
                     system(mb(i+1),md(j-4+1),ml(k+1)) = 1
                  else if(system(mb(i-2+1),md(j+2+1),ml(k+1)).eq.6) then
                     system(mb(i-2+1),md(j+2+1),ml(k+1)) = 1
                  else if(system(mb(i+2+1),md(j+2+1),ml(k+1)).eq.6) then
                     system(mb(i+2+1),md(j+2+1),ml(k+1)) = 1
                  else
                     write(*,*)'error in move 26 - 1 - 2'
                  end if
                  go to 333
               end if



 505           if(system(mb(i+1),md(j+2+1),k).eq.6) then
                  system(mb(i+1),md(j+2+1),k) = 1
                  if(system(mb(i+1),md(j+2+1),ml(k-1)).eq.5) then
                     system(mb(i+1),md(j+2+1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1+1),md(j+1+1),k).eq.6) then
                  system(mb(i+1+1),md(j+1+1),k) = 1
                  if(system(mb(i+1+1),md(j+1+1),ml(k-1)).eq.5) then
                     system(mb(i+1+1),md(j+1+1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1+1),md(j-1+1),k).eq.6) then
                  system(mb(i+1+1),md(j-1+1),k) = 1
                  if(system(mb(i+1+1),md(j-1+1),ml(k-1)).eq.5) then
                     system(mb(i+1+1),md(j-1+1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j-2+1),k).eq.6) then
                  system(mb(i+1),md(j-2+1),k) = 1
                  if(system(mb(i+1),md(j-2+1),ml(k-1)).eq.5) then
                     system(mb(i+1),md(j-2+1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1+1),md(j-1+1),k).eq.6) then
                  system(mb(i-1+1),md(j-1+1),k) = 1
                  if(system(mb(i-1+1),md(j-1+1),ml(k-1)).eq.5) then
                     system(mb(i-1+1),md(j-1+1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1+1),md(j+1+1),k).eq.6) then
                  system(mb(i-1+1),md(j+1+1),k) = 1
                  if(system(mb(i-1+1),md(j+1+1),ml(k-1)).eq.5) then
                     system(mb(i-1+1),md(j+1+1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 333
               end if

               if(up) then
                  if(system(mb(i+1),md(j+4+1),ml(k-1)).eq.5) then
                     system(mb(i+1),md(j+4+1),ml(k-1)) = 1
                  else if(system(mb(i+2+1),md(j-2+1),ml(k-1)).eq.5) then
                     system(mb(i+2+1),md(j-2+1),ml(k-1)) = 1
                  else if(system(mb(i-2+1),md(j-2+1),ml(k-1)).eq.5) then
                     system(mb(i-2+1),md(j-2+1),ml(k-1)) = 1
                  else
                     write(*,*)'error in move 26 - 2 - 1'
                  end if
                  go to 333
               else
                  if(system(mb(i+1),md(j-4+1),ml(k-1)).eq.5) then
                     system(mb(i+1),md(j-4+1),ml(k-1)) = 1
                  else if(system(mb(i-2+1),md(j+2+1),ml(k-1)).eq.5) then
                     system(mb(i-2+1),md(j+2+1),ml(k-1)) = 1
                  else if(system(mb(i+2+1),md(j+2+1),ml(k-1)).eq.5) then
                     system(mb(i+2+1),md(j+2+1),ml(k-1)) = 1
                  else
                     write(*,*)'error in move 26 - 2 - 2'
                  end if
                  go to 333
               end if               



            else if(system(mb(i+1),md(j+1),ml(k+1)).eq.2) then
               system(mb(i+1),md(j+1),ml(k+1)) = 1

c - - - - - - -close interlayer divac

               if(system(mb(i+1),md(j+2+1),ml(k+1)).eq.5) then
                  system(mb(i+1),md(j+2+1),ml(k+1)) = 1
                  if(system(mb(i+1),md(j+2+1),ml(k+2)).eq.6) then
                     system(mb(i+1),md(j+2+1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1+1),md(j+1+1),ml(k+1)).eq.5) then
                  system(mb(i+1+1),md(j+1+1),ml(k+1)) = 1
                  if(system(mb(i+1+1),md(j+1+1),ml(k+2)).eq.6) then
                     system(mb(i+1+1),md(j+1+1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1+1),md(j-1+1),ml(k+1)).eq.5) then
                  system(mb(i+1+1),md(j-1+1),ml(k+1)) = 1
                  if(system(mb(i+1+1),md(j-1+1),ml(k+2)).eq.6) then
                     system(mb(i+1+1),md(j-1+1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j-2+1),ml(k+1)).eq.5) then
                  system(mb(i+1),md(j-2+1),ml(k+1)) = 1
                  if(system(mb(i+1),md(j-2+1),ml(k+2)).eq.6) then
                     system(mb(i+1),md(j-2+1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1+1),md(j-1+1),ml(k+1)).eq.5) then
                  system(mb(i-1+1),md(j-1+1),ml(k+1)) = 1
                  if(system(mb(i-1+1),md(j-1+1),ml(k+2)).eq.6) then
                     system(mb(i-1+1),md(j-1+1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1+1),md(j+1+1),ml(k+1)).eq.5) then
                  system(mb(i-1+1),md(j+1+1),ml(k+1)) = 1
                  if(system(mb(i-1+1),md(j+1+1),ml(k+2)).eq.6) then
                     system(mb(i-1+1),md(j+1+1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 506
               end if

               if(up) then
                  if(system(mb(i+1),md(j+4+1),ml(k+2)).eq.6) then
                     system(mb(i+1),md(j+4+1),ml(k+2)) = 1
                  else if(system(mb(i+2+1),md(j-2+1),ml(k+2)).eq.6) then
                     system(mb(i+2+1),md(j-2+1),ml(k+2)) = 1
                  else if(system(mb(i-2+1),md(j-2+1),ml(k+2)).eq.6) then
                     system(mb(i-2+1),md(j-2+1),ml(k+2)) = 1
                  else
                     write(*,*)'error in move 26 - 3 - 1'
                  end if
                  go to 333
               else
                  if(system(mb(i+1),md(j-4+1),ml(k+2)).eq.6) then
                     system(mb(i+1),md(j-4+1),ml(k+2)) = 1
                  else if(system(mb(i-2+1),md(j+2+1),ml(k+2)).eq.6) then
                     system(mb(i-2+1),md(j+2+1),ml(k+2)) = 1
                  else if(system(mb(i+2+1),md(j+2+1),ml(k+2)).eq.6) then
                     system(mb(i+2+1),md(j+2+1),ml(k+2)) = 1
                  else
                     write(*,*)'error in move 26 - 3 - 2'
                  end if
                  go to 333
               end if



 506           if(system(mb(i+1),md(j+2+1),ml(k+1)).eq.6) then
                  system(mb(i+1),md(j+2+1),ml(k+1)) = 1
                  if(system(mb(i+1),md(j+2+1),ml(k)).eq.5) then
                     system(mb(i+1),md(j+2+1),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1+1),md(j+1+1),ml(k+1)).eq.6) then
                  system(mb(i+1+1),md(j+1+1),ml(k+1)) = 1
                  if(system(mb(i+1+1),md(j+1+1),ml(k)).eq.5) then
                     system(mb(i+1+1),md(j+1+1),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1+1),md(j-1+1),ml(k+1)).eq.6) then
                  system(mb(i+1+1),md(j-1+1),ml(k+1)) = 1
                  if(system(mb(i+1+1),md(j-1+1),ml(k)).eq.5) then
                     system(mb(i+1+1),md(j-1+1),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j-2+1),ml(k+1)).eq.6) then
                  system(mb(i+1),md(j-2+1),ml(k+1)) = 1
                  if(system(mb(i+1),md(j-2+1),ml(k)).eq.5) then
                     system(mb(i+1),md(j-2+1),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1+1),md(j-1+1),ml(k+1)).eq.6) then
                  system(mb(i-1+1),md(j-1+1),ml(k+1)) = 1
                  if(system(mb(i-1+1),md(j-1+1),ml(k)).eq.5) then
                     system(mb(i-1+1),md(j-1+1),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1+1),md(j+1+1),ml(k+1)).eq.6) then
                  system(mb(i-1+1),md(j+1+1),ml(k+1)) = 1
                  if(system(mb(i-1+1),md(j+1+1),ml(k)).eq.5) then
                     system(mb(i-1+1),md(j+1+1),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 333
               end if

               if(up) then
                  if(system(mb(i+1),md(j+4+1),ml(k)).eq.5) then
                     system(mb(i+1),md(j+4+1),ml(k)) = 1
                  else if(system(mb(i+2+1),md(j-2+1),ml(k)).eq.5) then
                     system(mb(i+2+1),md(j-2+1),ml(k)) = 1
                  else if(system(mb(i-2+1),md(j-2+1),ml(k)).eq.5) then
                     system(mb(i-2+1),md(j-2+1),ml(k)) = 1
                  else
                     write(*,*)'error in move 26 - 4 - 1'
                  end if
                  go to 333
               else
                  if(system(mb(i+1),md(j-4+1),ml(k)).eq.5) then
                     system(mb(i+1),md(j-4+1),ml(k)) = 1
                  else if(system(mb(i-2+1),md(j+2+1),ml(k)).eq.5) then
                     system(mb(i-2+1),md(j+2+1),ml(k)) = 1
                  else if(system(mb(i+2+1),md(j+2+1),ml(k)).eq.5) then
                     system(mb(i+2+1),md(j+2+1),ml(k)) = 1
                  else
                     write(*,*)'error in move 26 - 4 - 2'
                  end if
                  go to 333
               end if               

c..............intimate dir 3
            else if(system(mb(i+1),md(j-1),k).eq.2) then
               system(mb(i+1),md(j-1),k) = 1

c - - - - - - -close interlayer divac

               if(system(mb(i+1),md(j+2-1),k).eq.5) then
                  system(mb(i+1),md(j+2-1),k) = 1
                  if(system(mb(i+1),md(j+2-1),ml(k+1)).eq.6) then
                     system(mb(i+1),md(j+2-1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1+1),md(j+1-1),k).eq.5) then
                  system(mb(i+1+1),md(j+1-1),k) = 1
                  if(system(mb(i+1+1),md(j+1-1),ml(k+1)).eq.6) then
                     system(mb(i+1+1),md(j+1-1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1+1),md(j-1-1),k).eq.5) then
                  system(mb(i+1+1),md(j-1-1),k) = 1
                  if(system(mb(i+1+1),md(j-1-1),ml(k+1)).eq.6) then
                     system(mb(i+1+1),md(j-1-1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j-2-1),k).eq.5) then
                  system(mb(i+1),md(j-2-1),k) = 1
                  if(system(mb(i+1),md(j-2-1),ml(k+1)).eq.6) then
                     system(mb(i+1),md(j-2-1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1+1),md(j-1-1),k).eq.5) then
                  system(mb(i-1+1),md(j-1-1),k) = 1
                  if(system(mb(i-1+1),md(j-1-1),ml(k+1)).eq.6) then
                     system(mb(i-1+1),md(j-1-1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1+1),md(j+1-1),k).eq.5) then
                  system(mb(i-1+1),md(j+1-1),k) = 1
                  if(system(mb(i-1+1),md(j+1-1),ml(k+1)).eq.6) then
                     system(mb(i-1+1),md(j+1-1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 507
               end if

               if(up) then
                  if(system(mb(i+1),md(j+4-1),ml(k+1)).eq.6) then
                     system(mb(i+1),md(j+4-1),ml(k+1)) = 1
                  else if(system(mb(i+2+1),md(j-2-1),ml(k+1)).eq.6) then
                     system(mb(i+2+1),md(j-2-1),ml(k+1)) = 1
                  else if(system(mb(i-2+1),md(j-2-1),ml(k+1)).eq.6) then
                     system(mb(i-2+1),md(j-2-1),ml(k+1)) = 1
                  else
                     write(*,*)'error in move 26 - 1 - 1'
                  end if
                  go to 333
               else
                  if(system(mb(i+1),md(j-4-1),ml(k+1)).eq.6) then
                     system(mb(i+1),md(j-4-1),ml(k+1)) = 1
                  else if(system(mb(i-2+1),md(j+2-1),ml(k+1)).eq.6) then
                     system(mb(i-2+1),md(j+2-1),ml(k+1)) = 1
                  else if(system(mb(i+2+1),md(j+2-1),ml(k+1)).eq.6) then
                     system(mb(i+2+1),md(j+2-1),ml(k+1)) = 1
                  else
                     write(*,*)'error in move 26 - 1 - 2'
                  end if
                  go to 333
               end if



 507           if(system(mb(i+1),md(j+2-1),k).eq.6) then
                  system(mb(i+1),md(j+2-1),k) = 1
                  if(system(mb(i+1),md(j+2-1),ml(k-1)).eq.5) then
                     system(mb(i+1),md(j+2-1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1+1),md(j+1-1),k).eq.6) then
                  system(mb(i+1+1),md(j+1-1),k) = 1
                  if(system(mb(i+1+1),md(j+1-1),ml(k-1)).eq.5) then
                     system(mb(i+1+1),md(j+1-1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1+1),md(j-1-1),k).eq.6) then
                  system(mb(i+1+1),md(j-1-1),k) = 1
                  if(system(mb(i+1+1),md(j-1-1),ml(k-1)).eq.5) then
                     system(mb(i+1+1),md(j-1-1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j-2-1),k).eq.6) then
                  system(mb(i+1),md(j-2-1),k) = 1
                  if(system(mb(i+1),md(j-2-1),ml(k-1)).eq.5) then
                     system(mb(i+1),md(j-2-1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1+1),md(j-1-1),k).eq.6) then
                  system(mb(i-1+1),md(j-1-1),k) = 1
                  if(system(mb(i-1+1),md(j-1-1),ml(k-1)).eq.5) then
                     system(mb(i-1+1),md(j-1-1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1+1),md(j+1-1),k).eq.6) then
                  system(mb(i-1+1),md(j+1-1),k) = 1
                  if(system(mb(i-1+1),md(j+1-1),ml(k-1)).eq.5) then
                     system(mb(i-1+1),md(j+1-1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 333
               end if

               if(up) then
                  if(system(mb(i+1),md(j+4-1),ml(k-1)).eq.5) then
                     system(mb(i+1),md(j+4-1),ml(k-1)) = 1
                  else if(system(mb(i+2+1),md(j-2-1),ml(k-1)).eq.5) then
                     system(mb(i+2+1),md(j-2-1),ml(k-1)) = 1
                  else if(system(mb(i-2+1),md(j-2-1),ml(k-1)).eq.5) then
                     system(mb(i-2+1),md(j-2-1),ml(k-1)) = 1
                  else
                     write(*,*)'error in move 26 - 2 - 1'
                  end if
                  go to 333
               else
                  if(system(mb(i+1),md(j-4-1),ml(k-1)).eq.5) then
                     system(mb(i+1),md(j-4-1),ml(k-1)) = 1
                  else if(system(mb(i-2+1),md(j+2-1),ml(k-1)).eq.5) then
                     system(mb(i-2+1),md(j+2-1),ml(k-1)) = 1
                  else if(system(mb(i+2+1),md(j+2+1),ml(k-1)).eq.5) then
                     system(mb(i+2+1),md(j+2-1),ml(k-1)) = 1
                  else
                     write(*,*)'error in move 26 - 2 - 2'
                  end if
                  go to 333
               end if               

            else if(system(mb(i+1),md(j-1),ml(k+1)).eq.2) then
               system(mb(i+1),md(j-1),ml(k+1)) = 1

c - - - - - - -close interlayer divac

               if(system(mb(i+1),md(j+2-1),ml(k+1)).eq.5) then
                  system(mb(i+1),md(j+2-1),ml(k+1)) = 1
                  if(system(mb(i+1),md(j+2-1),ml(k+2)).eq.6) then
                     system(mb(i+1),md(j+2-1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1+1),md(j+1-1),ml(k+1)).eq.5) then
                  system(mb(i+1+1),md(j+1-1),ml(k+1)) = 1
                  if(system(mb(i+1+1),md(j+1-1),ml(k+2)).eq.6) then
                     system(mb(i+1+1),md(j+1-1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1+1),md(j-1-1),ml(k+1)).eq.5) then
                  system(mb(i+1+1),md(j-1-1),ml(k+1)) = 1
                  if(system(mb(i+1+1),md(j-1-1),ml(k+2)).eq.6) then
                     system(mb(i+1+1),md(j-1-1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j-2-1),ml(k+1)).eq.5) then
                  system(mb(i+1),md(j-2-1),ml(k+1)) = 1
                  if(system(mb(i+1),md(j-2-1),ml(k+2)).eq.6) then
                     system(mb(i+1),md(j-2-1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1+1),md(j-1-1),ml(k+1)).eq.5) then
                  system(mb(i-1+1),md(j-1-1),ml(k+1)) = 1
                  if(system(mb(i-1+1),md(j-1-1),ml(k+2)).eq.6) then
                     system(mb(i-1+1),md(j-1-1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1+1),md(j+1-1),ml(k+1)).eq.5) then
                  system(mb(i-1+1),md(j+1-1),ml(k+1)) = 1
                  if(system(mb(i-1+1),md(j+1-1),ml(k+2)).eq.6) then
                     system(mb(i-1+1),md(j+1-1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 508
               end if

               if(up) then
                  if(system(mb(i+1),md(j+4-1),ml(k+2)).eq.6) then
                     system(mb(i+1),md(j+4-1),ml(k+2)) = 1
                  else if(system(mb(i+2+1),md(j-2-1),ml(k+2)).eq.6) then
                     system(mb(i+2+1),md(j-2-1),ml(k+2)) = 1
                  else if(system(mb(i-2+1),md(j-2-1),ml(k+2)).eq.6) then
                     system(mb(i-2+1),md(j-2-1),ml(k+2)) = 1
                  else
                     write(*,*)'error in move 26 - 3 - 1'
                  end if
                  go to 333
               else
                  if(system(mb(i+1),md(j-4-1),ml(k+2)).eq.6) then
                     system(mb(i+1),md(j-4-1),ml(k+2)) = 1
                  else if(system(mb(i-2+1),md(j+2-1),ml(k+2)).eq.6) then
                     system(mb(i-2+1),md(j+2-1),ml(k+2)) = 1
                  else if(system(mb(i+2+1),md(j+2-1),ml(k+2)).eq.6) then
                     system(mb(i+2+1),md(j+2-1),ml(k+2)) = 1
                  else
                     write(*,*)'error in move 26 - 3 - 2'
                  end if
                  go to 333
               end if



 508           if(system(mb(i+1),md(j+2-1),ml(k+1)).eq.6) then
                  system(mb(i+1),md(j+2-1),ml(k+1)) = 1
                  if(system(mb(i+1),md(j+2-1),ml(k)).eq.5) then
                     system(mb(i+1),md(j+2-1),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1+1),md(j+1-1),ml(k+1)).eq.6) then
                  system(mb(i+1+1),md(j+1-1),ml(k+1)) = 1
                  if(system(mb(i+1+1),md(j+1-1),ml(k)).eq.5) then
                     system(mb(i+1+1),md(j+1-1),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1+1),md(j-1-1),ml(k+1)).eq.6) then
                  system(mb(i+1+1),md(j-1-1),ml(k+1)) = 1
                  if(system(mb(i+1+1),md(j-1-1),ml(k)).eq.5) then
                     system(mb(i+1+1),md(j-1-1),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j-2-1),ml(k+1)).eq.6) then
                  system(mb(i+1),md(j-2-1),ml(k+1)) = 1
                  if(system(mb(i+1),md(j-2-1),ml(k)).eq.5) then
                     system(mb(i+1),md(j-2-1),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1+1),md(j-1-1),ml(k+1)).eq.6) then
                  system(mb(i-1+1),md(j-1-1),ml(k+1)) = 1
                  if(system(mb(i-1+1),md(j-1-1),ml(k)).eq.5) then
                     system(mb(i-1+1),md(j-1-1),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1+1),md(j+1-1),ml(k+1)).eq.6) then
                  system(mb(i-1+1),md(j+1-1),ml(k+1)) = 1
                  if(system(mb(i-1+1),md(j+1-1),ml(k)).eq.5) then
                     system(mb(i-1+1),md(j+1-1),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 333
               end if

               if(up) then
                  if(system(mb(i+1),md(j+4-1),ml(k)).eq.5) then
                     system(mb(i+1),md(j+4-1),ml(k)) = 1
                  else if(system(mb(i+2+1),md(j-2-1),ml(k)).eq.5) then
                     system(mb(i+2+1),md(j-2-1),ml(k)) = 1
                  else if(system(mb(i-2+1),md(j-2-1),ml(k)).eq.5) then
                     system(mb(i-2+1),md(j-2-1),ml(k)) = 1
                  else
                     write(*,*)'error in move 26 - 4 - 1'
                  end if
                  go to 333
               else
                  if(system(mb(i+1),md(j-4-1),ml(k)).eq.5) then
                     system(mb(i+1),md(j-4-1),ml(k)) = 1
                  else if(system(mb(i-2+1),md(j+2-1),ml(k)).eq.5) then
                     system(mb(i-2+1),md(j+2-1),ml(k)) = 1
                  else if(system(mb(i+2+1),md(j+2-1),ml(k)).eq.5) then
                     system(mb(i+2+1),md(j+2-1),ml(k)) = 1
                  else
                     write(*,*)'error in move 26 - 4 - 2'
                  end if
                  go to 333
               end if               

c..............intimate dir 4
            else if(system(i,md(j-2),k).eq.2) then
               system(i,md(j-2),k) = 1

c - - - - - - -close interlayer divac
               if(system(i,md(j+2-2),k).eq.5) then
                  system(i,md(j+2-2),k) = 1
                  if(system(i,md(j+2-2),ml(k+1)).eq.6) then
                     system(i,md(j+2-2),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j+1-2),k).eq.5) then
                  system(mb(i+1),md(j+1-2),k) = 1
                  if(system(mb(i+1),md(j+1-2),ml(k+1)).eq.6) then
                     system(mb(i+1),md(j+1-2),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1),md(j-1-2),k).eq.5) then
                  system(mb(i+1),md(j-1-2),k) = 1
                  if(system(mb(i+1),md(j-1-2),ml(k+1)).eq.6) then
                     system(mb(i+1),md(j-1-2),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(i,md(j-2-2),k).eq.5) then
                  system(i,md(j-2-2),k) = 1
                  if(system(i,md(j-2-2),ml(k+1)).eq.6) then
                     system(i,md(j-2-2),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1),md(j-1-2),k).eq.5) then
                  system(mb(i-1),md(j-1-2),k) = 1
                  if(system(mb(i-1),md(j-1-2),ml(k+1)).eq.6) then
                     system(mb(i-1),md(j-1-2),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j+1-2),k).eq.5) then
                  system(mb(i-1),md(j+1-2),k) = 1
                  if(system(mb(i-1),md(j+1-2),ml(k+1)).eq.6) then
                     system(mb(i-1),md(j+1-2),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 509
               end if

               if(up) then
                  if(system(i,md(j+4-2),ml(k+1)).eq.6) then
                     system(i,md(j+4-2),ml(k+1)) = 1
                  else if(system(mb(i+2),md(j-2-2),ml(k+1)).eq.6) then
                     system(mb(i+2),md(j-2-2),ml(k+1)) = 1
                  else if(system(mb(i-2),md(j-2-2),ml(k+1)).eq.6) then
                     system(mb(i-2),md(j-2-2),ml(k+1)) = 1
                  else
                     write(*,*)'error in move 26 - 5 - 1'
                  end if
                  go to 333
               else
                  if(system(i,md(j-4-2),ml(k+1)).eq.6) then
                     system(i,md(j-4-2),ml(k+1)) = 1
                  else if(system(mb(i-2),md(j+2-2),ml(k+1)).eq.6) then
                     system(mb(i-2),md(j+2-2),ml(k+1)) = 1
                  else if(system(mb(i+2),md(j+2-2),ml(k+1)).eq.6) then
                     system(mb(i+2),md(j+2-2),ml(k+1)) = 1
                  else
                     write(*,*)'error in move 26 - 5 - 2'
                  end if
                  go to 333
               end if



 509           if(system(i,md(j+2-2),k).eq.6) then
                  system(i,md(j+2-2),k) = 1
                  if(system(i,md(j+2-2),ml(k-1)).eq.5) then
                     system(i,md(j+2-2),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j+1-2),k).eq.6) then
                  system(mb(i+1),md(j+1-2),k) = 1
                  if(system(mb(i+1),md(j+1-2),ml(k-1)).eq.5) then
                     system(mb(i+1),md(j+1-2),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1),md(j-1-2),k).eq.6) then
                  system(mb(i+1),md(j-1-2),k) = 1
                  if(system(mb(i+1),md(j-1-2),ml(k-1)).eq.5) then
                     system(mb(i+1),md(j-1-2),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(i,md(j-2-2),k).eq.6) then
                  system(i,md(j-2-2),k) = 1
                  if(system(i,md(j-2-2),ml(k-1)).eq.5) then
                     system(i,md(j-2-2),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1),md(j-1-2),k).eq.6) then
                  system(mb(i-1),md(j-1-2),k) = 1
                  if(system(mb(i-1),md(j-1-2),ml(k-1)).eq.5) then
                     system(mb(i-1),md(j-1-2),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j+1-2),k).eq.6) then
                  system(mb(i-1),md(j+1-2),k) = 1
                  if(system(mb(i-1),md(j+1-2),ml(k-1)).eq.5) then
                     system(mb(i-1),md(j+1-2),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 333
               end if

               if(up) then
                  if(system(i,md(j+4-2),ml(k-1)).eq.5) then
                     system(i,md(j+4-2),ml(k-1)) = 1
                  else if(system(mb(i+2),md(j-2-2),ml(k-1)).eq.5) then
                     system(mb(i+2),md(j-2-2),ml(k-1)) = 1
                  else if(system(mb(i-2),md(j-2-2),ml(k-1)).eq.5) then
                     system(mb(i-2),md(j-2-2),ml(k-1)) = 1
                  else
                     write(*,*)'error in move 26 - 6 - 1'
                  end if
                  go to 333
               else
                  if(system(i,md(j-4-2),ml(k-1)).eq.5) then
                     system(i,md(j-4-2),ml(k-1)) = 1
                  else if(system(mb(i-2),md(j+2-2),ml(k-1)).eq.5) then
                     system(mb(i-2),md(j+2-2),ml(k-1)) = 1
                  else if(system(mb(i+2),md(j+2-2),ml(k-1)).eq.5) then
                     system(mb(i+2),md(j+2-2),ml(k-1)) = 1
                  else
                     write(*,*)'error in move 26 - 6 - 2'
                  end if
                  go to 333
               end if               

            else if(system(i,md(j-2),ml(k+1)).eq.2) then
               system(i,md(j-2),ml(k+1)) = 1

c - - - - - - -close interlayer divac
               if(system(i,md(j+2-2),ml(k+1)).eq.5) then
                  system(i,md(j+2-2),ml(k+1)) = 1
                  if(system(i,md(j+2-2),ml(k+2)).eq.6) then
                     system(i,md(j+2-2),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j+1-2),ml(k+1)).eq.5) then
                  system(mb(i+1),md(j+1-2),ml(k+1)) = 1
                  if(system(mb(i+1),md(j+1-2),ml(k+2)).eq.6) then
                     system(mb(i+1),md(j+1-2),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1),md(j-1-2),ml(k+1)).eq.5) then
                  system(mb(i+1),md(j-1-2),ml(k+1)) = 1
                  if(system(mb(i+1),md(j-1-2),ml(k+2)).eq.6) then
                     system(mb(i+1),md(j-1-2),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(i,md(j-2-2),ml(k+1)).eq.5) then
                  system(i,md(j-2-2),ml(k+1)) = 1
                  if(system(i,md(j-2-2),ml(k+2)).eq.6) then
                     system(i,md(j-2-2),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1),md(j-1-2),ml(k+1)).eq.5) then
                  system(mb(i-1),md(j-1-2),ml(k+1)) = 1
                  if(system(mb(i-1),md(j-1-2),ml(k+2)).eq.6) then
                     system(mb(i-1),md(j-1-2),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j+1-2),ml(k+1)).eq.5) then
                  system(mb(i-1),md(j+1-2),ml(k+1)) = 1
                  if(system(mb(i-1),md(j+1-2),ml(k+2)).eq.6) then
                     system(mb(i-1),md(j+1-2),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 510
               end if

               if(up) then
                  if(system(i,md(j+4-2),ml(k+2)).eq.6) then
                     system(i,md(j+4-2),ml(k+2)) = 1
                  else if(system(mb(i+2),md(j-2-2),ml(k+2)).eq.6) then
                     system(mb(i+2),md(j-2-2),ml(k+2)) = 1
                  else if(system(mb(i-2),md(j-2-2),ml(k+2)).eq.6) then
                     system(mb(i-2),md(j-2-2),ml(k+2)) = 1
                  else
                     write(*,*)'error in move 26 - 3 - 1'
                  end if
                  go to 333
               else
                  if(system(i,md(j-4-2),ml(k+2)).eq.6) then
                     system(i,md(j-4-2),ml(k+2)) = 1
                  else if(system(mb(i-2),md(j+2-2),ml(k+2)).eq.6) then
                     system(mb(i-2),md(j+2-2),ml(k+2)) = 1
                  else if(system(mb(i+2),md(j+2-2),ml(k+2)).eq.6) then
                     system(mb(i+2),md(j+2-2),ml(k+2)) = 1
                  else
                     write(*,*)'error in move 26 - 3 - 2'
                  end if
                  go to 333
               end if



 510           if(system(i,md(j+2-2),ml(k+1)).eq.6) then
                  system(i,md(j+2-2),ml(k+1)) = 1
                  if(system(i,md(j+2-2),ml(k)).eq.5) then
                     system(i,md(j+2-2),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1),md(j+1-2),ml(k+1)).eq.6) then
                  system(mb(i+1),md(j+1-2),ml(k+1)) = 1
                  if(system(mb(i+1),md(j+1-2),ml(k)).eq.5) then
                     system(mb(i+1),md(j+1-2),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1),md(j-1-2),ml(k+1)).eq.6) then
                  system(mb(i+1),md(j-1-2),ml(k+1)) = 1
                  if(system(mb(i+1),md(j-1-2),ml(k)).eq.5) then
                     system(mb(i+1),md(j-1-2),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(i,md(j-2-2),ml(k+1)).eq.6) then
                  system(i,md(j-2-2),ml(k+1)) = 1
                  if(system(i,md(j-2-2),ml(k)).eq.5) then
                     system(i,md(j-2-2),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1),md(j-1-2),ml(k+1)).eq.6) then
                  system(mb(i-1),md(j-1-2),ml(k+1)) = 1
                  if(system(mb(i-1),md(j-1-2),ml(k)).eq.5) then
                     system(mb(i-1),md(j-1-2),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j+1-2),ml(k+1)).eq.6) then
                  system(mb(i-1),md(j+1-2),ml(k+1)) = 1
                  if(system(mb(i-1),md(j+1-2),ml(k)).eq.5) then
                     system(mb(i-1),md(j+1-2),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 333
               end if

               if(up) then
                  if(system(i,md(j+4-2),ml(k)).eq.5) then
                     system(i,md(j+4-2),ml(k)) = 1
                  else if(system(mb(i+2),md(j-2-2),ml(k)).eq.5) then
                     system(mb(i+2),md(j-2-2),ml(k)) = 1
                  else if(system(mb(i-2),md(j-2-2),ml(k)).eq.5) then
                     system(mb(i-2),md(j-2-2),ml(k)) = 1
                  else
                     write(*,*)'error in move 26 - 4 - 1'
                  end if
                  go to 333
               else
                  if(system(i,md(j-4-2),ml(k)).eq.5) then
                     system(i,md(j-4-2),ml(k)) = 1
                  else if(system(mb(i-2),md(j+2-2),ml(k)).eq.5) then
                     system(mb(i-2),md(j+2-2),ml(k)) = 1
                  else if(system(mb(i+2),md(j+2-2),ml(k)).eq.5) then
                     system(mb(i+2),md(j+2-2),ml(k)) = 1
                  else
                     write(*,*)'error in move 26 - 4 - 2'
                  end if
                  go to 333
               end if               

c..............intimate dir 5
            else if(system(mb(i-1),md(j-1),k).eq.2) then
               system(mb(i-1),md(j-1),k) = 1

c - - - - - - -close interlayer divac

               if(system(mb(i-1),md(j+2-1),k).eq.5) then
                  system(mb(i-1),md(j+2-1),k) = 1
                  if(system(mb(i-1),md(j+2-1),ml(k+1)).eq.6) then
                     system(mb(i-1),md(j+2-1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1-1),md(j+1-1),k).eq.5) then
                  system(mb(i+1-1),md(j+1-1),k) = 1
                  if(system(mb(i+1-1),md(j+1-1),ml(k+1)).eq.6) then
                     system(mb(i+1-1),md(j+1-1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1-1),md(j-1-1),k).eq.5) then
                  system(mb(i+1-1),md(j-1-1),k) = 1
                  if(system(mb(i+1-1),md(j-1-1),ml(k+1)).eq.6) then
                     system(mb(i+1-1),md(j-1-1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j-2-1),k).eq.5) then
                  system(mb(i-1),md(j-2-1),k) = 1
                  if(system(mb(i-1),md(j-2-1),ml(k+1)).eq.6) then
                     system(mb(i-1),md(j-2-1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1-1),md(j-1-1),k).eq.5) then
                  system(mb(i-1-1),md(j-1-1),k) = 1
                  if(system(mb(i-1-1),md(j-1-1),ml(k+1)).eq.6) then
                     system(mb(i-1-1),md(j-1-1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1-1),md(j+1-1),k).eq.5) then
                  system(mb(i-1-1),md(j+1-1),k) = 1
                  if(system(mb(i-1-1),md(j+1-1),ml(k+1)).eq.6) then
                     system(mb(i-1-1),md(j+1-1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 511
               end if

               if(up) then
                  if(system(mb(i-1),md(j+4-1),ml(k+1)).eq.6) then
                     system(mb(i-1),md(j+4-1),ml(k+1)) = 1
                  else if(system(mb(i+2-1),md(j-2-1),ml(k+1)).eq.6) then
                     system(mb(i+2-1),md(j-2-1),ml(k+1)) = 1
                  else if(system(mb(i-2-1),md(j-2-1),ml(k+1)).eq.6) then
                     system(mb(i-2-1),md(j-2-1),ml(k+1)) = 1
                  else
                     write(*,*)'error in move 26 - 1 - 1'
                  end if
                  go to 333
               else
                  if(system(mb(i-1),md(j-4-1),ml(k+1)).eq.6) then
                     system(mb(i-1),md(j-4-1),ml(k+1)) = 1
                  else if(system(mb(i-2-1),md(j+2-1),ml(k+1)).eq.6) then
                     system(mb(i-2-1),md(j+2-1),ml(k+1)) = 1
                  else if(system(mb(i+2-1),md(j+2-1),ml(k+1)).eq.6) then
                     system(mb(i+2-1),md(j+2-1),ml(k+1)) = 1
                  else
                     write(*,*)'error in move 26 - 1 - 2'
                  end if
                  go to 333
               end if



 511           if(system(mb(i-1),md(j+2-1),k).eq.6) then
                  system(mb(i-1),md(j+2-1),k) = 1
                  if(system(mb(i-1),md(j+2-1),ml(k-1)).eq.5) then
                     system(mb(i-1),md(j+2-1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1-1),md(j+1-1),k).eq.6) then
                  system(mb(i+1-1),md(j+1-1),k) = 1
                  if(system(mb(i+1-1),md(j+1-1),ml(k-1)).eq.5) then
                     system(mb(i+1-1),md(j+1-1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1-1),md(j-1-1),k).eq.6) then
                  system(mb(i+1-1),md(j-1-1),k) = 1
                  if(system(mb(i+1-1),md(j-1-1),ml(k-1)).eq.5) then
                     system(mb(i+1-1),md(j-1-1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j-2-1),k).eq.6) then
                  system(mb(i-1),md(j-2-1),k) = 1
                  if(system(mb(i-1),md(j-2-1),ml(k-1)).eq.5) then
                     system(mb(i-1),md(j-2-1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1-1),md(j-1-1),k).eq.6) then
                  system(mb(i-1-1),md(j-1-1),k) = 1
                  if(system(mb(i-1-1),md(j-1-1),ml(k-1)).eq.5) then
                     system(mb(i-1-1),md(j-1-1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1-1),md(j+1-1),k).eq.6) then
                  system(mb(i-1-1),md(j+1-1),k) = 1
                  if(system(mb(i-1-1),md(j+1-1),ml(k-1)).eq.5) then
                     system(mb(i-1-1),md(j+1-1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 333
               end if

               if(up) then
                  if(system(mb(i-1),md(j+4-1),ml(k-1)).eq.5) then
                     system(mb(i-1),md(j+4-1),ml(k-1)) = 1
                  else if(system(mb(i+2-1),md(j-2-1),ml(k-1)).eq.5) then
                     system(mb(i+2-1),md(j-2-1),ml(k-1)) = 1
                  else if(system(mb(i-2+1),md(j-2-1),ml(k-1)).eq.5) then
                     system(mb(i-2-1),md(j-2-1),ml(k-1)) = 1
                  else
                     write(*,*)'error in move 26 - 2 - 1'
                  end if
                  go to 333
               else
                  if(system(mb(i-1),md(j-4-1),ml(k-1)).eq.5) then
                     system(mb(i-1),md(j-4-1),ml(k-1)) = 1
                  else if(system(mb(i-2-1),md(j+2-1),ml(k-1)).eq.5) then
                     system(mb(i-2-1),md(j+2-1),ml(k-1)) = 1
                  else if(system(mb(i+2-1),md(j+2-1),ml(k-1)).eq.5) then
                     system(mb(i+2-1),md(j+2-1),ml(k-1)) = 1
                  else
                     write(*,*)'error in move 26 - 2 - 2'
                  end if
                  go to 333
               end if               

            else if(system(mb(i-1),md(j-1),ml(k+1)).eq.2) then
               system(mb(i-1),md(j-1),ml(k+1)) = 1

c - - - - - - -close interlayer divac

               if(system(mb(i-1),md(j+2-1),ml(k+1)).eq.5) then
                  system(mb(i-1),md(j+2-1),ml(k+1)) = 1
                  if(system(mb(i-1),md(j+2-1),ml(k+2)).eq.6) then
                     system(mb(i-1),md(j+2-1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1-1),md(j+1-1),ml(k+1)).eq.5) then
                  system(mb(i+1-1),md(j+1-1),ml(k+1)) = 1
                  if(system(mb(i+1-1),md(j+1-1),ml(k+2)).eq.6) then
                     system(mb(i+1-1),md(j+1-1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1-1),md(j-1-1),ml(k+1)).eq.5) then
                  system(mb(i+1-1),md(j-1-1),ml(k+1)) = 1
                  if(system(mb(i+1-1),md(j-1-1),ml(k+2)).eq.6) then
                     system(mb(i+1-1),md(j-1-1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j-2-1),ml(k+1)).eq.5) then
                  system(mb(i-1),md(j-2-1),ml(k+1)) = 1
                  if(system(mb(i-1),md(j-2-1),ml(k+2)).eq.6) then
                     system(mb(i-1),md(j-2-1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1-1),md(j-1-1),ml(k+1)).eq.5) then
                  system(mb(i-1-1),md(j-1-1),ml(k+1)) = 1
                  if(system(mb(i-1-1),md(j-1-1),ml(k+2)).eq.6) then
                     system(mb(i-1-1),md(j-1-1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1-1),md(j+1-1),ml(k+1)).eq.5) then
                  system(mb(i-1-1),md(j+1-1),ml(k+1)) = 1
                  if(system(mb(i-1-1),md(j+1-1),ml(k+2)).eq.6) then
                     system(mb(i-1-1),md(j+1-1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 512
               end if

               if(up) then
                  if(system(mb(i-1),md(j+4-1),ml(k+2)).eq.6) then
                     system(mb(i-1),md(j+4-1),ml(k+2)) = 1
                  else if(system(mb(i+2-1),md(j-2-1),ml(k+2)).eq.6) then
                     system(mb(i+2-1),md(j-2-1),ml(k+2)) = 1
                  else if(system(mb(i-2+1),md(j-2-1),ml(k+2)).eq.6) then
                     system(mb(i-2-1),md(j-2-1),ml(k+2)) = 1
                  else
                     write(*,*)'error in move 26 - 3 - 1'
                  end if
                  go to 333
               else
                  if(system(mb(i-1),md(j-4-1),ml(k+2)).eq.6) then
                     system(mb(i-1),md(j-4-1),ml(k+2)) = 1
                  else if(system(mb(i-2-1),md(j+2-1),ml(k+2)).eq.6) then
                     system(mb(i-2-1),md(j+2-1),ml(k+2)) = 1
                  else if(system(mb(i+2-1),md(j+2-1),ml(k+2)).eq.6) then
                     system(mb(i+2-1),md(j+2-1),ml(k+2)) = 1
                  else
                     write(*,*)'error in move 26 - 3 - 2'
                  end if
                  go to 333
               end if



 512           if(system(mb(i-1),md(j+2-1),ml(k+1)).eq.6) then
                  system(mb(i-1),md(j+2-1),ml(k+1)) = 1
                  if(system(mb(i-1),md(j+2-1),ml(k)).eq.5) then
                     system(mb(i-1),md(j+2-1),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1-1),md(j+1-1),ml(k+1)).eq.6) then
                  system(mb(i+1-1),md(j+1-1),ml(k+1)) = 1
                  if(system(mb(i+1-1),md(j+1-1),ml(k)).eq.5) then
                     system(mb(i+1-1),md(j+1-1),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1-1),md(j-1-1),ml(k+1)).eq.6) then
                  system(mb(i+1-1),md(j-1-1),ml(k+1)) = 1
                  if(system(mb(i+1-1),md(j-1-1),ml(k)).eq.5) then
                     system(mb(i+1-1),md(j-1-1),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j-2-1),ml(k+1)).eq.6) then
                  system(mb(i-1),md(j-2-1),ml(k+1)) = 1
                  if(system(mb(i-1),md(j-2-1),ml(k)).eq.5) then
                     system(mb(i-1),md(j-2-1),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1-1),md(j-1-1),ml(k+1)).eq.6) then
                  system(mb(i-1-1),md(j-1-1),ml(k+1)) = 1
                  if(system(mb(i-1-1),md(j-1-1),ml(k)).eq.5) then
                     system(mb(i-1-1),md(j-1-1),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1-1),md(j+1-1),ml(k+1)).eq.6) then
                  system(mb(i-1-1),md(j+1-1),ml(k+1)) = 1
                  if(system(mb(i-1-1),md(j+1-1),ml(k)).eq.5) then
                     system(mb(i-1-1),md(j+1-1),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 333
               end if

               if(up) then
                  if(system(mb(i-1),md(j+4-1),ml(k)).eq.5) then
                     system(mb(i-1),md(j+4-1),ml(k)) = 1
                  else if(system(mb(i+2-1),md(j-2-1),ml(k)).eq.5) then
                     system(mb(i+2-1),md(j-2-1),ml(k)) = 1
                  else if(system(mb(i-2-1),md(j-2-1),ml(k)).eq.5) then
                     system(mb(i-2-1),md(j-2-1),ml(k)) = 1
                  else
                     write(*,*)'error in move 26 - 4 - 1'
                  end if
                  go to 333
               else
                  if(system(mb(i-1),md(j-4-1),ml(k)).eq.5) then
                     system(mb(i-1),md(j-4-1),ml(k)) = 1
                  else if(system(mb(i-2-1),md(j+2-1),ml(k)).eq.5) then
                     system(mb(i-2-1),md(j+2-1),ml(k)) = 1
                  else if(system(mb(i+2-1),md(j+2-1),ml(k)).eq.5) then
                     system(mb(i+2-1),md(j+2-1),ml(k)) = 1
                  else
                     write(*,*)'error in move 26 - 4 - 2'
                  end if
                  go to 333
               end if               

c..............intimate dir 6
            else if(system(mb(i-1),md(j+1),k).eq.2) then
               system(mb(i-1),md(j+1),k) = 1

c - - - - - - -close interlayer divac

               if(system(mb(i-1),md(j+2+1),k).eq.5) then
                  system(mb(i-1),md(j+2+1),k) = 1
                  if(system(mb(i-1),md(j+2+1),ml(k+1)).eq.6) then
                     system(mb(i-1),md(j+2+1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1-1),md(j+1+1),k).eq.5) then
                  system(mb(i+1-1),md(j+1+1),k) = 1
                  if(system(mb(i+1-1),md(j+1+1),ml(k+1)).eq.6) then
                     system(mb(i+1-1),md(j+1+1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1-1),md(j-1+1),k).eq.5) then
                  system(mb(i+1-1),md(j-1+1),k) = 1
                  if(system(mb(i+1-1),md(j-1+1),ml(k+1)).eq.6) then
                     system(mb(i+1-1),md(j-1+1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j-2+1),k).eq.5) then
                  system(mb(i-1),md(j-2+1),k) = 1
                  if(system(mb(i-1),md(j-2+1),ml(k+1)).eq.6) then
                     system(mb(i-1),md(j-2+1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1-1),md(j-1+1),k).eq.5) then
                  system(mb(i-1-1),md(j-1+1),k) = 1
                  if(system(mb(i-1-1),md(j-1+1),ml(k+1)).eq.6) then
                     system(mb(i-1-1),md(j-1+1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1-1),md(j+1+1),k).eq.5) then
                  system(mb(i-1-1),md(j+1+1),k) = 1
                  if(system(mb(i-1-1),md(j+1+1),ml(k+1)).eq.6) then
                     system(mb(i-1-1),md(j+1+1),ml(k+1)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 513
               end if

               if(up) then
                  if(system(mb(i-1),md(j+4+1),ml(k+1)).eq.6) then
                     system(mb(i-1),md(j+4+1),ml(k+1)) = 1
                  else if(system(mb(i+2-1),md(j-2+1),ml(k+1)).eq.6) then
                     system(mb(i+2-1),md(j-2+1),ml(k+1)) = 1
                  else if(system(mb(i-2-1),md(j-2+1),ml(k+1)).eq.6) then
                     system(mb(i-2-1),md(j-2+1),ml(k+1)) = 1
                  else
                     write(*,*)'error in move 26 - 1 - 1'
                  end if
                  go to 333
               else
                  if(system(mb(i-1),md(j-4+1),ml(k+1)).eq.6) then
                     system(mb(i-1),md(j-4+1),ml(k+1)) = 1
                  else if(system(mb(i-2-1),md(j+2+1),ml(k+1)).eq.6) then
                     system(mb(i-2-1),md(j+2+1),ml(k+1)) = 1
                  else if(system(mb(i+2-1),md(j+2+1),ml(k+1)).eq.6) then
                     system(mb(i+2-1),md(j+2+1),ml(k+1)) = 1
                  else
                     write(*,*)'error in move 26 - 1 - 2'
                  end if
                  go to 333
               end if



 513           if(system(mb(i-1),md(j+2+1),k).eq.6) then
                  system(mb(i-1),md(j+2+1),k) = 1
                  if(system(mb(i-1),md(j+2+1),ml(k-1)).eq.5) then
                     system(mb(i-1),md(j+2+1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1-1),md(j+1+1),k).eq.6) then
                  system(mb(i+1-1),md(j+1+1),k) = 1
                  if(system(mb(i+1-1),md(j+1+1),ml(k-1)).eq.5) then
                     system(mb(i+1-1),md(j+1+1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1-1),md(j-1+1),k).eq.6) then
                  system(mb(i+1-1),md(j-1+1),k) = 1
                  if(system(mb(i+1-1),md(j-1+1),ml(k-1)).eq.5) then
                     system(mb(i+1-1),md(j-1+1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j-2+1),k).eq.6) then
                  system(mb(i-1),md(j-2+1),k) = 1
                  if(system(mb(i-1),md(j-2+1),ml(k-1)).eq.5) then
                     system(mb(i-1),md(j-2+1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1-1),md(j-1+1),k).eq.6) then
                  system(mb(i-1-1),md(j-1+1),k) = 1
                  if(system(mb(i-1-1),md(j-1+1),ml(k-1)).eq.5) then
                     system(mb(i-1-1),md(j-1+1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1-1),md(j+1+1),k).eq.6) then
                  system(mb(i-1-1),md(j+1+1),k) = 1
                  if(system(mb(i-1-1),md(j+1+1),ml(k-1)).eq.5) then
                     system(mb(i-1-1),md(j+1+1),ml(k-1)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 333
               end if

               if(up) then
                  if(system(mb(i-1),md(j+4+1),ml(k-1)).eq.5) then
                     system(mb(i-1),md(j+4+1),ml(k-1)) = 1
                  else if(system(mb(i+2-1),md(j-2+1),ml(k-1)).eq.5) then
                     system(mb(i+2-1),md(j-2+1),ml(k-1)) = 1
                  else if(system(mb(i-2+1),md(j-2+1),ml(k-1)).eq.5) then
                     system(mb(i-2-1),md(j-2+1),ml(k-1)) = 1
                  else
                     write(*,*)'error in move 26 - 2 - 1'
                  end if
                  go to 333
               else
                  if(system(mb(i-1),md(j-4+1),ml(k-1)).eq.5) then
                     system(mb(i-1),md(j-4+1),ml(k-1)) = 1
                  else if(system(mb(i-2-1),md(j+2+1),ml(k-1)).eq.5) then
                     system(mb(i-2-1),md(j+2+1),ml(k-1)) = 1
                  else if(system(mb(i+2-1),md(j+2+1),ml(k-1)).eq.5) then
                     system(mb(i+2-1),md(j+2+1),ml(k-1)) = 1
                  else
                     write(*,*)'error in move 26 - 2 - 2'
                  end if
                  go to 333
               end if               

            else if(system(mb(i-1),md(j+1),ml(k+1)).eq.2) then
               system(mb(i-1),md(j+1),ml(k+1)) = 1

c - - - - - - -close interlayer divac

               if(system(mb(i-1),md(j+2+1),ml(k+1)).eq.5) then
                  system(mb(i-1),md(j+2+1),ml(k+1)) = 1
                  if(system(mb(i-1),md(j+2+1),ml(k+2)).eq.6) then
                     system(mb(i-1),md(j+2+1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1-1),md(j+1+1),ml(k+1)).eq.5) then
                  system(mb(i+1-1),md(j+1+1),ml(k+1)) = 1
                  if(system(mb(i+1-1),md(j+1+1),ml(k+2)).eq.6) then
                     system(mb(i+1-1),md(j+1+1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1-1),md(j-1+1),ml(k+1)).eq.5) then
                  system(mb(i+1-1),md(j-1+1),ml(k+1)) = 1
                  if(system(mb(i+1-1),md(j-1+1),ml(k+2)).eq.6) then
                     system(mb(i+1-1),md(j-1+1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j-2+1),ml(k+1)).eq.5) then
                  system(mb(i-1),md(j-2+1),ml(k+1)) = 1
                  if(system(mb(i-1),md(j-2+1),ml(k+2)).eq.6) then
                     system(mb(i-1),md(j-2+1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1-1),md(j-1+1),ml(k+1)).eq.5) then
                  system(mb(i-1-1),md(j-1+1),ml(k+1)) = 1
                  if(system(mb(i-1-1),md(j-1+1),ml(k+2)).eq.6) then
                     system(mb(i-1-1),md(j-1+1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1-1),md(j+1+1),ml(k+1)).eq.5) then
                  system(mb(i-1-1),md(j+1+1),ml(k+1)) = 1
                  if(system(mb(i-1-1),md(j+1+1),ml(k+2)).eq.6) then
                     system(mb(i-1-1),md(j+1+1),ml(k+2)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 514
               end if

               if(up) then
                  if(system(mb(i-1),md(j+4+1),ml(k+2)).eq.6) then
                     system(mb(i-1),md(j+4+1),ml(k+2)) = 1
                  else if(system(mb(i+2-1),md(j-2+1),ml(k+2)).eq.6) then
                     system(mb(i+2-1),md(j-2+1),ml(k+2)) = 1
                  else if(system(mb(i-2+1),md(j-2+1),ml(k+2)).eq.6) then
                     system(mb(i-2-1),md(j-2+1),ml(k+2)) = 1
                  else
                     write(*,*)'error in move 26 - 3 - 1'
                  end if
                  go to 333
               else
                  if(system(mb(i-1),md(j-4+1),ml(k+2)).eq.6) then
                     system(mb(i-1),md(j-4+1),ml(k+2)) = 1
                  else if(system(mb(i-2-1),md(j+2+1),ml(k+2)).eq.6) then
                     system(mb(i-2-1),md(j+2+1),ml(k+2)) = 1
                  else if(system(mb(i+2-1),md(j+2+1),ml(k+2)).eq.6) then
                     system(mb(i+2-1),md(j+2+1),ml(k+2)) = 1
                  else
                     write(*,*)'error in move 26 - 3 - 2'
                  end if
                  go to 333
               end if



 514           if(system(mb(i-1),md(j+2+1),ml(k+1)).eq.6) then
                  system(mb(i-1),md(j+2+1),ml(k+1)) = 1
                  if(system(mb(i-1),md(j+2+1),ml(k)).eq.5) then
                     system(mb(i-1),md(j+2+1),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i+1-1),md(j+1+1),ml(k+1)).eq.6) then
                  system(mb(i+1-1),md(j+1+1),ml(k+1)) = 1
                  if(system(mb(i+1-1),md(j+1+1),ml(k)).eq.5) then
                     system(mb(i+1-1),md(j+1+1),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i+1-1),md(j-1+1),ml(k+1)).eq.6) then
                  system(mb(i+1-1),md(j-1+1),ml(k+1)) = 1
                  if(system(mb(i+1-1),md(j-1+1),ml(k)).eq.5) then
                     system(mb(i+1-1),md(j-1+1),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1),md(j-2+1),ml(k+1)).eq.6) then
                  system(mb(i-1),md(j-2+1),ml(k+1)) = 1
                  if(system(mb(i-1),md(j-2+1),ml(k)).eq.5) then
                     system(mb(i-1),md(j-2+1),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else if(system(mb(i-1-1),md(j-1+1),ml(k+1)).eq.6) then
                  system(mb(i-1-1),md(j-1+1),ml(k+1)) = 1
                  if(system(mb(i-1-1),md(j-1+1),ml(k)).eq.5) then
                     system(mb(i-1-1),md(j-1+1),ml(k)) = 1
                     go to 333
                  end if
                  up = .true.
               else if(system(mb(i-1-1),md(j+1+1),ml(k+1)).eq.6) then
                  system(mb(i-1-1),md(j+1+1),ml(k+1)) = 1
                  if(system(mb(i-1-1),md(j+1+1),ml(k)).eq.5) then
                     system(mb(i-1-1),md(j+1+1),ml(k)) = 1
                     go to 333
                  end if
                  up = .false.
               else
                  go to 333
               end if

               if(up) then
                  if(system(mb(i-1),md(j+4+1),ml(k)).eq.5) then
                     system(mb(i-1),md(j+4+1),ml(k)) = 1
                  else if(system(mb(i+2-1),md(j-2+1),ml(k)).eq.5) then
                     system(mb(i+2-1),md(j-2+1),ml(k)) = 1
                  else if(system(mb(i-2-1),md(j-2+1),ml(k)).eq.5) then
                     system(mb(i-2-1),md(j-2+1),ml(k)) = 1
                  else
                     write(*,*)'error in move 26 - 4 - 1'
                  end if
                  go to 333
               else
                  if(system(mb(i-1),md(j-4+1),ml(k)).eq.5) then
                     system(mb(i-1),md(j-4+1),ml(k)) = 1
                  else if(system(mb(i-2-1),md(j+2+1),ml(k)).eq.5) then
                     system(mb(i-2-1),md(j+2+1),ml(k)) = 1
                  else if(system(mb(i+2-1),md(j+2+1),ml(k)).eq.5) then
                     system(mb(i+2-1),md(j+2+1),ml(k)) = 1
                  else
                     write(*,*)'error in move 26 - 4 - 2'
                  end if
                  go to 333
               end if               

            end if
 333        isystem(i,j,k) = 8
         end if

c******lower grafted
      else if(isystem(i,j,k).ge.11 .and.isystem(i,j,k).le.16) then
         
         if(move(4).eq.21) then
            isystem(i,j,k) = isystem(i,j,k) + 2
            if(isystem(i,j,k).gt.16) then
               isystem(i,j,k) = isystem(i,j,k) - 6
            end if
         else if(move(4).eq.22) then
            isystem(i,j,k) = isystem(i,j,k) - 2
            if(isystem(i,j,k).lt.11) then
               isystem(i,j,k) = isystem(i,j,k) + 6
            end if
         else if(move(4).eq.23) then

            if(isystem(i,j,k).eq.11 .and. system(i,md(j+2),k).ne.0) then
               if(isystem(mb(i+1),md(j+3),k).eq.8) then
                  isystem(mb(i+1),md(j+3),k) = 15
               else
                  write(*,*)'error in itrans grafted 23 1'
                  stop
               end if
            else if(isystem(i,j,k).eq.13 .and. 
     &              system(i,md(j+2),k).ne.0) then
               if(isystem(mb(i+1),md(j-3),k).eq.8) then
                  isystem(mb(i+1),md(j-3),k) = 11
               else
                  write(*,*)'error in itrans grafted 23 3'
                  stop
               end if
            else if(isystem(i,j,k).eq.15 .and. 
     &              system(i,md(j+2),k).ne.0) then
               if(isystem(mb(i-2),j,k).eq.8) then
                  isystem(mb(i-2),j,k) = 13
               else
                  write(*,*)'error in itrans grafted 23 5'
                  stop
               end if
               
            else if(isystem(i,j,k).eq.12 .and. 
     &              system(i,md(j+2),k).eq.0) then
               if(isystem(mb(i+2),j,k).eq.8) then
                  isystem(mb(i+2),j,k) = 16
               else
                  write(*,*)'error in itrans grafted 23 2'
                  stop
               end if
            else if(isystem(i,j,k).eq.14 .and. 
     &              system(i,md(j+2),k).eq.0) then
               if(isystem(mb(i-1),md(j-3),k).eq.8) then
                  isystem(mb(i-1),md(j-3),k) = 12
               else
                  write(*,*)'error in itrans grafted 23 4'
                  stop
               end if
            else if(isystem(i,j,k).eq.16 .and. 
     &              system(i,md(j+2),k).eq.0) then
               if(isystem(mb(i-1),md(j+3),k).eq.8) then
                  isystem(mb(i-1),md(j+3),k) = 14
               else
                  write(*,*)'error in itrans grafted 23 6'
                  stop
               end if
            else
               write(*,*)'error in itrans grafted down 23'
               stop
            end if
            isystem(i,j,k) = 8
         else if(move(4).eq.24) then
            isystem(i,j,k) = isystem(i,j,k) - 10
         else if(move(4).eq.25) then
            if(system(i,j,k).eq.1) then
               system(i,j,k) = 9
            else
               write(*,*)'error in itrans grafted 25 down'
               stop
            end if
            isystem(i,j,k) = 8
         else if(move(4).eq.26) then
c........if up
            if(isystem(i,j,k).eq.11 .and. 
     &           system(i,md(j+2),k).ne.0) then
               system(i,md(j+2),k) = 9
            else if(isystem(i,j,k).eq.13 .and. 
     &           system(i,md(j+2),k).ne.0) then
               system(mb(i+1),md(j-1),k) = 9
            else if(isystem(i,j,k).eq.15 .and. 
     &              system(i,md(j+2),k).ne.0) then
               system(mb(i-1),md(j-1),k) = 9
            else if(isystem(i,j,k).eq.12 .and. 
     &              system(i,md(j+2),k).eq.0) then
               system(mb(i+1),md(j+1),k) = 9
            else if(isystem(i,j,k).eq.14 .and. 
     &              system(i,md(j+2),k).eq.0) then
               system(i,md(j-2),k) = 9
            else if(isystem(i,j,k).eq.16 .and. 
     &              system(i,md(j+2),k).eq.0) then
               system(mb(i-1),md(j+1),k) = 9
            else
               write(*,*)'error in itrans grafted 26 down'
               stop
            end if   
            isystem(i,j,k) = 8
         end if

c******upper grafted
      else if(isystem(i,j,k).ge.21 .and.isystem(i,j,k).le.26) then
         
         if(move(4).eq.21) then
            isystem(i,j,k) = isystem(i,j,k) + 2
            if(isystem(i,j,k).gt.26) then
               isystem(i,j,k) = isystem(i,j,k) - 6
            end if
         else if(move(4).eq.22) then
            isystem(i,j,k) = isystem(i,j,k) - 2
            if(isystem(i,j,k).lt.21) then
               isystem(i,j,k) = isystem(i,j,k) + 6
            end if
         else if(move(4).eq.23) then

            if(isystem(i,j,k).eq.21 .and. 
     &           system(i,md(j+2),ml(k+1)).ne.0) then
               if(isystem(mb(i+1),md(j+3),k).eq.8) then
                  isystem(mb(i+1),md(j+3),k) = 25
               else
                  write(*,*)'error in itrans grafted 23 1 u'
                  stop
               end if
            else if(isystem(i,j,k).eq.23 .and. 
     &              system(i,md(j+2),ml(k+1)).ne.0) then
               if(isystem(mb(i+1),md(j-3),k).eq.8) then
                  isystem(mb(i+1),md(j-3),k) = 21
               else
                  write(*,*)'error in itrans grafted 23 3 u'
                  stop
               end if
            else if(isystem(i,j,k).eq.25 .and. 
     &              system(i,md(j+2),ml(k+1)).ne.0) then
               if(isystem(mb(i-2),j,k).eq.8) then
                  isystem(mb(i-2),j,k) = 23
               else
                  write(*,*)'error in itrans grafted 23 5 u'
                  stop
               end if
               
            else if(isystem(i,j,k).eq.22 .and. 
     &              system(i,md(j+2),ml(k+1)).eq.0) then
               if(isystem(mb(i+2),j,k).eq.8) then
                  isystem(mb(i+2),j,k) = 26
               else
                  write(*,*)'error in itrans grafted 23 2 u'
                  stop
               end if
            else if(isystem(i,j,k).eq.24 .and. 
     &              system(i,md(j+2),ml(k+1)).eq.0) then
               if(isystem(mb(i-1),md(j-3),k).eq.8) then
                  isystem(mb(i-1),md(j-3),k) = 22
               else
                  write(*,*)isystem(mb(i-1),md(j-3),k),i-1,j-3,k
                  write(*,*)'error in itrans grafted 23 4 u'
                  stop
               end if
            else if(isystem(i,j,k).eq.26 .and. 
     &              system(i,md(j+2),ml(k+1)).eq.0) then
               if(isystem(mb(i-1),md(j+3),k).eq.8) then
                  isystem(mb(i-1),md(j+3),k) = 24
               else
                  write(*,*)'error in itrans grafted 23 6 u'
                  stop
               end if
            else
               write(*,*)'error in itrans grafted up 23'
               stop
            end if
            isystem(i,j,k) = 8
         else if(move(4).eq.24) then
            isystem(i,j,k) = isystem(i,j,k) - 20
         else if(move(4).eq.25) then
            if(system(i,j,ml(k+1)).eq.1) then
               system(i,j,ml(k+1)) = 9
            else
               write(*,*)'error in itrans grafted 25 up'
               stop
            end if
            isystem(i,j,k) = 8
         else if(move(4).eq.26) then
c........if up
            if(isystem(i,j,k).eq.21 .and. 
     &           system(i,md(j+2),ml(k+1)).ne.0) then
               system(i,md(j+2),ml(k+1)) = 9
            else if(isystem(i,j,k).eq.23 .and. 
     &           system(i,md(j+2),ml(k+1)).ne.0) then
               system(mb(i+1),md(j-1),ml(k+1)) = 9
            else if(isystem(i,j,k).eq.25 .and. 
     &              system(i,md(j+2),ml(k+1)).ne.0) then
               system(mb(i-1),md(j-1),ml(k+1)) = 9
            else if(isystem(i,j,k).eq.22 .and. 
     &              system(i,md(j+2),ml(k+1)).eq.0) then
               system(mb(i+1),md(j+1),ml(k+1)) = 9
            else if(isystem(i,j,k).eq.24 .and. 
     &              system(i,md(j+2),ml(k+1)).eq.0) then
               system(i,md(j-2),ml(k+1)) = 9
            else if(isystem(i,j,k).eq.26 .and. 
     &              system(i,md(j+2),ml(k+1)).eq.0) then
               system(mb(i-1),md(j+1),ml(k+1)) = 9
            else
               write(*,*)'error in itrans grafted 26 down'
               stop
            end if   
            isystem(i,j,k) = 8
         end if

      end if

      end 


      subroutine maketrans(system,move)
c---------------------------------------------------------
c     move atoms acording to the site and kmc direction  -  
c     move(1-3) : site indices                           -
c     move(4) : direction                                -
c             : 1-6 in plane, 7 up, 8 down               -
c             : 9 move atom to upper plane               -
c             : 10 move atom to lower plane              -
c---------------------------------------------------------   
      integer*1 system(800,800,199)
      integer move(4)
      integer*2 atm1,atm2,atm3

      i = move(1)
      j = move(2)
      k = move(3)

      if(move(4).eq.1) then
c........do dimer slide to avoid undercoordinated atom
         if(system(mb(i+2),md(j),k).eq.2 .or.  
     &        system(mb(i+1),md(j-3),k).eq.2) then
            system(mb(i),md(j+2),k) = 1
            system(mb(i+1),md(j-1),k) = 2
         else if(system(mb(i-2),md(j),k).eq.2 .or.  
     &        system(mb(i-1),md(j-3),k).eq.2) then
            system(mb(i),md(j+2),k) = 1
            system(mb(i-1),md(j-1),k) = 2
         else
            system(mb(i),md(j+2),k) = 1
            system(mb(i),md(j),k) = 2
         end if
      else if(move(4).eq.2) then
         if(system(mb(i+1),md(j-3),k).eq.2 .or.  
     &        system(mb(i-1),md(j-3),k).eq.2) then
            system(mb(i+1),md(j+1),k) = 1
            system(mb(i),md(j-2),k) = 2
         else if(system(mb(i-1),md(j+3),k).eq.2 .or.  
     &           system(mb(i-2),md(j),k).eq.2) then
            system(mb(i+1),md(j+1),k) = 1
            system(mb(i-1),md(j+1),k) = 2
         else
            system(mb(i+1),md(j+1),k) = 1
            system(mb(i),md(j),k) = 2
         end if
      else if(move(4).eq.3) then
         if(system(mb(i-1),md(j-3),k).eq.2 .or.  
     &        system(mb(i-2),md(j),k).eq.2) then
            system(mb(i+1),md(j-1),k) = 1
            system(mb(i-1),md(j-1),k) = 2
         else if(system(mb(i+1),md(j+3),k).eq.2 .or.  
     &           system(mb(i-1),md(j+3),k).eq.2) then
            system(mb(i+1),md(j-1),k) = 1
            system(mb(i),md(j+2),k) = 2
         else
            system(mb(i+1),md(j-1),k) = 1
            system(mb(i),md(j),k) = 2
         end if

      else if(move(4).eq.4) then
c........do dimer slide to avoid undercoordinated atom
         if(system(mb(i-2),md(j),k).eq.2 .or.  
     &        system(mb(i-1),md(j+3),k).eq.2) then
            system(mb(i),md(j-2),k) = 1
            system(mb(i-1),md(j+1),k) = 2
         else if(system(mb(i+2),md(j),k).eq.2 .or.  
     &        system(mb(i+1),md(j+3),k).eq.2) then
            system(mb(i),md(j-2),k) = 1
            system(mb(i+1),md(j+1),k) = 2
         else
            system(mb(i),md(j-2),k) = 1
            system(mb(i),md(j),k) = 2
         end if
      else if(move(4).eq.5) then
         if(system(mb(i-1),md(j+3),k).eq.2 .or.  
     &        system(mb(i+1),md(j+3),k).eq.2) then
            system(mb(i-1),md(j-1),k) = 1
            system(mb(i),md(j+2),k) = 2
         else if(system(mb(i+1),md(j-3),k).eq.2 .or.  
     &           system(mb(i+2),md(j),k).eq.2) then
            system(mb(i-1),md(j-1),k) = 1
            system(mb(i+1),md(j-1),k) = 2
         else
            system(mb(i-1),md(j-1),k) = 1
            system(mb(i),md(j),k) = 2
         end if
      else if(move(4).eq.6) then
         if(system(mb(i+1),md(j+3),k).eq.2 .or.  
     &        system(mb(i+2),md(j),k).eq.2) then
            system(mb(i-1),md(j+1),k) = 1
            system(mb(i+1),md(j+1),k) = 2
         else if(system(mb(i-1),md(j-3),k).eq.2 .or.  
     &           system(mb(i+1),md(j-3),k).eq.2) then
            system(mb(i-1),md(j+1),k) = 1
            system(mb(i),md(j-2),k) = 2
         else
            system(mb(i-1),md(j+1),k) = 1
            system(mb(i),md(j),k) = 2
         end if

c-----making and breaking interlayer bonds


c.....move up
      else if(move(4).eq.7) then
c........dv1 and dv2
         if(system(i,j,k).eq.1 .and. system(i,j,ml(k+1)).eq.1) then
            system(i,j,k) = 5
            system(i,j,ml(k+1)) = 6
         else if(system(i,j,k).eq.6 .and. system(i,j,ml(k-1)).eq.5)then
            system(i,j,k) = 1
            system(i,j,ml(k-1)) = 1
c........dv3 direction 1
         else if(system(i,j,k).eq.1 .and. system(i,md(j+2),k).eq.2) then
            if(system(i,j,ml(k+1)).eq.2) then
               system(i,j,k) = 5
               system(i,md(j-2),ml(k+1)) = 6
            else if(system(i,md(j+4),ml(k+1)).eq.2) then
               system(i,j,k) = 5
               system(i,md(j+2),ml(k+1)) = 6
            else
               write(*,*)'Error in move 7 - 1'
               stop
            endif
c........dv3 direction 2
         else if(system(i,j,k).eq.1 .and. 
     &           system(mb(i+1),md(j+1),k).eq.2)then
            if(system(i,j,ml(k+1)).eq.2) then
               system(i,j,k) = 5
               system(mb(i-1),md(j-1),ml(k+1)) = 6
            else if(system(mb(i+2),md(j+2),ml(k+1)).eq.2) then
               system(i,j,k) = 5
               system(mb(i+1),md(j+1),ml(k+1)) = 6
            else
               write(*,*)'Error in move 7 - 2'
               stop
            endif
c........dv3 direction 3
         else if(system(i,j,k).eq.1 .and. 
     &           system(mb(i+1),md(j-1),k).eq.2)then
            if(system(i,j,ml(k+1)).eq.2) then
               system(i,j,k) = 5
               system(mb(i-1),md(j+1),ml(k+1)) = 6
            else if(system(mb(i+2),md(j-2),ml(k+1)).eq.2) then
               system(i,j,k) = 5
               system(mb(i+1),md(j-1),ml(k+1)) = 6
            else
               write(*,*)'Error in move 7 - 3'
               stop
            endif
c........dv3 direction 4
         else if(system(i,j,k).eq.1 .and. system(i,md(j-2),k).eq.2) then
            if(system(i,j,ml(k+1)).eq.2) then
               system(i,j,k) = 5
               system(i,md(j+2),ml(k+1)) = 6
            else if(system(i,md(j-4),ml(k+1)).eq.2) then
               system(i,j,k) = 5
               system(i,md(j-2),ml(k+1)) = 6
            else
               write(*,*)'Error in move 7 - 4'
               stop
            endif
c........dv3 direction 5
         else if(system(i,j,k).eq.1 .and. 
     &           system(mb(i-1),md(j-1),k).eq.2)then
            if(system(i,j,ml(k+1)).eq.2) then
               system(i,j,k) = 5
               system(mb(i+1),md(j+1),ml(k+1)) = 6
            else if(system(mb(i-2),md(j-2),ml(k+1)).eq.2) then
               system(i,j,k) = 5
               system(mb(i-1),md(j-1),ml(k+1)) = 6
            else
               write(*,*)'Error in move 7 - 5'
               stop
            endif
c........dv3 direction 6
         else if(system(i,j,k).eq.1 .and. 
     &           system(mb(i-1),md(j+1),k).eq.2)then
            if(system(i,j,ml(k+1)).eq.2) then
               system(i,j,k) = 5
               system(mb(i+1),md(j-1),ml(k+1)) = 6
            else if(system(mb(i-2),md(j+2),ml(k+1)).eq.2) then
               system(i,j,k) = 5
               system(mb(i-1),md(j+1),ml(k+1)) = 6
            else
               write(*,*)'Error in move 7 - 6'
               stop
            endif

c--------breaking bond
c........dv3 direction 1
         else if(system(i,j,k).eq.6 .and. system(i,md(j+2),k).eq.2) then
            if(system(i,j,ml(k-1)).eq.2) then
               system(i,j,k) = 1
               system(i,md(j-2),ml(k-1)) = 1
            else if(system(i,md(j+4),ml(k-1)).eq.2) then
               system(i,j,k) = 1
               system(i,md(j+2),ml(k-1)) = 1
            else
               write(*,*)'Error in move 7 - 1 - 2'
               stop
            endif
c........dv3 direction 2
         else if(system(i,j,k).eq.6 .and. 
     &           system(mb(i+1),md(j+1),k).eq.2)then
            if(system(i,j,ml(k-1)).eq.2) then
               system(i,j,k) = 1
               system(mb(i-1),md(j-1),ml(k-1)) = 1
            else if(system(mb(i+2),md(j+2),ml(k-1)).eq.2) then
               system(i,j,k) = 1
               system(mb(i+1),md(j+1),ml(k-1)) = 1
            else
               write(*,*)'Error in move 7 - 2 - 2'
               stop
            endif
c........dv3 direction 3
         else if(system(i,j,k).eq.6 .and. 
     &           system(mb(i+1),md(j-1),k).eq.2)then
            if(system(i,j,ml(k-1)).eq.2) then
               system(i,j,k) = 1
               system(mb(i-1),md(j+1),ml(k-1)) = 1
            else if(system(mb(i+2),md(j-2),ml(k-1)).eq.2) then
               system(i,j,k) = 1
               system(mb(i+1),md(j-1),ml(k-1)) = 1
            else
               write(*,*)'Error in move 7 - 3 - 2'
               stop
            endif
c........dv3 direction 4
         else if(system(i,j,k).eq.6 .and. system(i,md(j-2),k).eq.2) then
            if(system(i,j,ml(k-1)).eq.2) then
               system(i,j,k) = 1
               system(i,md(j+2),ml(k-1)) = 1
            else if(system(i,md(j-4),ml(k-1)).eq.2) then
               system(i,j,k) = 1
               system(i,md(j-2),ml(k-1)) = 1
            else
               write(*,*)'Error in move 7 - 4 - 2'
               stop
            endif
c........dv3 direction 5
         else if(system(i,j,k).eq.6 .and. 
     &           system(mb(i-1),md(j-1),k).eq.2)then
            if(system(i,j,ml(k-1)).eq.2) then
               system(i,j,k) = 1
               system(mb(i+1),md(j+1),ml(k-1)) = 1
            else if(system(mb(i-2),md(j-2),ml(k-1)).eq.2) then
               system(i,j,k) = 1
               system(mb(i-1),md(j-1),ml(k-1)) = 1
            else
               write(*,*)'Error in move 7 - 5 - 2'
               stop
            endif
c........dv3 direction 6
         else if(system(i,j,k).eq.6 .and. 
     &           system(mb(i-1),md(j+1),k).eq.2)then
            if(system(i,j,ml(k-1)).eq.2) then
               system(i,j,k) = 1
               system(mb(i+1),md(j-1),ml(k-1)) = 1
            else if(system(mb(i-2),md(j+2),ml(k-1)).eq.2) then
               system(i,j,k) = 1
               system(mb(i-1),md(j+1),ml(k-1)) = 1
            else
               write(*,*)'Error in move 7 - 6 - 2'
               stop
            endif
         else
            write(*,*)'Error in move 7'
            stop
         endif

c.....move down


      else if(move(4).eq.8) then
         if(system(i,j,k).eq.1 .and. system(i,j,ml(k-1)).eq.1) then
            system(i,j,k) = 6
            system(i,j,ml(k-1)) = 5
         else if(system(i,j,k).eq.5 .and. system(i,j,ml(k+1)).eq.6) then
            system(i,j,k) = 1
            system(i,j,ml(k+1)) = 1
c........dv3 direction 1
         else if(system(i,j,k).eq.1 .and. system(i,md(j+2),k).eq.2) then
            if(system(i,j,ml(k-1)).eq.2) then
               system(i,j,k) = 6
               system(i,md(j-2),ml(k-1)) = 5
            else if(system(i,md(j+4),ml(k-1)).eq.2) then
               system(i,j,k) = 6
               system(i,md(j+2),ml(k-1)) = 5
            else
               write(*,*)'Error in move 8 - 1'
               stop
            endif
c........dv3 direction 2
         else if(system(i,j,k).eq.1 .and. 
     &           system(mb(i+1),md(j+1),k).eq.2)then
            if(system(i,j,ml(k-1)).eq.2) then
               system(i,j,k) = 6
               system(mb(i-1),md(j-1),ml(k-1)) = 5
            else if(system(mb(i+2),md(j+2),ml(k-1)).eq.2) then
               system(i,j,k) = 6
               system(mb(i+1),md(j+1),ml(k-1)) = 5
            else
               write(*,*)'Error in move 8 - 2'
               stop
            endif
c........dv3 direction 3
         else if(system(i,j,k).eq.1 .and. 
     &           system(mb(i+1),md(j-1),k).eq.2)then
            if(system(i,j,ml(k-1)).eq.2) then
               system(i,j,k) = 6
               system(mb(i-1),md(j+1),ml(k-1)) = 5
            else if(system(mb(i+2),md(j-2),ml(k-1)).eq.2) then
               system(i,j,k) = 6
               system(mb(i+1),md(j-1),ml(k-1)) = 5
            else
               write(*,*)'Error in move 8 - 3'
               stop
            endif
c........dv3 direction 4
         else if(system(i,j,k).eq.1 .and. system(i,md(j-2),k).eq.2) then
            if(system(i,j,ml(k-1)).eq.2) then
               system(i,j,k) = 6
               system(i,md(j+2),ml(k-1)) = 5
            else if(system(i,md(j-4),ml(k-1)).eq.2) then
               system(i,j,k) = 6
               system(i,md(j-2),ml(k-1)) = 5
            else
               write(*,*)'Error in move 8 - 4'
               stop
            endif
c........dv3 direction 5
         else if(system(i,j,k).eq.1 .and. 
     &           system(mb(i-1),md(j-1),k).eq.2)then
            if(system(i,j,ml(k-1)).eq.2) then
               system(i,j,k) = 6
               system(mb(i+1),md(j+1),ml(k-1)) = 5
            else if(system(mb(i-2),md(j-2),ml(k-1)).eq.2) then
               system(i,j,k) = 6
               system(mb(i-1),md(j-1),ml(k-1)) = 5
            else
               write(*,*)'Error in move 8 - 5'
               stop
            endif
c........dv3 direction 6
         else if(system(i,j,k).eq.1 .and. 
     &           system(mb(i-1),md(j+1),k).eq.2)then
            if(system(i,j,ml(k-1)).eq.2) then
               system(i,j,k) = 6
               system(mb(i+1),md(j-1),ml(k-1)) = 5
            else if(system(mb(i-2),md(j+2),ml(k-1)).eq.2) then
               system(i,j,k) = 6
               system(mb(i-1),md(j+1),ml(k-1)) = 5
            else
               write(*,*)'Error in move 8 - 6'
               stop
            endif

c--------breaking bond
c........dv3 direction 1
         else if(system(i,j,k).eq.5 .and. system(i,md(j+2),k).eq.2) then
            if(system(i,j,ml(k+1)).eq.2) then
               system(i,j,k) = 1
               system(i,md(j-2),ml(k+1)) = 1
            else if(system(i,md(j+4),ml(k+1)).eq.2) then
               system(i,j,k) = 1
               system(i,md(j+2),ml(k+1)) = 1
            else
               write(*,*)'Error in move 8 - 1 - 2'
               stop
            endif
c........dv3 direction 2
         else if(system(i,j,k).eq.5 .and. 
     &           system(mb(i+1),md(j+1),k).eq.2)then
            if(system(i,j,ml(k+1)).eq.2) then
               system(i,j,k) = 1
               system(mb(i-1),md(j-1),ml(k+1)) = 1
            else if(system(mb(i+2),md(j+2),ml(k+1)).eq.2) then
               system(i,j,k) = 1
               system(mb(i+1),md(j+1),ml(k+1)) = 1
            else
               write(*,*)'Error in move 8 - 2 - 2'
               stop
            endif
c........dv3 direction 3
         else if(system(i,j,k).eq.5 .and. 
     &           system(mb(i+1),md(j-1),k).eq.2)then
            if(system(i,j,ml(k+1)).eq.2) then
               system(i,j,k) = 1
               system(mb(i-1),md(j+1),ml(k+1)) = 1
            else if(system(mb(i+2),md(j-2),ml(k+1)).eq.2) then
               system(i,j,k) = 1
               system(mb(i+1),md(j-1),ml(k+1)) = 1
            else
               write(*,*)'Error in move 8 - 3 - 2'
               stop
            endif
c........dv3 direction 4
         else if(system(i,j,k).eq.5 .and. system(i,md(j-2),k).eq.2) then
            if(system(i,j,ml(k+1)).eq.2) then
               system(i,j,k) = 1
               system(i,md(j+2),ml(k+1)) = 1
            else if(system(i,md(j-4),ml(k+1)).eq.2) then
               system(i,j,k) = 1
               system(i,md(j-2),ml(k+1)) = 1
            else
               write(*,*)'Error in move 8 - 4 - 2'
               stop
            endif
c........dv3 direction 5
         else if(system(i,j,k).eq.5 .and. 
     &           system(mb(i-1),md(j-1),k).eq.2)then
            if(system(i,j,ml(k+1)).eq.2) then
               system(i,j,k) = 1
               system(mb(i+1),md(j+1),ml(k+1)) = 1
            else if(system(mb(i-2),md(j-2),ml(k+1)).eq.2) then
               system(i,j,k) = 1
               system(mb(i-1),md(j-1),ml(k+1)) = 1
            else
               write(*,*)'Error in move 8 - 5 - 2'
               stop
            endif
c........dv3 direction 6
         else if(system(i,j,k).eq.5 .and. 
     &           system(mb(i-1),md(j+1),k).eq.2)then
            if(system(i,j,ml(k+1)).eq.2) then
               system(i,j,k) = 1
               system(mb(i+1),md(j-1),ml(k+1)) = 1
            else if(system(mb(i-2),md(j+2),ml(k+1)).eq.2) then
               system(i,j,k) = 1
               system(mb(i-1),md(j+1),ml(k+1)) = 1
            else
               write(*,*)'Error in move 8 - 6 - 2'
               stop
            endif
         else
            write(*,*)'Error in move 8'
            stop
         end if


c............................................
c.....move to the upper layer (atom transfer)
c............................................
      else if(move(4).eq.9) then
c--------DVIL1 and DVIL2
         if(system(i,j,k).eq.5 .and. system(i,j,ml(k+1)).eq.6) then
             system(i,j,ml(k+1)) = 1
c-----------DVIL1 direction 1
            if(system(i,md(j+2),k).eq.2 .and. 
     &           system(i,md(j-2),ml(k+1)).eq.2) then
               if(system(mb(i-1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k+1)) = 1
               else if(system(mb(i-1),md(j-3),k).eq.2 .and. 
     &                 system(mb(i+1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k+1)) = 1
                  system(mb(i-1),md(j-1),k) = 2
                  system(mb(i-1),md(j-3),k) = 1
               else if(system(mb(i-1),md(j-3),k).eq.1 .and. 
     &                 system(mb(i+1),md(j-3),k).eq.2) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k+1)) = 1
                  system(mb(i+1),md(j-1),k) = 2
                  system(mb(i+1),md(j-3),k) = 1
               else
                  write(*,*)'Error in jump up 1 - 1'
                  stop
               end if
c-----------DVIL2 direction 1 right 
            else if(system(i,md(j+2),k).eq.2 .and. 
     &              system(mb(i+1),md(j+1),ml(k+1)).eq.2) then
               if(system(mb(i-1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k+1)) = 1
               else if(system(mb(i-1),md(j-3),k).eq.2 .and. 
     &                 system(mb(i+1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k+1)) = 1
                  system(mb(i-1),md(j-1),k) = 2
                  system(mb(i-1),md(j-3),k) = 1
               else if(system(mb(i-1),md(j-3),k).eq.1 .and. 
     &                 system(mb(i+1),md(j-3),k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k+1)) = 1
                  system(mb(i+1),md(j-1),k) = 2
                  system(mb(i+1),md(j-3),k) = 1
               else
                  write(*,*)'Error in jump up 1 - 2'
                  stop
               end if
c-----------DVIL2 direction 1 left
            else if(system(i,md(j+2),k).eq.2 .and. 
     &              system(mb(i-1),md(j+1),ml(k+1)).eq.2) then
               if(system(mb(i-1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k+1)) = 1
               else if(system(mb(i-1),md(j-3),k).eq.2 .and. 
     &                 system(mb(i+1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k+1)) = 1
                  system(mb(i-1),md(j-1),k) = 2
                  system(mb(i-1),md(j-3),k) = 1
               else if(system(mb(i-1),md(j-3),k).eq.1 .and. 
     &                 system(mb(i+1),md(j-3),k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k+1)) = 1
                  system(mb(i+1),md(j-1),k) = 2
                  system(mb(i+1),md(j-3),k) = 1
               else
                  write(*,*)'Error in jump up 1 - 3'
                  stop
               end if
c-----------DVIL1 direction 2
            else if(system(mb(i+1),md(j+1),k).eq.2 .and. 
     &              system(mb(i-1),md(j-1),ml(k+1)).eq.2) then
               if(system(mb(i-2),j,k).eq.1 .and. 
     &              system(mb(i-1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k+1)) = 1
               else if(system(mb(i-2),j,k).eq.2 .and. 
     &              system(mb(i-1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k+1)) = 1
                  system(mb(i-1),md(j+1),k) = 2
                  system(mb(i-2),j,k) = 1
               else if(system(mb(i-2),j,k).eq.1 .and. 
     &              system(mb(i-1),md(j-3),k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k+1)) = 1
                  system(i,md(j-2),k) = 2
                  system(mb(i-1),md(j-3),k) = 1
               else
                  write(*,*)'Error in jump up 2 - 1'
                  stop
               end if
c-----------DVIL2 direction 2 right 
            else if(system(mb(i+1),md(j+1),k).eq.2 .and. 
     &              system(mb(i+1),md(j-1),ml(k+1)).eq.2) then
               if(system(mb(i-2),j,k).eq.1 .and. 
     &              system(mb(i-1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k+1)) = 1
               else if(system(mb(i-2),j,k).eq.2 .and. 
     &              system(mb(i-1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k+1)) = 1
                  system(mb(i-1),md(j+1),k) = 2
                  system(mb(i-2),j,k) = 1
               else if(system(mb(i-2),j,k).eq.1 .and. 
     &              system(mb(i-1),md(j-3),k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k+1)) = 1
                  system(i,md(j-2),k) = 2
                  system(mb(i-1),md(j-3),k) = 1
               else
                  write(*,*)'Error in jump up 2 - 2'
                  stop
               end if
c-----------DVIL2 direction 2 left
            else if(system(mb(i+1),md(j+1),k).eq.2 .and. 
     &              system(i,md(j+2),ml(k+1)).eq.2) then
               if(system(mb(i-2),j,k).eq.1 .and. 
     &              system(mb(i-1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k+1)) = 1
               else if(system(mb(i-2),j,k).eq.2 .and. 
     &              system(mb(i-1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k+1)) = 1
                  system(mb(i-1),md(j+1),k) = 2
                  system(mb(i-2),j,k) = 1
               else if(system(mb(i-2),j,k).eq.1 .and. 
     &              system(mb(i-1),md(j-3),k).eq.2) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k+1)) = 1
                  system(i,md(j-2),k) = 2
                  system(mb(i-1),md(j-3),k) = 1
               else
                  write(*,*)'Error in jump up 2 - 2'
                  stop
               end if

c-----------DVIL1 direction 3
            else if(system(mb(i+1),md(j-1),k).eq.2 .and. 
     &              system(mb(i-1),md(j+1),ml(k+1)).eq.2) then
               if(system(mb(i-1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k+1)) = 1
               else if(system(mb(i-1),md(j+3),k).eq.2 .and. 
     &              system(mb(i-2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k+1)) = 1
                  system(i,md(j+2),k) = 2
                  system(mb(i-1),md(j+3),k) = 1
               else if(system(mb(i-1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-2),j,k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k+1)) = 1
                  system(mb(i-1),md(j-1),k) = 2
                  system(mb(i-2),j,k) = 1
               else
                  write(*,*)'Error in jump up 3 - 1'
                  stop
               end if

c-----------DVIL2 direction 3 right
            else if(system(mb(i+1),md(j-1),k).eq.2 .and. 
     &              system(i,md(j-2),ml(k+1)).eq.2) then
               if(system(mb(i-1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k+1)) = 1
               else if(system(mb(i-1),md(j+3),k).eq.2 .and. 
     &              system(mb(i-2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k+1)) = 1
                  system(i,md(j+2),k) = 2
                  system(mb(i-1),md(j+3),k) = 1
               else if(system(mb(i-1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-2),j,k).eq.2) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k+1)) = 1
                  system(mb(i-1),md(j-1),k) = 2
                  system(mb(i-2),j,k) = 1
               else
                  write(*,*)'Error in jump up 3 - 2'
                  stop
               end if
c-----------DVIL2 direction 3 left
            else if(system(mb(i+1),md(j-1),k).eq.2 .and. 
     &              system(mb(i+1),md(j+1),ml(k+1)).eq.2) then
               if(system(mb(i-1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k+1)) = 1
               else if(system(mb(i-1),md(j+3),k).eq.2 .and. 
     &              system(mb(i-2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k+1)) = 1
                  system(i,md(j+2),k) = 2
                  system(mb(i-1),md(j+3),k) = 1
               else if(system(mb(i-1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-2),j,k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k+1)) = 1
                  system(mb(i-1),md(j-1),k) = 2
                  system(mb(i-2),j,k) = 1
               else
                  write(*,*)'Error in jump up 3 - 3'
                  stop
               end if

c+++++++++++DVIL1 direction 4
            else if(system(i,md(j-2),k).eq.2 .and. 
     &           system(i,md(j+2),ml(k+1)).eq.2) then
               if(system(mb(i+1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k+1)) = 1
               else if(system(mb(i+1),md(j+3),k).eq.2 .and. 
     &                 system(mb(i-1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k+1)) = 1
                  system(mb(i+1),md(j+1),k) = 2
                  system(mb(i+1),md(j+3),k) = 1
               else if(system(mb(i+1),md(j+3),k).eq.1 .and. 
     &                 system(mb(i-1),md(j+3),k).eq.2) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k+1)) = 1
                  system(mb(i-1),md(j+1),k) = 2
                  system(mb(i-1),md(j+3),k) = 1
               else
                  write(*,*)'Error in jump up 1 + 1'
                  stop
               end if
c+++++++++++DVIL2 direction 4 right 
            else if(system(i,md(j-2),k).eq.2 .and. 
     &              system(mb(i-1),md(j-1),ml(k+1)).eq.2) then
               if(system(mb(i+1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k+1)) = 1
               else if(system(mb(i+1),md(j+3),k).eq.2 .and. 
     &                 system(mb(i-1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k+1)) = 1
                  system(mb(i+1),md(j+1),k) = 2
                  system(mb(i+1),md(j+3),k) = 1
               else if(system(mb(i+1),md(j+3),k).eq.1 .and. 
     &                 system(mb(i-1),md(j+3),k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k+1)) = 1
                  system(mb(i-1),md(j+1),k) = 2
                  system(mb(i-1),md(j+3),k) = 1
               else
                  write(*,*)'Error in jump up 1 + 2'
                  stop
               end if
c+++++++++++DVIL2 direction 4 left
            else if(system(i,md(j-2),k).eq.2 .and. 
     &              system(mb(i+1),md(j-1),ml(k+1)).eq.2) then
               if(system(mb(i+1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k+1)) = 1
               else if(system(mb(i+1),md(j+3),k).eq.2 .and. 
     &                 system(mb(i-1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k+1)) = 1
                  system(mb(i+1),md(j+1),k) = 2
                  system(mb(i+1),md(j+3),k) = 1
               else if(system(mb(i+1),md(j+3),k).eq.1 .and. 
     &                 system(mb(i-1),md(j+3),k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k+1)) = 1
                  system(mb(i-1),md(j+1),k) = 2
                  system(mb(i-1),md(j+3),k) = 1
               else
                  write(*,*)'Error in jump up 1 + 3'
                  stop
               end if
c+++++++++++DVIL1 direction 5
            else if(system(mb(i-1),md(j-1),k).eq.2 .and. 
     &              system(mb(i+1),md(j+1),ml(k+1)).eq.2) then
               if(system(mb(i+2),j,k).eq.1 .and. 
     &              system(mb(i+1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k+1)) = 1
               else if(system(mb(i+2),j,k).eq.2 .and. 
     &              system(mb(i+1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k+1)) = 1
                  system(mb(i+1),md(j-1),k) = 2
                  system(mb(i+2),j,k) = 1
               else if(system(mb(i+2),j,k).eq.1 .and. 
     &              system(mb(i+1),md(j+3),k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k+1)) = 1
                  system(i,md(j+2),k) = 2
                  system(mb(i+1),md(j+3),k) = 1
               else
                  write(*,*)'Error in jump up 2 + 1'
                  stop
               end if
c+++++++++++DVIL2 direction 5 right 
            else if(system(mb(i-1),md(j-1),k).eq.2 .and. 
     &              system(mb(i-1),md(j+1),ml(k+1)).eq.2) then
               if(system(mb(i+2),j,k).eq.1 .and. 
     &              system(mb(i+1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k+1)) = 1
               else if(system(mb(i+2),j,k).eq.2 .and. 
     &              system(mb(i+1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k+1)) = 1
                  system(mb(i+1),md(j-1),k) = 2
                  system(mb(i+2),j,k) = 1
               else if(system(mb(i+2),j,k).eq.1 .and. 
     &              system(mb(i+1),md(j+3),k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k+1)) = 1
                  system(i,md(j+2),k) = 2
                  system(mb(i+1),md(j+3),k) = 1
               else
                  write(*,*)'Error in jump up 2 + 2'
                  stop
               end if
c+++++++++++DVIL2 direction 5 left
            else if(system(mb(i-1),md(j-1),k).eq.2 .and. 
     &              system(i,md(j-2),ml(k+1)).eq.2) then
               if(system(mb(i+2),j,k).eq.1 .and. 
     &              system(mb(i+1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k+1)) = 1
               else if(system(mb(i+2),j,k).eq.2 .and. 
     &              system(mb(i+1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k+1)) = 1
                  system(mb(i+1),md(j-1),k) = 2
                  system(mb(i+2),j,k) = 1
               else if(system(mb(i+2),j,k).eq.1 .and. 
     &              system(mb(i+1),md(j+3),k).eq.2) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k+1)) = 1
                  system(i,md(j+2),k) = 2
                  system(mb(i+1),md(j+3),k) = 1
               else
                  write(*,*)'Error in jump up 2 + 2'
                  stop
               end if

c+++++++++++DVIL1 direction 6
            else if(system(mb(i-1),md(j+1),k).eq.2 .and. 
     &              system(mb(i+1),md(j-1),ml(k+1)).eq.2) then
               if(system(mb(i+1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k+1)) = 1
               else if(system(mb(i+1),md(j-3),k).eq.2 .and. 
     &              system(mb(i+2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k+1)) = 1
                  system(i,md(j-2),k) = 2
                  system(mb(i+1),md(j-3),k) = 1
               else if(system(mb(i+1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+2),j,k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k+1)) = 1
                  system(mb(i+1),md(j+1),k) = 2
                  system(mb(i+2),j,k) = 1
               else
                  write(*,*)'Error in jump up 3 + 1'
                  stop
               end if

c+++++++++++DVIL2 direction 6 right
            else if(system(mb(i-1),md(j+1),k).eq.2 .and. 
     &              system(i,md(j+2),ml(k+1)).eq.2) then
               if(system(mb(i+1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k+1)) = 1
               else if(system(mb(i+1),md(j-3),k).eq.2 .and. 
     &              system(mb(i+2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k+1)) = 1
                  system(i,md(j-2),k) = 2
                  system(mb(i+1),md(j-3),k) = 1
               else if(system(mb(i+1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+2),j,k).eq.2) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k+1)) = 1
                  system(mb(i+1),md(j+1),k) = 2
                  system(mb(i+2),j,k) = 1
               else
                  write(*,*)'Error in jump up 3 + 2'
                  stop
               end if
c+++++++++++DVIL2 direction 6 left
            else if(system(mb(i-1),md(j+1),k).eq.2 .and. 
     &              system(mb(i-1),md(j-1),ml(k+1)).eq.2) then
               if(system(mb(i+1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k+1)) = 1
               else if(system(mb(i+1),md(j-3),k).eq.2 .and. 
     &              system(mb(i+2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k+1)) = 1
                  system(i,md(j-2),k) = 2
                  system(mb(i+1),md(j-3),k) = 1
               else if(system(mb(i+1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+2),j,k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k+1)) = 1
                  system(mb(i+1),md(j+1),k) = 2
                  system(mb(i+2),j,k) = 1
               else
                  write(*,*)'Error in jump up 3 + 3'
                  stop
               end if
            else
               write(*,*)'Error in DVIL1 and DVIL2 trasfer atom up'
               stop
            end if
c...........DVIL3
         else if(system(i,j,k).eq.5 .and. 
     &           system(i,j,ml(k+1)).eq.2) then
            system(i,j,k) = 2
            system(i,j,ml(k+1)) = 1
            if(system(i,md(j+2),k).eq.2) then
               system(i,md(j-2),ml(k+1)) = 1
               if(system(mb(i-1),md(j-3),k).eq.2 .and. 
     &              system(mb(i+1),md(j-3),k).eq.1) then
                  system(mb(i-1),md(j-1),k) = 2
                  system(mb(i-1),md(j-3),k) = 1
               else if(system(mb(i-1),md(j-3),k).eq.1 .and. 
     &                 system(mb(i+1),md(j-3),k).eq.2) then
                  system(mb(i+1),md(j-1),k) = 2
                  system(mb(i+1),md(j-3),k) = 1
               end if   
            else if(system(mb(i+1),md(j+1),k).eq.2) then
               system(mb(i-1),md(j-1),ml(k+1)) = 1
               if(system(mb(i-2),j,k).eq.2 .and. 
     &              system(mb(i-1),md(j-3),k).eq.1) then
                  system(mb(i-1),md(j+1),k) = 2
                  system(mb(i-2),j,k) = 1
               else if(system(mb(i-2),j,k).eq.1 .and. 
     &                 system(mb(i-1),md(j-3),k).eq.2) then
                  system(i,md(j-2),k) = 2
                  system(mb(i-1),md(j-3),k) = 1
               end if
            else if(system(mb(i+1),md(j-1),k).eq.2) then
               system(mb(i-1),md(j+1),ml(k+1)) = 1
               if(system(mb(i-1),md(j+3),k).eq.2 .and. 
     &              system(mb(i-2),j,k).eq.1) then
                  system(i,md(j+2),k) = 2
                  system(mb(i-1),md(j+3),k) = 1
               else if(system(mb(i-1),md(j+3),k).eq.1 .and. 
     &                 system(mb(i-2),j,k).eq.2) then
                  system(mb(i-1),md(j-1),k) = 2
                  system(mb(i-2),j,k) = 1
               end if
            else if(system(i,md(j-2),k).eq.2) then
               system(i,md(j+2),ml(k+1)) = 1
               if(system(mb(i+1),md(j+3),k).eq.2 .and. 
     &              system(mb(i-1),md(j+3),k).eq.1) then
                  system(mb(i+1),md(j+1),k) = 2
                  system(mb(i+1),md(j+3),k) = 1
               else if(system(mb(i+1),md(j+3),k).eq.1 .and. 
     &                 system(mb(i-1),md(j+3),k).eq.2) then
                  system(mb(i-1),md(j+1),k) = 2
                  system(mb(i-1),md(j+3),k) = 1
               end if   
            else if(system(mb(i-1),md(j-1),k).eq.2) then
               system(mb(i+1),md(j+1),ml(k+1)) = 1
               if(system(mb(i+2),j,k).eq.2 .and. 
     &              system(mb(i+1),md(j+3),k).eq.1) then
                  system(mb(i+1),md(j-1),k) = 2
                  system(mb(i+2),j,k) = 1
               else if(system(mb(i+2),j,k).eq.1 .and. 
     &                 system(mb(i+1),md(j+3),k).eq.2) then
                  system(i,md(j+2),k) = 2
                  system(mb(i+1),md(j+3),k) = 1
               end if
            else if(system(mb(i-1),md(j+1),k).eq.2) then
               system(mb(i+1),md(j-1),ml(k+1)) = 1
               if(system(mb(i+1),md(j-3),k).eq.2 .and. 
     &              system(mb(i+2),j,k).eq.1) then
                  system(i,md(j-2),k) = 2
                  system(mb(i+1),md(j-3),k) = 1
               else if(system(mb(i+1),md(j-3),k).eq.1 .and. 
     &                 system(mb(i+2),j,k).eq.2) then
                  system(mb(i+1),md(j+1),k) = 2
                  system(mb(i+2),j,k) = 1
               end if
            else
               write(*,*)'Error in DVIL3 trasfer atom up'
               stop
            end if
         else
            write(*,*)'Error in trasfer atom up'
            stop
         end if

c............................................
c.....move to the lower layer (atom transfer)
c............................................
      else if(move(4).eq.10) then
c--------DVIL1 and DVIL2
         if(system(i,j,k).eq.6 .and. system(i,j,ml(k-1)).eq.5) then
            system(i,j,ml(k-1)) = 1
c-----------DVIL1 direction 1
            if(system(i,md(j+2),k).eq.2 .and. 
     &           system(i,md(j-2),ml(k-1)).eq.2) then
               if(system(mb(i-1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k-1)) = 1
               else if(system(mb(i-1),md(j-3),k).eq.2 .and. 
     &                 system(mb(i+1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k-1)) = 1
                  system(mb(i-1),md(j-1),k) = 2
                  system(mb(i-1),md(j-3),k) = 1
               else if(system(mb(i-1),md(j-3),k).eq.1 .and. 
     &                 system(mb(i+1),md(j-3),k).eq.2) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k-1)) = 1
                  system(mb(i+1),md(j-1),k) = 2
                  system(mb(i+1),md(j-3),k) = 1
               else
                  write(*,*)'Error in jump up 1 - 1'
                  stop
               end if
c-----------DVIL2 direction 1 right 
            else if(system(i,md(j+2),k).eq.2 .and. 
     &              system(mb(i+1),md(j+1),ml(k-1)).eq.2) then
               if(system(mb(i-1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k-1)) = 1
               else if(system(mb(i-1),md(j-3),k).eq.2 .and. 
     &                 system(mb(i+1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k-1)) = 1
                  system(mb(i-1),md(j-1),k) = 2
                  system(mb(i-1),md(j-3),k) = 1
               else if(system(mb(i-1),md(j-3),k).eq.1 .and. 
     &                 system(mb(i+1),md(j-3),k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k-1)) = 1
                  system(mb(i+1),md(j-1),k) = 2
                  system(mb(i+1),md(j-3),k) = 1
               else
                  write(*,*)'Error in jump up 1 - 2'
                  stop
               end if
c-----------DVIL2 direction 1 left
            else if(system(i,md(j+2),k).eq.2 .and. 
     &              system(mb(i-1),md(j+1),ml(k-1)).eq.2) then
               if(system(mb(i-1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k-1)) = 1
               else if(system(mb(i-1),md(j-3),k).eq.2 .and. 
     &                 system(mb(i+1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k-1)) = 1
                  system(mb(i-1),md(j-1),k) = 2
                  system(mb(i-1),md(j-3),k) = 1
               else if(system(mb(i-1),md(j-3),k).eq.1 .and. 
     &                 system(mb(i+1),md(j-3),k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k-1)) = 1
                  system(mb(i+1),md(j-1),k) = 2
                  system(mb(i+1),md(j-3),k) = 1
               else
                  write(*,*)'Error in jump up 1 - 3'
                  stop
               end if
c-----------DVIL1 direction 2
            else if(system(mb(i+1),md(j+1),k).eq.2 .and. 
     &              system(mb(i-1),md(j-1),ml(k-1)).eq.2) then
               if(system(mb(i-2),j,k).eq.1 .and. 
     &              system(mb(i-1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k-1)) = 1
               else if(system(mb(i-2),j,k).eq.2 .and. 
     &              system(mb(i-1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k-1)) = 1
                  system(mb(i-1),md(j+1),k) = 2
                  system(mb(i-2),j,k) = 1
               else if(system(mb(i-2),j,k).eq.1 .and. 
     &              system(mb(i-1),md(j-3),k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k-1)) = 1
                  system(i,md(j-2),k) = 2
                  system(mb(i-1),md(j-3),k) = 1
               else
                  write(*,*)'Error in jump up 2 - 1'
                  stop
               end if
c-----------DVIL2 direction 2 right 
            else if(system(mb(i+1),md(j+1),k).eq.2 .and. 
     &              system(mb(i+1),md(j-1),ml(k-1)).eq.2) then
               if(system(mb(i-2),j,k).eq.1 .and. 
     &              system(mb(i-1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k-1)) = 1
               else if(system(mb(i-2),j,k).eq.2 .and. 
     &              system(mb(i-1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k-1)) = 1
                  system(mb(i-1),md(j+1),k) = 2
                  system(mb(i-2),j,k) = 1
               else if(system(mb(i-2),j,k).eq.1 .and. 
     &              system(mb(i-1),md(j-3),k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k-1)) = 1
                  system(i,md(j-2),k) = 2
                  system(mb(i-1),md(j-3),k) = 1
               else
                  write(*,*)'Error in jump up 2 - 2'
                  stop
               end if
c-----------DVIL2 direction 2 left
            else if(system(mb(i+1),md(j+1),k).eq.2 .and. 
     &              system(i,md(j+2),ml(k-1)).eq.2) then
               if(system(mb(i-2),j,k).eq.1 .and. 
     &              system(mb(i-1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k-1)) = 1
               else if(system(mb(i-2),j,k).eq.2 .and. 
     &              system(mb(i-1),md(j-3),k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k-1)) = 1
                  system(mb(i-1),md(j+1),k) = 2
                  system(mb(i-2),j,k) = 1
               else if(system(mb(i-2),j,k).eq.1 .and. 
     &              system(mb(i-1),md(j-3),k).eq.2) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k-1)) = 1
                  system(i,md(j-2),k) = 2
                  system(mb(i-1),md(j-3),k) = 1
               else
                  write(*,*)'Error in jump up 2 - 2'
                  stop
               end if

c-----------DVIL1 direction 3
            else if(system(mb(i+1),md(j-1),k).eq.2 .and. 
     &              system(mb(i-1),md(j+1),ml(k-1)).eq.2) then
               if(system(mb(i-1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k-1)) = 1
               else if(system(mb(i-1),md(j+3),k).eq.2 .and. 
     &              system(mb(i-2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k-1)) = 1
                  system(i,md(j+2),k) = 2
                  system(mb(i-1),md(j+3),k) = 1
               else if(system(mb(i-1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-2),j,k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k-1)) = 1
                  system(mb(i-1),md(j-1),k) = 2
                  system(mb(i-2),j,k) = 1
               else
                  write(*,*)'Error in jump up 3 - 1'
                  stop
               end if

c-----------DVIL2 direction 3 right
            else if(system(mb(i+1),md(j-1),k).eq.2 .and. 
     &              system(i,md(j-2),ml(k-1)).eq.2) then
               if(system(mb(i-1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k-1)) = 1
               else if(system(mb(i-1),md(j+3),k).eq.2 .and. 
     &              system(mb(i-2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k-1)) = 1
                  system(i,md(j+2),k) = 2
                  system(mb(i-1),md(j+3),k) = 1
               else if(system(mb(i-1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-2),j,k).eq.2) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k-1)) = 1
                  system(mb(i-1),md(j-1),k) = 2
                  system(mb(i-2),j,k) = 1
               else
                  write(*,*)'Error in jump up 3 - 2'
                  stop
               end if
c-----------DVIL2 direction 3 left
            else if(system(mb(i+1),md(j-1),k).eq.2 .and. 
     &              system(mb(i+1),md(j+1),ml(k-1)).eq.2) then
               if(system(mb(i-1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k-1)) = 1
               else if(system(mb(i-1),md(j+3),k).eq.2 .and. 
     &              system(mb(i-2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k-1)) = 1
                  system(i,md(j+2),k) = 2
                  system(mb(i-1),md(j+3),k) = 1
               else if(system(mb(i-1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-2),j,k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k-1)) = 1
                  system(mb(i-1),md(j-1),k) = 2
                  system(mb(i-2),j,k) = 1
               else
                  write(*,*)'Error in jump up 3 - 3'
                  stop
               end if

c+++++++++++DVIL1 direction 4
            else if(system(i,md(j-2),k).eq.2 .and. 
     &           system(i,md(j+2),ml(k-1)).eq.2) then
               if(system(mb(i+1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k-1)) = 1
               else if(system(mb(i+1),md(j+3),k).eq.2 .and. 
     &                 system(mb(i-1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k-1)) = 1
                  system(mb(i+1),md(j+1),k) = 2
                  system(mb(i+1),md(j+3),k) = 1
               else if(system(mb(i+1),md(j+3),k).eq.1 .and. 
     &                 system(mb(i-1),md(j+3),k).eq.2) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k-1)) = 1
                  system(mb(i-1),md(j+1),k) = 2
                  system(mb(i-1),md(j+3),k) = 1
               else
                  write(*,*)'Error in jump up 1 + 1'
                  stop
               end if
c+++++++++++DVIL2 direction 4 right 
            else if(system(i,md(j-2),k).eq.2 .and. 
     &              system(mb(i-1),md(j-1),ml(k-1)).eq.2) then
               if(system(mb(i+1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k-1)) = 1
               else if(system(mb(i+1),md(j+3),k).eq.2 .and. 
     &                 system(mb(i-1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k-1)) = 1
                  system(mb(i+1),md(j+1),k) = 2
                  system(mb(i+1),md(j+3),k) = 1
               else if(system(mb(i+1),md(j+3),k).eq.1 .and. 
     &                 system(mb(i-1),md(j+3),k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k-1)) = 1
                  system(mb(i-1),md(j+1),k) = 2
                  system(mb(i-1),md(j+3),k) = 1
               else
                  write(*,*)'Error in jump up 1 + 2'
                  stop
               end if
c+++++++++++DVIL2 direction 4 left
            else if(system(i,md(j-2),k).eq.2 .and. 
     &              system(mb(i+1),md(j-1),ml(k-1)).eq.2) then
               if(system(mb(i+1),md(j+3),k).eq.1 .and. 
     &              system(mb(i-1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k-1)) = 1
               else if(system(mb(i+1),md(j+3),k).eq.2 .and. 
     &                 system(mb(i-1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k-1)) = 1
                  system(mb(i+1),md(j+1),k) = 2
                  system(mb(i+1),md(j+3),k) = 1
               else if(system(mb(i+1),md(j+3),k).eq.1 .and. 
     &                 system(mb(i-1),md(j+3),k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k-1)) = 1
                  system(mb(i-1),md(j+1),k) = 2
                  system(mb(i-1),md(j+3),k) = 1
               else
                  write(*,*)'Error in jump up 1 + 3'
                  stop
               end if
c+++++++++++DVIL1 direction 5
            else if(system(mb(i-1),md(j-1),k).eq.2 .and. 
     &              system(mb(i+1),md(j+1),ml(k-1)).eq.2) then
               if(system(mb(i+2),j,k).eq.1 .and. 
     &              system(mb(i+1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k-1)) = 1
               else if(system(mb(i+2),j,k).eq.2 .and. 
     &              system(mb(i+1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k-1)) = 1
                  system(mb(i+1),md(j-1),k) = 2
                  system(mb(i+2),j,k) = 1
               else if(system(mb(i+2),j,k).eq.1 .and. 
     &              system(mb(i+1),md(j+3),k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j+1),ml(k-1)) = 1
                  system(i,md(j+2),k) = 2
                  system(mb(i+1),md(j+3),k) = 1
               else
                  write(*,*)'Error in jump up 2 + 1'
                  stop
               end if
c+++++++++++DVIL2 direction 5 right 
            else if(system(mb(i-1),md(j-1),k).eq.2 .and. 
     &              system(mb(i-1),md(j+1),ml(k-1)).eq.2) then
               if(system(mb(i+2),j,k).eq.1 .and. 
     &              system(mb(i+1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k-1)) = 1
               else if(system(mb(i+2),j,k).eq.2 .and. 
     &              system(mb(i+1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k-1)) = 1
                  system(mb(i+1),md(j-1),k) = 2
                  system(mb(i+2),j,k) = 1
               else if(system(mb(i+2),j,k).eq.1 .and. 
     &              system(mb(i+1),md(j+3),k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j+1),ml(k-1)) = 1
                  system(i,md(j+2),k) = 2
                  system(mb(i+1),md(j+3),k) = 1
               else
                  write(*,*)'Error in jump up 2 + 2'
                  stop
               end if
c+++++++++++DVIL2 direction 5 left
            else if(system(mb(i-1),md(j-1),k).eq.2 .and. 
     &              system(i,md(j-2),ml(k-1)).eq.2) then
               if(system(mb(i+2),j,k).eq.1 .and. 
     &              system(mb(i+1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k-1)) = 1
               else if(system(mb(i+2),j,k).eq.2 .and. 
     &              system(mb(i+1),md(j+3),k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k-1)) = 1
                  system(mb(i+1),md(j-1),k) = 2
                  system(mb(i+2),j,k) = 1
               else if(system(mb(i+2),j,k).eq.1 .and. 
     &              system(mb(i+1),md(j+3),k).eq.2) then
                  system(i,j,k) = 2
                  system(i,md(j-2),ml(k-1)) = 1
                  system(i,md(j+2),k) = 2
                  system(mb(i+1),md(j+3),k) = 1
               else
                  write(*,*)'Error in jump up 2 + 2'
                  stop
               end if

c+++++++++++DVIL1 direction 6
            else if(system(mb(i-1),md(j+1),k).eq.2 .and. 
     &              system(mb(i+1),md(j-1),ml(k-1)).eq.2) then
               if(system(mb(i+1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k-1)) = 1
               else if(system(mb(i+1),md(j-3),k).eq.2 .and. 
     &              system(mb(i+2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k-1)) = 1
                  system(i,md(j-2),k) = 2
                  system(mb(i+1),md(j-3),k) = 1
               else if(system(mb(i+1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+2),j,k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i+1),md(j-1),ml(k-1)) = 1
                  system(mb(i+1),md(j+1),k) = 2
                  system(mb(i+2),j,k) = 1
               else
                  write(*,*)'Error in jump up 3 + 1'
                  stop
               end if

c+++++++++++DVIL2 direction 6 right
            else if(system(mb(i-1),md(j+1),k).eq.2 .and. 
     &              system(i,md(j+2),ml(k-1)).eq.2) then
               if(system(mb(i+1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k-1)) = 1
               else if(system(mb(i+1),md(j-3),k).eq.2 .and. 
     &              system(mb(i+2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k-1)) = 1
                  system(i,md(j-2),k) = 2
                  system(mb(i+1),md(j-3),k) = 1
               else if(system(mb(i+1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+2),j,k).eq.2) then
                  system(i,j,k) = 2
                  system(i,md(j+2),ml(k-1)) = 1
                  system(mb(i+1),md(j+1),k) = 2
                  system(mb(i+2),j,k) = 1
               else
                  write(*,*)'Error in jump up 3 + 2'
                  stop
               end if
c+++++++++++DVIL2 direction 6 left
            else if(system(mb(i-1),md(j+1),k).eq.2 .and. 
     &              system(mb(i-1),md(j-1),ml(k-1)).eq.2) then
               if(system(mb(i+1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k-1)) = 1
               else if(system(mb(i+1),md(j-3),k).eq.2 .and. 
     &              system(mb(i+2),j,k).eq.1) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k-1)) = 1
                  system(i,md(j-2),k) = 2
                  system(mb(i+1),md(j-3),k) = 1
               else if(system(mb(i+1),md(j-3),k).eq.1 .and. 
     &              system(mb(i+2),j,k).eq.2) then
                  system(i,j,k) = 2
                  system(mb(i-1),md(j-1),ml(k-1)) = 1
                  system(mb(i+1),md(j+1),k) = 2
                  system(mb(i+2),j,k) = 1
               else
                  write(*,*)'Error in jump up 3 + 3'
                  stop
               end if
            else
               write(*,*)'Error in DVIL1 and DVIL2 trasfer atom down'
               stop
            end if
c...........DVIL3
         else if(system(i,j,k).eq.6 .and. 
     &           system(i,j,ml(k-1)).eq.2) then
            system(i,j,k) = 2
            system(i,j,ml(k-1)) = 1
            if(system(i,md(j+2),k).eq.2) then
               system(i,md(j-2),ml(k-1)) = 1
               if(system(mb(i-1),md(j-3),k).eq.2 .and. 
     &              system(mb(i+1),md(j-3),k).eq.1) then
                  system(mb(i-1),md(j-1),k) = 2
                  system(mb(i-1),md(j-3),k) = 1
               else if(system(mb(i-1),md(j-3),k).eq.1 .and. 
     &                 system(mb(i+1),md(j-3),k).eq.2) then
                  system(mb(i+1),md(j-1),k) = 2
                  system(mb(i+1),md(j-3),k) = 1
               end if   
            else if(system(mb(i+1),md(j+1),k).eq.2) then
               system(mb(i-1),md(j-1),ml(k-1)) = 1
               if(system(mb(i-2),j,k).eq.2 .and. 
     &              system(mb(i-1),md(j-3),k).eq.1) then
                  system(mb(i-1),md(j+1),k) = 2
                  system(mb(i-2),j,k) = 1
               else if(system(mb(i-2),j,k).eq.1 .and. 
     &                 system(mb(i-1),md(j-3),k).eq.2) then
                  system(i,md(j-2),k) = 2
                  system(mb(i-1),md(j-3),k) = 1
               end if
            else if(system(mb(i+1),md(j-1),k).eq.2) then
               system(mb(i-1),md(j+1),ml(k-1)) = 1
               if(system(mb(i-1),md(j+3),k).eq.2 .and. 
     &              system(mb(i-2),j,k).eq.1) then
                  system(i,md(j+2),k) = 2
                  system(mb(i-1),md(j+3),k) = 1
               else if(system(mb(i-1),md(j+3),k).eq.1 .and. 
     &                 system(mb(i-2),j,k).eq.2) then
                  system(mb(i-1),md(j-1),k) = 2
                  system(mb(i-2),j,k) = 1
               end if
            else if(system(i,md(j-2),k).eq.2) then
               system(i,md(j+2),ml(k-1)) = 1
               if(system(mb(i+1),md(j+3),k).eq.2 .and. 
     &              system(mb(i-1),md(j+3),k).eq.1) then
                  system(mb(i+1),md(j+1),k) = 2
                  system(mb(i+1),md(j+3),k) = 1
               else if(system(mb(i+1),md(j+3),k).eq.1 .and. 
     &                 system(mb(i-1),md(j+3),k).eq.2) then
                  system(mb(i-1),md(j+1),k) = 2
                  system(mb(i-1),md(j+3),k) = 1
               end if   
            else if(system(mb(i-1),md(j-1),k).eq.2) then
               system(mb(i+1),md(j+1),ml(k-1)) = 1
               if(system(mb(i+2),j,k).eq.2 .and. 
     &              system(mb(i+1),md(j+3),k).eq.1) then
                  system(mb(i+1),md(j-1),k) = 2
                  system(mb(i+2),j,k) = 1
               else if(system(mb(i+2),j,k).eq.1 .and. 
     &                 system(mb(i+1),md(j+3),k).eq.2) then
                  system(i,md(j+2),k) = 2
                  system(mb(i+1),md(j+3),k) = 1
               end if
            else if(system(mb(i-1),md(j+1),k).eq.2) then
               system(mb(i+1),md(j-1),ml(k-1)) = 1
               if(system(mb(i+1),md(j-3),k).eq.2 .and. 
     &              system(mb(i+2),j,k).eq.1) then
                  system(i,md(j-2),k) = 2
                  system(mb(i+1),md(j-3),k) = 1
               else if(system(mb(i+1),md(j-3),k).eq.1 .and. 
     &                 system(mb(i+2),j,k).eq.2) then
                  system(mb(i+1),md(j+1),k) = 2
                  system(mb(i+2),j,k) = 1
               end if
            else
               write(*,*)'Error in DVIL3 trasfer atom down'
               stop
            end if
         else
            write(*,*)'Error in trasfer atom down'
            stop
         end if
      else
         write(*,*)move(4)
         write(*,*)'Error in move end'
         stop
      end if

      end 

