       subroutine CutStr(Line,NumLin,LinPos,LinEnd,NFIL,iErr)
c.................................................................
c....It reads the line from the file UNIT=NFIL as a string   .....
c.... Line*80 (if NFIL.ne.0) and replies with a set of       .....
c.... substrings. Their number is NumLin, while starting and .....
c.... ending positions are LinPos and LinEnd.                .....
c.................................................................
c Up to 30 substrings is implied to be in the initial string Line.
c.................................................................
c  If NFIL=0 it is supposed that the string Line*80 exists and
c  must not be read from the external source, but must be taken
c  as an input from the parent program
c.................................................................
c  iErr=1 in the case of error input
c  iErr=2 in the case of end of file
c  iErr=0 otherwise
c.................................................................
c
      character Line*120,cha
      integer LinEnd(30),LinPos(30)
      iErr=0
c
c............ get the line skipping comments
      if(NFIL.ne.0) then
         call SkipL(NFIL,Line,iErr)
         if(iErr.ne.0) return
      end if
      NumLin=0
      icase=0
c............ Loop over all characters in the Line(1:80)............
       do 5 i=1,120
       cha=Line(i:i)
c________ if the 1st not blank character is found
       if(cha.ne.' '.and.icase.eq.0) then
          NumLin=NumLin+1
          if(NumLin.gt.30) go to 10
          LinPos(NumLin)=i
          icase=1
       end if
c_________ if the 1st blank character is found after the substring
       if(cha.eq.' '.and.icase.eq.1) then
          LinEnd(NumLin)=i-1
          icase=0
       end if
 5     end do
       if(icase.eq.1) LinEnd(NumLin)=120
       return
10     iErr=1
       return
       end

      subroutine SkipL(NFIL,Line,iErr)
c.......................................................................
c...... Program skips all "empty" string containing comments, etc.......
c...... All lines, containing **,--,==,*- and -* symbols are     .......
c...... recognized as the comments lines and are skipped. The    .......
c...... driver is set up at the beginning of the first "nonempty".......
c...... record                                                   .......
c.......................................................................
c  iErr=1 in the case of error input
c  iErr=2 in the case of end of file
c  iErr=0 otherwise
c.................................................................
c
      parameter (nComm=6)
      character Line*120,Comm(nComm)*2
      data Comm/'*-','-*','--','==','**','##'/
c.......... read the current line from the driver NFIL
1     read (NFIL,'(a)',err=10,end=20) Line
c_________ analyze Line if it is empty
      if(Line(1:1).eq.CHAR(0)) go to 1
c_________ analyze Line if it is filled by simple blanks
        do i=1,120
        if(Line(i:i).ne.CHAR(0).and.Line(i:i).ne.' ') go to 2
        end do
      go to 1
c_________ analyze Line for the comment symbols; in the case if any
c          one from the list Comm() was found it reads the next line
2       do i=1,nComm
        i0=INDEX(Line,Comm(i))
        if(i0.ne.0) go to 1
        end do
      iErr=0
      return
c......... in the case of an error
10    iErr=1
      return
 20   iErr=2
      return
      end

      subroutine find_string(string,l,NFIL,iErr)
      character*(*) string
      character line*120
      logical memo
c
      iErr=0
 10   read (NFIL,'(a)',END=40,ERR=40) Line
      i0=index(Line,string(1:l))
      if(i0.eq.0) go to 10
      write(9,'(/a)')
     & '.........O.K.! String '//string(1:l)//' is found !!! .........'
      return
 40   iErr=1
      write(9,'(/a)')
     & '.........BAD! String '//string(1:l)//' NOT found !!! .........'
      end

      subroutine right(char,len,len00)
      character*(*)  char
      do i=len,1,-1
         j=i
         if(char(i:i).ne.' ') go to 1
      end do
 1    jj=len-j
      do i=j,1,-1
         char(jj+i:jj+i) = char(i:i)
      end do
      do i=1,jj
         char(i:i)=' '
      end do
      len00=jj+1
      return
      end
