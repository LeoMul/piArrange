program parArrange 
    !arranges pi files
    !leo mulholland 13.04.25

    !implicit none 
    integer,parameter :: max_iter = 9999
    integer :: i,j
    character*10 :: currentfile
    logical :: ex 

    integer :: existantFiles(10000)
    integer :: numFiles,fileindex
    integer,parameter :: fileoffset = 100

    integer :: nelec,nzed,nast,masterNumpoints,headInt,numbound
    integer :: localnumpoints 

    real*8,allocatable :: ionizedBoundstates(:)
    integer,allocatable :: list2J(:)
    integer,allocatable :: indexpointer(:)

    real*8 :: xpisum 
    real*8,allocatable :: energyRange(:)
    real*8,allocatable :: crosssection(:,:)
    integer :: offset 
    integer :: temp 

    numFiles = 0
    masterNumpoints = 0
    do i = 0,max_iter
        currentfile = file_name(i)
        inquire(file=currentfile,exist=ex)
        if(ex) then 
            numFiles = numFiles + 1
            existantFiles(numFiles) = i 
            open(numFiles,file=currentfile,form='formatted')
            call readHeaderInFile(&
                numFiles,&
                nelec,&
                nzed,&
                nast,&
                localnumpoints,&
                headInt,&
                numbound&
            )
            masterNumpoints = masterNumpoints + localnumpoints
        end if 
    end do

    write(0,*),'Found ',numFiles,' files.'
    write(0,*),'Found ',numbound,' N+1 bound states.'
    write(0,*),'Commencing loop over bound states.'
    allocate(energyRange(masterNumpoints))
    allocate(crosssection(masterNumpoints,nast))
    allocate(indexpointer(masterNumpoints))
    energyRange = 0.0d0 
    crosssection = 0.0d0
    indexpointer = 0 

    do i = 1,numbound
        write(0,*) "reading bound state",i
        offset = 0
        !collect all of the data 
        do j = 1,numFiles
            call readCrossSectionInFile(&
                j,&
                i,&
                energyRange,&
                crosssection,&
                offset&
            )
        end do 
        !stop 
        !now call a quick sort 
        do j =1,masterNumpoints
            indexpointer(j) = j
        end do 

        !call qsort(energyRange,masterNumpoints,indexpointer)

        !dump out this bound state
        temp = 100000
        open(temp,file='dump'//file_name(i))
        do j = 1,masterNumpoints    
            xpisum = sum(crosssection(indexpointer(j),:))

            write(temp,1) energyRange(j),xpisum, (crosssection(indexpointer(j),k),k=1,nast)
        end do 
        close(temp)

    end do 

    !close all the files
    do i = 1,numFiles
        fileindex = existantFiles(i)
        close(fileindex)
    end do 
    1   format(1PE14.8,10000(1PE11.3))

    contains  

    subroutine readHeaderInFile(&
        filenumber,&
        nelec,&
        nzed,&
        nast,&
        masterNumpoints,&
        headInt,&
        numbound&
        )
        integer :: nelec,nzed,nast,masterNumpoints,headInt,numbound
        integer :: filenumber
        integer :: j,i

        !read
        read(filenumber,*) nzed,nelec
        read(filenumber,*) nast,masterNumpoints,headInt,numbound
        !to do, make a consitency check. 

        !only allocate the first time
        !could code this better but its 4pm on a friday
        if(filenumber.eq.1) then 
            allocate(ionizedBoundstates(nast))
            allocate(list2J(nast))
        end if

        read(filenumber,*) (j,list2J(i),i=1,nast)
        read(filenumber,*) (ionizedBoundstates(i),i=1,nast)

    end subroutine

    subroutine readCrossSectionInFile(&
        filenumber,&
        boundnumber,&
        energyRange,&
        crosssection,&
        offset&
        )   
        implicit none 
        integer,intent(in) ::filenumber, boundnumber
        real*8 :: energyRange(:),crosssection(:,:)
        integer :: l,m,n,o,k
        real*8 :: ip,ee,ef,diff
        integer :: mynum,localStates,offset
        real*8,allocatable :: diffarray(:)
        read(filenumber,*) l,m,n,o
        read(filenumber,*) ip,mynum
        allocate(diffarray(nast))
        diffarray = 0.0d0 
        do l=1,mynum

            read(filenumber,1,advance='no') ef
            ee = ef - abs(ip)
            energyRange(l+offset) = ef 
            localStates = 0 
            !determine how many CS to read for this energy.
            do k = 1,nast 
                diff = ionizedBoundstates(k)*(nzed-nelec)**2 - ee 
                diffarray(k) = diff 
                if (diff.le.1e-6) then 
                    localStates = localStates + 1
                    if (diff.gt.0.0d0) then 
                        write(0,*) "    possible energy close to threshold."
                        write(0,*) "     Bound energy: ",ionizedBoundstates(k)*(nzed-nelec)**2
                        write(0,*) "    Photon energy: ",ee

                    end if 
                else 
                    exit 
                end if 
            end do 
            read(filenumber,2) (crosssection(l+offset,k),k=1,localStates)

        end do
        !Read them.
        offset = offset + mynum 

        1   format(1PE14.8)
        2   format(6(1PE11.3)/(14X,6(E11.3)))
    end subroutine

    function file_name(number)
        !finds omega file name from an integer.
        character*10 :: file_name
        integer :: number
        integer :: i1,i2,i3,i4
        character*1 NUM(0:9),filec,filed,fileu,filev
        DATA NUM /'0','1','2','3','4','5','6','7','8','9'/

        i1= number/1000
        i2=(number-1000*(number/1000))/100
        i3=(number-(1000*(number/1000))-i2*100)/10
        i4=(number-(1000*(number/1000))-i2*100)-i3*10

        filec = NUM(i1)
        filed = NUM(i2)
        fileu = NUM(i3)
        filev = NUM(i4)

        file_name = 'XPIPAR'//filec//filed//fileu//filev
        
    end function

    SUBROUTINE qsort(a, n, t)

        !I copy and pasted this from the original code. L

        !    NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
        !    BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.
        !    REAL*8 PRECISION, ALSO CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.
                IMPLICIT NONE
                INTEGER, INTENT(IN)    :: n
                REAL*8, INTENT(INOUT)    :: a(n)
                INTEGER, INTENT(INOUT) :: t(n)
        !    Local Variables
                INTEGER    :: i, j, k, l, r, s, stackl(15), stackr(15), ww
                REAL*8    :: w, x
                s = 1
                stackl(1) = 1
                stackr(1) = n
        !    KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.
         10     CONTINUE
                l = stackl(s)
                r = stackr(s)
                s = s - 1
        !    KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.
         20     CONTINUE
                i = l
                j = r
                k = (l+r) / 2
                x = a(k)
        !    REPEAT UNTIL I > J.
              DO
              DO
               IF (a(i).LT.x) THEN    ! Search from lower end
              i = i + 1
              CYCLE
               ELSE
               EXIT
               END IF
              END DO
              DO
              IF (x.LT.a(j)) THEN    ! Search from upper end
              j = j - 1
              CYCLE
              ELSE
              EXIT
              END IF
              END DO
              IF (i.LE.j) THEN    ! Swap positions i & j
              w = a(i)
              ww = t(i)
              a(i) = a(j)
              t(i) = t(j)
              a(j) = w
              t(j) = ww
              i = i + 1
              j = j - 1
              IF (i.GT.j) EXIT
              ELSE
              EXIT
              END IF
              END DO
              IF (j-l.GE.r-i) THEN
              IF (l.LT.j) THEN
              s = s + 1
              stackl(s) = l
              stackr(s) = j
              END IF
              l = i
              ELSE
              IF (i.LT.r) THEN
              s = s + 1
              stackl(s) = i
              stackr(s) = r
              END IF
              r = j
              END IF
              IF (l.LT.r) GO TO 20
              IF (s.NE.0) GO TO 10
              RETURN
    END SUBROUTINE qsort

end program 