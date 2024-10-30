!$Header: /usr/local/ollincvs/Codes/OllinAxis-BiB/src/base/functions.f90,v 1.2 2021/03/08 01:00:01 malcubi Exp $

! ************************************
! ***   VARIOUS USEFUL FUNCTIONS   ***
! ************************************

! Here I define various useful functions:
!
! contains        Checks if a string contains a specific pattern (in a list).
! message         Created a properly formatted message to screen.



! *****************************
! ***   FUNCTION CONTAINS   ***
! *****************************

  logical function contains(string,pattern)

! This function checks if a given string contains a specific pattern
! in a comma or space separated list, that is, the pattern must be
! separated by commas and/or spaces from other elements in the list.

  integer i,l1,l2

  character(len=*) string,pattern

! Initialize.

  contains = .false.

! Find out length of strings and check if the
! pattern is smaller than the string.

  l1 = len(string)
  l2 = len(pattern)

! print *, pattern,string

  if (l2>l1) return

! Check if there are commas or spaces in the pattern.
! This is not allowed so the code complains and stops.

  do i=1,l2
     if ((pattern(i:i)==",").or.(pattern(i:i)==" ")) then
       stop 'Error in function "contains". The pattern is not allowed to have commas or spaces.\n'
     end if
  end do

! If the pattern is equal to the full string then
! the list has only one element.

  if (string==pattern) then
     contains = .true.
     return
  end if

! Now check if the pattern is contained in string.

  do i=1,l1-l2+1

     if (string(i:i+l2-1)==pattern) then
 
!       If the pattern is in the string check for commas or spaces.

        if ((i==1).and. &
            ((i+l2-1==l1).or.(string(i+l2:i+l2)==",").or.(string(i+l2:i+l2)==" "))) then
           contains = .true.
        else if (((string(i-1:i-1)==",").or.(string(i-1:i-1)==" ")).and. &
            ((i+l2-1==l1).or.(string(i+l2:i+l2)==",").or.(string(i+l2:i+l2)==" "))) then
           contains = .true.
        end if

     end if

  end do

! End.

  end function contains






! ****************************
! ***   FUNCTION MESSAGE   ***
! ****************************

  character(len=*) function message(string1,string2)

! This function creates a properly formatted message
! to be sent out to the screen taking the two input
! strings (string1,string2).

  implicit none

  logical space
  integer i,j,l

  character(len=*) string1,string2
  character(60) work,work2

! Flag for extra space.

  space = .false.

! Find out length of full concatenated string.

  if (string1==" ") then
     l = len(trim(string2))
  else if (string2==" ") then
     l = len(trim(string1))
  else
     l = len(trim(string1)//" "//trim(string2))
  end if

! Length is even, string is OK.

  if (real(l/2)==real(l)/2.d0) then

     j = l/2

     if (string1==" ") then
        work = trim(string2)
     else if (string2==" ") then
        work = trim(string1)
     else
        work = trim(string1)//" "//trim(string2)
     end if

! Length is odd, string needs extra space.

  else

     j = l/2 + 1

     if (string1==" ") then
        space = .true.
        work = trim(string2)//" "
     else if (string2==" ") then
        space = .true.
        work = trim(string1)//" "
     else
        work = trim(string1)//"  "//trim(string2)
     end if

  end if

! Construct message. Notice that the message is constructed
! from right to left.

  work2 = "***"

  do i=1,18-j
     work2 = " "//work2
  end do

  if (space) then
     work2 = trim(work)//" "//work2
  else
     work2 = trim(work)//work2
  end if

  do i=1,18-j
     work2 = " "//work2
  end do

  message = trim("***"//work2)

! End.

  end function message


