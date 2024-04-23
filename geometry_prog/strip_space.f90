subroutine strip_space(string)

character(len=*) :: string
integer :: stringLen, last, actual

stringLen = len(string)
last = 1
actual = 1

do while (actual < stringLen) 
   if (string(last:last) == ' ') then
      actual = actual + 1
      string(last:last) = string(actual:actual)
      string(actual:actual) = ' '
   else
      last = last + 1 
      if (actual < last) actual = last
   endif
enddo

end subroutine
