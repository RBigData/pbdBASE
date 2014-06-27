function fpmod(a, b) &
result(ret)
  double precision, intent(in) :: a, b
  double precision :: ret
  double precision, parameter :: zero = 0.0d0
  
  
  if (b == 0) then
    ret = nan
  else if (b > 0) then
    if (a >= 0) then
      ret = dmod(a, b)
    else
      ret = b - dmod(-a, b)
    end if
  else
    if (a == 0) then
      ret = zero
    else if (a > 0) then
      ret = b + dmod(a, -b)
    else
      ret = -dmod(-a, -b)
    end if
  end if
  
  return
end function

