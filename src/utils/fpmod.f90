function fpmod(a, b) result(ret)
  double precision, intent(in) :: a, b
  double precision :: ret
  
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
      ret = 0.0d0
    else if (a > 0) then
      ret = b + dmod(a, -b)
    else
      ret = -dmod(-a, -b)
    end if
  end if
  
  return
end function
