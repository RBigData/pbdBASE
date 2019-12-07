function fpmod(a, b) &
result(ret)
  ! use, intrinsic :: ieee_arithmetic
  double precision, intent(in) :: a, b
  double precision :: ret
  double precision, parameter :: zero = 0.0d0
  
  
  if (abs(b) < 1.0d-8) then
    ret = nan !ieee_value(1.0d0, ieee_quiet_nan)
  else if (b > 0) then
    if (a >= 0) then
      ret = dmod(a, b)
    else
      ret = b - dmod(-a, b)
    end if
  else
    if (abs(a) < 1.0d-8) then
      ret = zero
    else if (a > 0) then
      ret = b + dmod(a, -b)
    else
      ret = -dmod(-a, -b)
    end if
  end if
  
  return
end function
