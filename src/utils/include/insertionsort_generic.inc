! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


do i = 2, xlen
  tmp = x(i)
  j = i
  1 continue
    if (j > 1 .and. x(j-1) > tmp) then
      x(j) = x(j-1)
      j = j - 1
      goto 1
    end if
    
    x(j) = tmp
end do
