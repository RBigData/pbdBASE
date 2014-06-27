! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


ret = l

pvt = x(ind)

call swap(x, ind, r)

do i = l, r
  if (x(i).lt.pvt) then
    call swap(x, i, ret)
    ret = ret + 1
  end if
end do

call swap(x, ret, r)

ret = ret - 1 ! for "<=" return
