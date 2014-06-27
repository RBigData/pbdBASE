! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


! bubble sort with a poor man's while loop
1 continue
  do i = 1, 2
    if (tmp(i) > tmp(i+1)) then
        call swap(tmp, i, i+1)
        goto 1
    end if
  end do

med = tmp(2)
