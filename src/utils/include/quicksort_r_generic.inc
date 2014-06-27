! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013-2014, Schmidt


if (l < r) then
  ! choose median of l, r, and (l+r)/2 as pivot
  tmp = (/ l, (l+r)/2, r /)
  pvt = quicksort_median_of_3(tmp)
  
  ! partition x by pvt
  if (tmp(1) == pvt) then
    ind = l 
  else if (tmp(2) == pvt) then
    ind = (l+r)/2
  else
    ind = r
  end if
  
  call swap(x, ind, r)
  
  call quicksort_partition(x, l, r, ind, itmp)

