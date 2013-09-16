c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     +                                                                      +
c     + This file is part of DrNUM.                                          +
c     +                                                                      +
c     + Copyright 2013 numrax GmbH, enGits GmbH                              +
c     +                                                                      +
c     + DrNUM is free software: you can redistribute it and/or modify        +
c     + it under the terms of the GNU General Public License as published by +
c     + the Free Software Foundation, either version 3 of the License, or    +
c     + (at your option) any later version.                                  +
c     +                                                                      +
c     + DrNUM is distributed in the hope that it will be useful,             +
c     + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
c     + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
c     + GNU General Public License for more details.                         +
c     +                                                                      +
c     + You should have received a copy of the GNU General Public License    +
c     + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
c     +                                                                      +
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     This is a simple test program

      PROGRAM ShmTest
      IMPLICIT NONE

      INTEGER   shm, mutex
      INTEGER   i
      REAL*8    real_data(21)
      REAL*8    real_sum
      INTEGER   int_data(21)
      INTEGER   int_sum
      INTEGER   int_value

      CHARACTER*200 text

      CALL shm_new(shm, 1, 10000, 1)
      CALL mutex_new(mutex, 1, 1)

      real_sum = 1
      real_data(1) = 1
      int_sum = 10
      int_data(1) = 10
      do 100 i = 2, 20
        real_data(i) = real_data(i-1) + 1
        real_sum = real_sum + real_data(i)
        int_data(i) = int_data(i-1) + 1
        int_sum = int_sum + int_data(i)
100   continue
      write (*,*) 'real_sum = ', real_sum
      write (*,*) 'int_sum  = ', int_sum
      text = 'Hello World!'

      CALL mutex_lock(mutex)
      CALL shm_write_real_array(shm, 'real_data.', 20, real_data)
      CALL shm_write_integer_array(shm, 'int_data.', 20, int_data)
      CALL shm_write_real_value(shm, 'real_sum.', real_sum)
      CALL shm_write_integer_value(shm, 'int_sum.', int_sum)
      CALL shm_write_character_array(shm, 'text1.', 200, text)
      CALL shm_write_integer_value(shm, 'int_value.', 1)
      CALL mutex_unlock(mutex)

200   read (*,*), int_value
      CALL mutex_lock(mutex)
      CALL shm_write_integer_value(shm, 'int_value.', int_value)
      CALL mutex_unlock(mutex)
      if (int_value.EQ.0) goto 300
      goto 200

300   CALL shm_delete(shm)
      CALL mutex_delete(mutex)
      STOP
      END
