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
c     + enGrid is distributed in the hope that it will be useful,            +
c     + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
c     + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
c     + GNU General Public License for more details.                         +
c     +                                                                      +
c     + You should have received a copy of the GNU General Public License    +
c     + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
c     +                                                                      +
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     This is a simple simulation of a shock tube.
c     Goal of this program is to provide test data for the DrNUM shared memory library

      PROGRAM shock_tube
      IMPLICIT NONE
      INTEGER NUM_NODES, i, file
      INTEGER shm, mutex, barrier
      PARAMETER (NUM_NODES = 10000)
      REAL*8 L, DX, CV, CP, R, GAMMA, CFL, DT, max_time, time, V
      REAL*8 write_time, write_interv
      REAL*8 flux_rho_f, flux_rhou_f, flux_rhoE_f, a_f
      REAL*8 flux_rho_t, flux_rhou_t, flux_rhoE_t, a_t
      REAL*8 flux_rho, flux_rhou, flux_rhoE
      REAL*8 psi_f, psi_t, alpha, omega, pflux
      PARAMETER (L = 10.0)
      PARAMETER (CP = 1004.5)
      PARAMETER (R  = 287.0)
      PARAMETER (CFL = 0.1)

      REAL*8 X(NUM_NODES)
      REAL*8 rho(NUM_NODES)
      REAL*8 rhou(NUM_NODES)
      REAL*8 rhoE(NUM_NODES)
      REAL*8 res_rho(NUM_NODES)
      REAL*8 res_rhou(NUM_NODES)
      REAL*8 res_rhoE(NUM_NODES)
      REAL*8 Ma(NUM_NODES)
      REAL*8 p(NUM_NODES)
      REAL*8 u(NUM_NODES)
      REAL*8 T(NUM_NODES)

      DX = L/(NUM_NODES-1)
      CV = CP - R
      GAMMA = CP/CV
      DT = CFL*DX/SQRT(GAMMA*R*300.0)

      write (*,*) 'DT=', DT
      write (*,*) 'sound wave travel time=', L/SQRT(gamma*R*300.0)

      CALL shm_new(shm, 1, 6*(NUM_NODES*8 + 1000), 1)
      CALL mutex_new(mutex, 2, 1)
      CALL barrier_new(barrier, 3, 1)

      DO 100 i = 1, NUM_NODES
        X(i) = (i-1)*DX
        p(i) = 1e6
        IF (i.GT.NUM_NODES/2) p(i) = 1e5
        T(i) = 300.0
        Ma(i) = 0.0
        u(i) = 0.0
        rho(i) = p(i)/(R*T(i))
        rhou(i) = 0.0
        rhoE(i) = rho(i)*(CV*T(i) + 0.5*U(i)**2)
100   CONTINUE

      time = 0.0
      write_time = 0.0
      CALL mutex_lock(mutex)
      CALL shm_write_real_value(shm, 'time.', time)
      CALL shm_write_real_array(shm, 'x.', NUM_NODES, X)
      CALL shm_write_real_array(shm, 'p.', NUM_NODES, p)
      CALL shm_write_real_array(shm, 'u.', NUM_NODES, u)
      CALL shm_write_real_array(shm, 'Ma.', NUM_NODES, Ma)
      CALL shm_write_real_array(shm, 'density.', NUM_NODES, rho)
      CALL shm_write_real_array(shm, 'T.', NUM_NODES, T)
      CALL mutex_unlock(mutex)
      CALL barrier_wait(barrier)
      write (*,*) 'Please enter max time!'
      read (*,*), max_time
      write (*,*) 'Please write interval!'
      read (*,*), write_interv

200   DO 300 i = 1, NUM_NODES
        res_rho(i)  = 0.0
        res_rhou(i) = 0.0
        res_rhoE(i) = 0.0
        u(i) = rhou(i)/rho(i)
        T(i) = (rhoE(i)/rho(i) - 0.5*u(i)**2)/CV
        p(i) = rho(i)*T(i)*R
300   CONTINUE
      res_rhou(1) = p(1)
      res_rhou(NUM_NODES) = -p(NUM_NODES)
      DO 400 i = 1, NUM_NODES - 1
        a_f = SQRT(GAMMA*R*T(i))
        a_t = SQRT(GAMMA*R*T(i+1))
        flux_rho_f = rho(i)
        flux_rhou_f = rhou(i)
        flux_rhoE_f = rhoE(i) + p(i)
        flux_rho_t = rho(i+1)
        flux_rhou_t = rhou(i+1)
        flux_rhoE_t = rhoE(i+1) + p(i+1)
        psi_f = MAX(MAX(a_f + u(i), a_t + u(i+1)), 0.0)
        psi_t = MAX(MAX(a_f - u(i), a_t - u(i+1)), 0.0)
        alpha = 0.5
        omega = alpha*MAX(psi_f, psi_t)
        pflux = alpha*p(i) + (1.0-alpha)*p(i+1);
        flux_rho  =   alpha*u(i)*flux_rho_f
     #              + (1.0-alpha)*u(i+1)*flux_rho_t
     #              - omega*(flux_rho_t  - flux_rho_f)
        flux_rhou =   alpha*u(i)*flux_rhou_f
     #              + (1.0-alpha)*u(i+1)*flux_rhou_t
     #              - omega*(flux_rhou_t - flux_rhou_f) + pflux
        flux_rhoE =   alpha*u(i)*flux_rhoE_f
     #              + (1.0-alpha)*u(i+1)*flux_rhoE_t
     #              - omega*(flux_rhoE_t - flux_rhoE_f)
        res_rho(i)  = res_rho(i)  - flux_rho
        res_rhou(i) = res_rhou(i) - flux_rhou
        res_rhoE(i) = res_rhoE(i) - flux_rhoE
        res_rho(i+1)  = res_rho(i+1)  + flux_rho
        res_rhou(i+1) = res_rhou(i+1) + flux_rhou
        res_rhoE(i+1) = res_rhoE(i+1) + flux_rhoE
400   CONTINUE
      DO 500 i = 1, NUM_NODES
        V = DX
        if ((i.EQ.1).OR.(i.EQ.NUM_NODES)) V = 0.5*DX
        rho(i)  = rho(i)  + DT/V*res_rho(i)
        rhou(i) = rhou(i) + DT/V*res_rhou(i)
        rhoE(i) = rhoE(i) + DT/V*res_rhoE(i)
        u(i) = rhou(i)/rho(i)
        T(i) = (rhoE(i)/rho(i) - 0.5*u(i)**2)/CV
        p(i) = rho(i)*T(i)*R
        Ma(i) = u(i)/SQRT(GAMMA*R*T(i))
500   CONTINUE
      time = time + DT
      write_time = write_time + DT
      if (write_time.LT.write_interv) goto 600
        write (*,*) 'time=', time
        CALL mutex_lock(mutex)
        CALL shm_write_real_value(shm, 'time.', time)
        CALL shm_write_real_array(shm, 'x.', NUM_NODES, X)
        CALL shm_write_real_array(shm, 'p.', NUM_NODES, p)
        CALL shm_write_real_array(shm, 'u.', NUM_NODES, u)
        CALL shm_write_real_array(shm, 'Ma.', NUM_NODES, Ma)
        CALL shm_write_real_array(shm, 'density.', NUM_NODES, rho)
        CALL shm_write_real_array(shm, 'T.', NUM_NODES, T)
        CALL mutex_unlock(mutex)
        CALL barrier_wait(barrier)
        write_time = 0.0
600   if (time.LT.max_time) goto 200
      write (*,*) 'Please enter max time!'
      read (*,*), max_time
      if (time.LT.max_time) goto 200

      CALL shm_delete(shm)
      CALL mutex_delete(mutex)
      CALL barrier_delete(barrier)

      STOP
      END







