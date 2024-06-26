!-------------------------------------------------
! bcflag = 1  => dirichlet bc
! bcflag = 2  => neumann bc
! bcflag = 3  => pulsatile bc
!
!     |-------|-------|-------|-------|-------
!     |       |       |       |       |       
!     |   o   |   o   |   o   |   o   |   o   
!     |       |  bcy  |  bcy  |       |       
!     |-------|---+---|-------|-------|-------
!     |       |*******|*******|       |       
!     |   obcx+*******|*******|bcxo   |   o   
!     |       |*******|*******|       |       
!     |-------|-------|-------|-------|-------
!
!
!-------------------------------------------------
   SUBROUTINE set_outer_velocity_bc()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE MPI_module
    
    IMPLICIT NONE

    INTEGER             :: i,j,ii,jj,kk,k
    REAL(KIND=CGREAL)   :: temp


!--------------------------
! Outer boundary conditions
!--------------------------


!---
! left boundary
!---
        ii= 0
        i = 1
        SELECT CASE (bcx1)
          CASE (BC_TYPE_DIRICHLET)      ! dirichlet bc
             DO k=zb1,zb2   !0,nz  
             DO j=0,ny 
                bcxu(ii,j,k) = ux1!*ramp_factor
                bcxv(ii,j,k) = vx1!*ramp_factor
                bcxw(ii,j,k) = wx1!*ramp_factor
             ENDDO
             ENDDO
          
          CASE (BC_TYPE_ZERO_GRADIENT)  ! outflow bc ( zero gradient ; explicit)  
             DO k=zb1,zb2   !0,nz  
             DO j=0,ny 
                bcxu(ii,j,k) = u(i,j,k)
                bcxv(ii,j,k) = v(i,j,k)
                bcxw(ii,j,k) = w(i,j,k)
             ENDDO
             ENDDO

          CASE (BC_TYPE_PULSATILE_INFLOW)
             DO k=zb1,zb2   !0,nz  
             DO j=0,ny   
                bcxu(ii,j,k) = ux1*sin(2.0_CGREAL*pi*freq_ux1*time)
                bcxv(ii,j,k) = vx1*sin(2.0_CGREAL*pi*freq_vx1*time)
                bcxw(ii,j,k) = wx1*sin(2.0_CGREAL*pi*freq_wx1*time)
             ENDDO
             ENDDO
            
          CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
             DO k=zb1,zb2   !0,nz
             DO j=0,ny             
                bcxu(ii,j,k) = 0.0_CGREAL
                bcxv(ii,j,k) = v(i,j,k)
                bcxw(ii,j,k) = w(i,j,k)
             ENDDO
             ENDDO

          CASE (BC_TYPE_PERIODIC)       ! periodic bc  (explicit & dirty implementation)
             DO k=zb1,zb2   !0,nz  
             DO j=0,ny 
                bcxu(ii,j,k) = 0.5_CGREAL*( u(i,j,k) + u(nx-1,j,k) )
                bcxv(ii,j,k) = 0.5_CGREAL*( v(i,j,k) + v(nx-1,j,k) )
                bcxw(ii,j,k) = 0.5_CGREAL*( w(i,j,k) + w(nx-1,j,k) )
             ENDDO
             ENDDO

          CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
             DO k=zb1,zb2   !0,nz   
             DO j=0,ny
                ! By F.-B. Tian, use the channel flow
                bcxu(ii,j,k) = ux1*6.0*(y(j)-yOrigin)*(yout-y(j))/(yout-yOrigin)**2
                bcxv(ii,j,k) = vx1 
                bcxw(ii,j,k) = wx1
             ENDDO
             ENDDO
          case(BC_TYPE_TAYLOR_GREEN)
             DO k=zb1, zb2
                DO j=0, ny
                  bcxu(ii, j, k) = -0.5*(dcos(xc(i-1)) * dsin(yc(j)) * exp(-2*reinv*time) + dcos(xc(i)) * dsin(yc(j)) * exp(-2*reinv*time))
                  bcxv(ii, j, k) = 0.5*(dsin(xc(i-1)) * dcos(yc(j)) * exp(-2*reinv*time) + dsin(xc(i)) * dcos(yc(j)) * exp(-2*reinv*time))
                  bcxw(ii, j, k) = 0.0_CGREAL
                ENDDO
             ENDDO
          CASE (BC_TYPE_SHEAR)       ! shear bc 
            print*, 'Not implemented!'
            stop
            !DO k=zb1,zb2   !0,nz 
            !DO j=0,ny
               !            bcxu(ii,j,k) = ux1*y(j)/yout
               !            bcxv(ii,j,k) = vx1 
               !            bcxw(ii,j,k) = wx1  
            !ENDDO
            !ENDDO            
        END SELECT 

!---
! right boundary
!---
        ii= 1
        i = nx-1
        SELECT CASE (bcx2)        
        CASE (BC_TYPE_DIRICHLET)      ! dirichlet bc              
           DO k=zb1,zb2   !0,nz   
           DO j=0,ny
              bcxu(ii,j,k) = ux2 
              bcxv(ii,j,k) = vx2 
              bcxw(ii,j,k) = wx2 
           ENDDO
           ENDDO
        
        CASE (BC_TYPE_ZERO_GRADIENT)  ! outflow bc ( zero gradient ; explicit)
           DO k=zb1,zb2   !0,nz   
           DO j=0,ny
              bcxu(ii,j,k) = u(i,j,k)
              bcxv(ii,j,k) = v(i,j,k)
              bcxw(ii,j,k) = w(i,j,k)
           ENDDO
           ENDDO
                
        CASE (BC_TYPE_PULSATILE_INFLOW)         
           DO k=zb1,zb2   !0,nz   
           DO j=0,ny
              bcxu(ii,j,k) = ux2*sin(2.0_CGREAL*pi*freq_ux2*time)
              bcxv(ii,j,k) = vx2*sin(2.0_CGREAL*pi*freq_vx2*time)
              bcxw(ii,j,k) = wx2*sin(2.0_CGREAL*pi*freq_wx2*time)
           ENDDO
           ENDDO
        
        CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
           DO k=zb1,zb2   !0,nz   
           DO j=0,ny
              bcxu(ii,j,k) = 0.0_CGREAL
              bcxv(ii,j,k) = v(i,j,k)
              bcxw(ii,j,k) = w(i,j,k)
           ENDDO
           ENDDO

        CASE (BC_TYPE_PERIODIC)       ! periodic bc
           DO k=zb1,zb2   !0,nz   
           DO j=0,ny
              bcxu(ii,j,k) = 0.5_CGREAL*( u(i,j,k) + u(1,j,k) )
              bcxv(ii,j,k) = 0.5_CGREAL*( v(i,j,k) + v(1,j,k) )
              bcxw(ii,j,k) = 0.5_CGREAL*( w(i,j,k) + w(1,j,k) )
           ENDDO
           ENDDO

        CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
           DO k=zb1,zb2   !0,nz   
           DO j=0,ny
              ! By F.-B. Tian, use the channel flow
              bcxu(ii,j,k) = ux2*6.0*(y(j)-yOrigin)*(yout-y(j))/(yout-yOrigin)**2
              bcxv(ii,j,k) = vx2 
              bcxw(ii,j,k) = wx2
           ENDDO
           ENDDO
        case(BC_TYPE_TAYLOR_GREEN)
             DO k=zb1, zb2
                DO j=0, ny
                  bcxu(ii, j, k) = -0.5*(dcos(xc(i+1)) * dsin(yc(j)) * exp(-2*reinv*time) + dcos(xc(i)) * dsin(yc(j)) * exp(-2*reinv*time))
                  bcxv(ii, j, k) = 0.5*(dsin(xc(i+1)) * dcos(yc(j)) * exp(-2*reinv*time) + dsin(xc(i)) * dcos(yc(j)) * exp(-2*reinv*time))
                  bcxw(ii, j, k) = 0.0_CGREAL
                ENDDO
             ENDDO
        CASE (BC_TYPE_SHEAR)          ! shear bc
           print*, 'Not implemented!'
           stop
           !DO k=zb1,zb2   !0,nz   
           !DO j=0,ny
              !            bcxu(ii,j,k) = ux2*y(j)/yout 
              !            bcxv(ii,j,k) = vx2 
              !            bcxw(ii,j,k) = wx2
           !ENDDO
           !ENDDO
        END SELECT
 
!---
! front boundary
!---
        jj= 0
        j = 1
        SELECT CASE (bcy1)
        CASE (BC_TYPE_DIRICHLET)             ! dirichlet bc
           DO k=zb1,zb2
           DO i=0,nx        
              bcyu(i,jj,k) = uy1 
              bcyv(i,jj,k) = vy1 
              bcyw(i,jj,k) = wy1 
           ENDDO
           ENDDO
 
        CASE (BC_TYPE_ZERO_GRADIENT)         ! outflow bc ( zero gradient ; explicit)
           DO k=zb1,zb2
           DO i=0,nx
              bcyu(i,jj,k) = u(i,j,k)
              bcyv(i,jj,k) = v(i,j,k)
              bcyw(i,jj,k) = w(i,j,k)
           ENDDO
           ENDDO
        
        CASE (BC_TYPE_PULSATILE_INFLOW)
           DO k=zb1,zb2
           DO i=0,nx
              bcyu(i,jj,k) = uy1*sin(2.0_CGREAL*pi*freq_uy1*time)
              bcyv(i,jj,k) = vy1*sin(2.0_CGREAL*pi*freq_vy1*time)
              bcyw(i,jj,k) = wy1*sin(2.0_CGREAL*pi*freq_wy1*time)
           ENDDO
           ENDDO        
            
        CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
           DO k=zb1,zb2
           DO i=0,nx
              bcyu(i,jj,k) = u(i,j,k)
              bcyv(i,jj,k) = 0.0_CGREAL
              bcyw(i,jj,k) = w(i,j,k)
           ENDDO
           ENDDO

        CASE (BC_TYPE_PERIODIC)       ! periodic bc 
           DO k=zb1,zb2
           DO i=0,nx
              bcyu(i,jj,k) = 0.5_CGREAL*( u(i,j,k) + u(i,ny-1,k) )
              bcyv(i,jj,k) = 0.5_CGREAL*( v(i,j,k) + v(i,ny-1,k) )
              bcyw(i,jj,k) = 0.5_CGREAL*( w(i,j,k) + w(i,ny-1,k) )
           ENDDO
           ENDDO
        
        CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
           DO k=zb1,zb2
           DO i=0,nx
              ! By F.-B. Tian, use the channel flow
              bcxu(ii,j,k) = uy1
              bcxv(ii,j,k) = vy1*6.0*(z(k)-zOrigin)*(zout-z(k))/(zout-zOrigin)**2 
              bcxw(ii,j,k) = wy1
           ENDDO
           ENDDO
        case(BC_TYPE_TAYLOR_GREEN)
             DO k=zb1, zb2
                DO i=0, nx
                  bcyu(i, jj, k) = -0.5*(dcos(xc(i)) * dsin(yc(j-1)) * exp(-2*reinv*time) + dcos(xc(i)) * dsin(yc(j)) * exp(-2*reinv*time))
                  bcyv(i, jj, k) = 0.5*(dsin(xc(i)) * dcos(yc(j-1)) * exp(-2*reinv*time) + dsin(xc(i)) * dcos(yc(j)) * exp(-2*reinv*time))
                  bcyw(i, jj, k) = 0.0_CGREAL
                ENDDO
             ENDDO
        CASE (BC_TYPE_SHEAR)                ! shear bc
            print*, 'Not implemented!'
            stop
            !DO k=zb1,zb2
            !DO i=0,nx
              !            bcyu(i,jj,k) = uy1
              !            bcyv(i,jj,k) = vy1*x(i)/xout
              !            bcyw(i,jj,k) = wy1
            !ENDDO
            !ENDDO
        
        END SELECT 

!---
! back boundary
!---
        jj= 1
        j = ny-1
        SELECT CASE (bcy2)
          CASE (BC_TYPE_DIRICHLET)             ! dirichlet bc
           DO k=zb1,zb2
           DO i=0,nx
              bcyu(i,jj,k) = uy2 
              bcyv(i,jj,k) = vy2 
              bcyw(i,jj,k) = wy2 
           ENDDO
           ENDDO

        CASE (BC_TYPE_ZERO_GRADIENT)         ! outflow bc ( zero gradient ; explicit)
           DO k=zb1,zb2
           DO i=0,nx
              bcyu(i,jj,k) = u(i,j,k)
              bcyv(i,jj,k) = v(i,j,k)
              bcyw(i,jj,k) = w(i,j,k)
           ENDDO
           ENDDO

        CASE (BC_TYPE_PULSATILE_INFLOW)
           DO k=zb1,zb2
           DO i=0,nx
              bcyu(i,jj,k) = uy2*sin(2.0_CGREAL*pi*freq_uy2*time)
              bcyv(i,jj,k) = vy2*sin(2.0_CGREAL*pi*freq_vy2*time)
              bcyw(i,jj,k) = wy2*sin(2.0_CGREAL*pi*freq_wy2*time)
           ENDDO
           ENDDO

        CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
           DO k=zb1,zb2
           DO i=0,nx
              bcyu(i,jj,k) = u(i,j,k)
              bcyv(i,jj,k) = 0.0_CGREAL
              bcyw(i,jj,k) = w(i,j,k)
           ENDDO
           ENDDO

        CASE (BC_TYPE_PERIODIC)       ! periodic bc 
           DO k=zb1,zb2
           DO i=0,nx
              bcyu(i,jj,k) = 0.5_CGREAL*( u(i,j,k) + u(i,1,k) )
              bcyv(i,jj,k) = 0.5_CGREAL*( v(i,j,k) + v(i,1,k) )
              bcyw(i,jj,k) = 0.5_CGREAL*( w(i,j,k) + w(i,1,k) )
           ENDDO
           ENDDO

        CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
           DO k=zb1,zb2
           DO i=0,nx
              ! By F.-B. Tian, use the channel flow
              bcxu(ii,j,k) = uy2
              bcxv(ii,j,k) = vy2*6.0*(z(k)-zOrigin)*(zout-z(k))/(zout-zOrigin)**2 
              bcxw(ii,j,k) = wy2
           ENDDO
           ENDDO
        case(BC_TYPE_TAYLOR_GREEN)
             DO k=zb1, zb2
                DO i=0, nx
                  bcyu(i, jj, k) = -0.5*(dcos(xc(i)) * dsin(yc(j+1)) * exp(-2*reinv*time) + dcos(xc(i)) * dsin(yc(j)) * exp(-2*reinv*time))
                  bcyv(i, jj, k) = 0.5*(dsin(xc(i)) * dcos(yc(j+1)) * exp(-2*reinv*time) + dsin(xc(i)) * dcos(yc(j)) * exp(-2*reinv*time))
                  bcyw(i, jj, k) = 0.0_CGREAL
                ENDDO
             ENDDO
        CASE (BC_TYPE_SHEAR)             ! shear bc
            print*, 'Not implemented!'
            stop
            !DO k=zb1,zb2
            !DO i=0,nx
               !            bcyu(i,jj,k) = uy2
               !            bcyv(i,jj,k) = vy2*x(i)/xout 
               !            bcyw(i,jj,k) = wy2
            !ENDDO
            !ENDDO

        END SELECT 

!---
! bottom boundary
!---
        if(iProc .eq. iProc_leftmost) then    ! Only for the 1st subdomain      
        
        kk= 0
        k = 1
        SELECT CASE (bcz1)
        CASE (BC_TYPE_DIRICHLET)             ! dirichlet bc
           DO i=0,nx
           DO j=0,ny
              bczu(i,j,kk) = uz1 
              bczv(i,j,kk) = vz1 
              bczw(i,j,kk) = wz1               
           ENDDO
           ENDDO
 
        CASE (BC_TYPE_ZERO_GRADIENT)         ! outflow bc ( zero gradient ; explicit)
           DO i=0,nx
           DO j=0,ny
              bczu(i,j,kk) = u(i,j,k)
              bczv(i,j,kk) = v(i,j,k)
              bczw(i,j,kk) = w(i,j,k)              
           ENDDO
           ENDDO

        CASE (BC_TYPE_PULSATILE_INFLOW)
           DO i=0,nx
           DO j=0,ny
              bczu(i,j,kk) = uz1*cos(2.0_CGREAL*pi*freq_uz1*time)
              bczv(i,j,kk) = vz1*sin(2.0_CGREAL*pi*freq_vz1*time)
              bczw(i,j,kk) = wz1*sin(2.0_CGREAL*pi*freq_wz1*time)              
           ENDDO
           ENDDO
            
        CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
           DO i=0,nx
           DO j=0,ny
              bczu(i,j,kk) = u(i,j,k)
              bczw(i,j,kk) = 0.0_CGREAL
              bczv(i,j,kk) = v(i,j,k)              
           ENDDO
           ENDDO

        CASE (BC_TYPE_PERIODIC)       ! periodic bc 
           DO i=0,nx
           DO j=0,ny
              bczu(i,j,kk) = 0.5_CGREAL*( u(i,j,k) + u(i,j,nz-1) )
              bczv(i,j,kk) = 0.5_CGREAL*( v(i,j,k) + v(i,j,nz-1) )
              bczw(i,j,kk) = 0.5_CGREAL*( w(i,j,k) + w(i,j,nz-1) )              
           ENDDO
           ENDDO
 
        CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
           DO i=0,nx
           DO j=0,ny
              ! By F.-B. Tian, use the channel flow
              bcxu(ii,j,k) = uz1
              bcxv(ii,j,k) = vz1 
              bcxw(ii,j,k) = wz1*6.0*(y(j)-yOrigin)*(yout-y(j))/(yout-yOrigin)**2              
           ENDDO
           ENDDO
        case(BC_TYPE_TAYLOR_GREEN)
             DO i=0, nx
                DO j=0, ny
                  bczu(i, j, kk) = -dcos(xc(i)) * dsin(yc(j)) * exp(-2*reinv*time)
                  bczv(i, j, kk) = dsin(xc(i)) * dcos(yc(j)) * exp(-2*reinv*time)
                  bczw(i, j, kk) = 0.0_CGREAL
                ENDDO
             ENDDO
        CASE (BC_TYPE_SHEAR)                ! shear bc
           print*, 'Not implemented!'
           stop
           !DO i=0,nx
           !DO j=0,ny
              !            bczu(i,j,kk) = uz1
              !            bczv(i,j,kk) = vz1 
              !            bczw(i,j,kk) = wz1*y(j)/yout
              
           !ENDDO
           !ENDDO

        END SELECT 
        endif    ! if(iProc .eq. 0)
     
!---
! top boundary
!---
        if(iProc .eq. iProc_rightmost) then   ! only for the last subdomain
           
        kk= 1
        k = nz-1
        SELECT CASE (bcz2)
        CASE (BC_TYPE_DIRICHLET)             ! dirichlet bc
           DO i=0,nx
           DO j=0,ny
              bczu(i,j,kk) = uz2 
              bczv(i,j,kk) = vz2 
              bczw(i,j,kk) = wz2               
           ENDDO
           ENDDO

        CASE (BC_TYPE_ZERO_GRADIENT)         ! outflow bc ( zero gradient ; explicit)
           DO i=0,nx
           DO j=0,ny
              bczu(i,j,kk) = u(i,j,k)
              bczv(i,j,kk) = v(i,j,k)
              bczw(i,j,kk) = w(i,j,k)              
           ENDDO
           ENDDO

        CASE (BC_TYPE_PULSATILE_INFLOW)
           DO i=0,nx
           DO j=0,ny
              bczu(i,j,kk) = uz2*sin(2.0_CGREAL*pi*freq_uz2*time)
              bczv(i,j,kk) = vz2*sin(2.0_CGREAL*pi*freq_vz2*time)
              bczw(i,j,kk) = wz2*sin(2.0_CGREAL*pi*freq_wz2*time)              
           ENDDO
           ENDDO
           
        CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
           DO i=0,nx
           DO j=0,ny
              bczu(i,j,kk) = u(i,j,k)
              bczw(i,j,kk) = 0.0_CGREAL
              bczv(i,j,kk) = v(i,j,k)
           ENDDO
           ENDDO

        CASE (BC_TYPE_PERIODIC)       ! periodic bc 
           DO i=0,nx
           DO j=0,ny
              bczu(i,j,kk) = 0.5_CGREAL*( u(i,j,k) + u(i,j,1) )
              bczv(i,j,kk) = 0.5_CGREAL*( v(i,j,k) + v(i,j,1) )
              bczw(i,j,kk) = 0.5_CGREAL*( w(i,j,k) + w(i,j,1) )
           ENDDO
           ENDDO

        CASE (BC_TYPE_USER_SPECIFIED)       !  user specified 
           DO i=0,nx
           DO j=0,ny
              ! By F.-B. Tian, use the channel flow
              bcxu(ii,j,k) = uz2
              bcxv(ii,j,k) = vz2 
              bcxw(ii,j,k) = wz2*6.0*(y(j)-yOrigin)*(yout-y(j))/(yout-yOrigin)**2
           ENDDO
           ENDDO
        case(BC_TYPE_TAYLOR_GREEN)
             DO i=0, nx
                DO j=0, ny
                  bczu(i, j, kk) = -dcos(xc(i)) * dsin(yc(j)) * exp(-2*reinv*time)
                  bczv(i, j, kk) = dsin(xc(i)) * dcos(yc(j)) * exp(-2*reinv*time)
                  bczw(i, j, kk) = 0.0_CGREAL
                ENDDO
             ENDDO
        CASE (BC_TYPE_SHEAR)             ! shear bc
           print*, 'Not implemented!'
           stop
           !DO i=0,nx
           !DO j=0,ny
              !            bczu(ii,j,k) = uz2
              !            bczv(ii,j,k) = vz2 
              !            bczw(ii,j,k) = wz2*y(j)/yout
           !ENDDO
           !ENDDO
        END SELECT
        endif  !if (iProc .eq. iProc_rightmost)
      
!      ! adjust BC to satisfy global mas conservation
!      CALL enforce_global_mass_consv()

!      ! set values for outer ghost points
!      CALL set_outer_ghost_velocity()

   END SUBROUTINE set_outer_velocity_bc
!-------------------------------------------------------------------------

! compute mass flux at all boundaries and adjust outflow BC so as to satisfy
! global mass conservation.

   SUBROUTINE enforce_global_mass_consv()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE GCM_arrays
    USE MPI_module
    
    IMPLICIT NONE

    INTEGER             :: i,j,k,n
    REAL(KIND=CGREAL)   :: massflux,massflux_i,correction_vel,fluxBody
    REAL(KIND=CGREAL)   :: flag,rib1, rib2
    REAL(KIND=CGREAL)   :: area_i

    massflux   = 0.0_CGREAL
    massflux_i = 0.0_CGREAL
    
    IF ( bcx1 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcx2 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcy1 == BC_TYPE_ZERO_GRADIENT .OR. & 
         bcy2 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcz1 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcz2 == BC_TYPE_ZERO_GRADIENT  )   THEN

       !DO j=1,ny-1
       !   rib1 = 1.0_CGREAL-REAL(iblank(1   ,j),KIND=CGREAL)
       !   rib2 = 1.0_CGREAL-REAL(iblank(nx-1,j),KIND=CGREAL)
       !   massflux = massflux + (-bcxu(0,j)*rib1 + bcxu(1,j)*rib2)*dy(j)
       !ENDDO
       !DO i=1,nx-1
       !   rib1 = 1.0_CGREAL-REAL(iblank(i,   1),KIND=CGREAL)
       !   rib2 = 1.0_CGREAL-REAL(iblank(j,ny-1),KIND=CGREAL)
       !   massflux = massflux + (-bcyv(i,0)*rib1 + bcyv(i,1)*rib2)*dx(i)
       !ENDDO
       
       ! Note that face_vel needs to be updated outside of this subroutine.
       !CALL face_vel()

       ! compute the overall mass flux.  Ghost-cells and blank cells are
       ! excluded.

       DO k=zc_start,zc_end  !1,nz-1
       DO j=1,ny-1
       DO i=1,nx-1
         flag = (1.0_CGREAL-REAL(iblank       (i,j,k),KIND=CGREAL))  &
              * (1.0_CGREAL-REAL(ghostcellMark(i,j,k),KIND=CGREAL))

         massflux_i = massflux_i +  flag *                      &
                  ( (-face_u(i,j,k) + face_u(i+1,j,k))*dy(j)*dz(k) +  &
                    (-face_v(i,j,k) + face_v(i,j+1,k))*dx(i)*dz(k) +  &
                    (-face_w(i,j,k) + face_w(i,j,k+1))*dx(i)*dy(j)    &
                  )
       ENDDO ! i
       ENDDO ! j
       ENDDO ! k
       !print*, 'Massflux_i and massflux of iProc ', iProc, 'are ', massflux_i, massflux

       IF(mix_GC_form == 1) THEN
         CALL correct_ghost_massflux(massflux_i)
         !print*,iProc, 'massflux = ', massflux_i
       ENDIF

       ! Sum up massflux from all subdomains, and all processors will get the final sum
       CALL MPI_ALLREDUCE(massflux_i, massflux, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, flow_comm, ierr)       
       !print*, 'Massflux_i and massflux of iProc ', iProc, 'are ', massflux_i, massflux

       ! compute areas of external boundaries
       area_left    = 0.0_CGREAL
       area_right   = 0.0_CGREAL
       area_bot     = 0.0_CGREAL
       area_top     = 0.0_CGREAL
       area_front   = 0.0_CGREAL
       area_back    = 0.0_CGREAL
       !outflow_area= 0.0_CGREAL
       area_i       = 0.0_CGREAL

       DO k=zc_start,zc_end   !1,nz-1
       DO j=1,ny-1
          area_left  = area_left  + dy(j)*dz(k)*(1.0_CGREAL-iblank(1,   j,k) )
          area_right = area_right + dy(j)*dz(k)*(1.0_CGREAL-iblank(nx-1,j,k) )
       ENDDO
       ENDDO
       DO k=zc_start,zc_end   !1,nz-1
       DO i=1,nx-1
          area_bot   = area_bot   + dx(i)*dz(k)*(1.0_CGREAL-iblank(i,1   ,k) )
          area_top   = area_top   + dx(i)*dz(k)*(1.0_CGREAL-iblank(i,ny-1,k) )
       ENDDO
       ENDDO

       IF(iProc==iProc_leftmost) then
         DO j=1,ny-1
         DO i=1,nx-1
            area_back  = area_back  + dx(i)*dy(j)*(1.0_CGREAL-iblank(i,j   ,1) )
         ENDDO
         ENDDO
       ENDIF

       IF(iProc==iProc_rightmost) then
         DO j=1,ny-1
         DO i=1,nx-1
            area_front = area_front + dx(i)*dy(j)*(1.0_CGREAL-iblank(i,j,nz-1) )
         ENDDO
         ENDDO
       ENDIF

       !area_back and area_front remains zero for middle processors
       IF (bcx1 == BC_TYPE_ZERO_GRADIENT) area_i = area_i + area_left
       IF (bcx2 == BC_TYPE_ZERO_GRADIENT) area_i = area_i + area_right
       IF (bcy1 == BC_TYPE_ZERO_GRADIENT) area_i = area_i + area_bot
       IF (bcy2 == BC_TYPE_ZERO_GRADIENT) area_i = area_i + area_top
       IF (bcz1 == BC_TYPE_ZERO_GRADIENT) area_i = area_i + area_back                
       IF (bcz2 == BC_TYPE_ZERO_GRADIENT) area_i = area_i + area_front    

       CALL MPI_ALLREDUCE(area_i, outflow_area, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, FLOW_COMM, ierr)
       correction_vel =-massflux/outflow_area  

       !IF ( MOD(ntime,nmonitor) == 0 .and. iProc==PROC_M) THEN
       !  PRINT*,' SET_BC:massflux       = ',massflux
       !  PRINT*,' SET_BC:outflow_area   = ',outflow_area
       !  PRINT*,' SET_BC:correction_vel = ',correction_vel
       !END IF ! ntime

       DO k=zc_start-1,zc_end+1   !0,nz
       DO j=0,ny
         IF (bcx1 == BC_TYPE_ZERO_GRADIENT) bcxu(0,j,k) = bcxu(0,j,k) - correction_vel
         IF (bcx2 == BC_TYPE_ZERO_GRADIENT) bcxu(1,j,k) = bcxu(1,j,k) + correction_vel
       ENDDO ! j
       ENDDO ! k

       DO k=zc_start-1,zc_end+1   !0,nz
       DO i=0,nx
         IF (bcy1 == BC_TYPE_ZERO_GRADIENT) bcyv(i,0,k) = bcyv(i,0,k) - correction_vel
         IF (bcy2 == BC_TYPE_ZERO_GRADIENT) bcyv(i,1,k) = bcyv(i,1,k) + correction_vel
       ENDDO ! i
       ENDDO ! k

       IF(iProc==iProc_leftmost) then
         DO j=0,ny
         DO i=0,nx
           IF (bcz1 == BC_TYPE_ZERO_GRADIENT) bczw(i,j,0) = bczw(i,j,0) - correction_vel
         ENDDO ! i
         ENDDO ! j
       ENDIF
       IF(iProc==iProc_rightmost) then
         DO j=0,ny
         DO i=0,nx
           IF (bcz2 == BC_TYPE_ZERO_GRADIENT) bczw(i,j,1) = bczw(i,j,1) + correction_vel
         ENDDO ! i
         ENDDO ! j
       ENDIF

       DO k = zc_start,zc_end   !1,nz-1
       DO j = 1,ny-1
          IF (bcx1 == BC_TYPE_ZERO_GRADIENT) face_u(1, j,k) = face_u(1, j,k) - correction_vel
          IF (bcx2 == BC_TYPE_ZERO_GRADIENT) face_u(nx,j,k) = face_u(nx,j,k) + correction_vel
       ENDDO
       ENDDO 

       DO k = zc_start,zc_end   !1,nz-1
       DO i = 1,nx-1
          IF (bcy1 == BC_TYPE_ZERO_GRADIENT) face_v(i,1,k ) = face_v(i,1 ,k) - correction_vel
          IF (bcy2 == BC_TYPE_ZERO_GRADIENT) face_v(i,ny,k) = face_v(i,ny,k) + correction_vel
       ENDDO
       ENDDO

       IF(iProc==iProc_leftmost) then
         DO j = 1,ny-1
         DO i = 1,nx-1
            IF (bcz1 == BC_TYPE_ZERO_GRADIENT) face_w(i,j,1 ) = face_w(i,j,1)  - correction_vel
         ENDDO
         ENDDO
       ENDIF
       IF(iProc==iProc_rightmost) then
         DO j = 1,ny-1
         DO i = 1,nx-1
            IF (bcz2 == BC_TYPE_ZERO_GRADIENT) face_w(i,j,nz) = face_w(i,j,nz) + correction_vel
         ENDDO
         ENDDO
       ENDIF

    ENDIF ! bcx1

   END SUBROUTINE enforce_global_mass_consv

!-------------------------------------------------------------------------

   SUBROUTINE set_outer_ghost_velocity()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE blasius_profile
    USE MPI_module

    IMPLICIT NONE

    INTEGER             :: i,j,k,ii,jj,kk

!--------
! outer ghost boundary conditions for velocities
!--------

      DO k=zb1,zb2   !0,nz
      DO j=0,ny

! left boundary
        ii= 0
        i = 0
        u(i,j,k) = 2.0_CGREAL*bcxu(ii,j,k) - u(i+1,j,k)
        v(i,j,k) = 2.0_CGREAL*bcxv(ii,j,k) - v(i+1,j,k)
        w(i,j,k) = 2.0_CGREAL*bcxw(ii,j,k) - w(i+1,j,k)

! right boundary
        ii= 1
        i = nx
        u(i,j,k) = 2.0_CGREAL*bcxu(ii,j,k) - u(i-1,j,k)
        v(i,j,k) = 2.0_CGREAL*bcxv(ii,j,k) - v(i-1,j,k)
        w(i,j,k) = 2.0_CGREAL*bcxw(ii,j,k) - w(i-1,j,k)

      ENDDO ! j
      ENDDO ! k

      DO k=zb1,zb2   !0,nz
      DO i=0,nx

! front boundary
        jj= 0
        j = 0
        u(i,j,k) = 2.0_CGREAL*bcyu(i,jj,k) - u(i,j+1,k)
        v(i,j,k) = 2.0_CGREAL*bcyv(i,jj,k) - v(i,j+1,k)
        w(i,j,k) = 2.0_CGREAL*bcyw(i,jj,k) - w(i,j+1,k)

! back boundary
        jj= 1
        j = ny
        u(i,j,k) = 2.0_CGREAL*bcyu(i,jj,k) - u(i,j-1,k)
        v(i,j,k) = 2.0_CGREAL*bcyv(i,jj,k) - v(i,j-1,k)
        w(i,j,k) = 2.0_CGREAL*bcyw(i,jj,k) - w(i,j-1,k)
      ENDDO ! i
      ENDDO ! k

! bottom boundary
      if(iProc .eq. iProc_leftmost) then      
      DO i=0,nx
      DO j=0,ny
        kk= 0
        k = 0
        u(i,j,k) = 2.0_CGREAL*bczu(i,j,kk) - u(i,j,k+1)
        v(i,j,k) = 2.0_CGREAL*bczv(i,j,kk) - v(i,j,k+1)
        w(i,j,k) = 2.0_CGREAL*bczw(i,j,kk) - w(i,j,k+1)
      ENDDO ! i
      ENDDO ! j
      endif
   
! top boundary
      if(iProc .eq. iProc_rightmost) then 
      DO i=0,nx
      DO j=0,ny
        kk= 1
        k = nz
        u(i,j,k) = 2.0_CGREAL*bczu(i,j,kk) - u(i,j,k-1)
        v(i,j,k) = 2.0_CGREAL*bczv(i,j,kk) - v(i,j,k-1)
        w(i,j,k) = 2.0_CGREAL*bczw(i,j,kk) - w(i,j,k-1)
      ENDDO ! i
      ENDDO ! j
      endif
   
   END SUBROUTINE set_outer_ghost_velocity
!-------------------------------------------------------------------------
! This subroutine sets up the B.C. for the pressure solver at the outer boundary.
! Note that the size of the actual matrix may vary during the multigrid iterations.
! The actual size is defined as (iin, jin).
! Also note that the B.C. may be homogeneous if the equation is solved on the 
! coarse grids.

   SUBROUTINE set_outer_pressure_bc(pres,nx_mg,ny_mg,zb1g,zb2g,nlev)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE mg_module
    USE MPI_module
    
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nx_mg,ny_mg,zb1g,zb2g,nlev
    REAL(KIND=CGREAL),DIMENSION(0:nx_mg+1,0:ny_mg+1,zb1g:zb2g),INTENT (INOUT) ::pres

   
    REAL(KIND=CGREAL):: ppx1,ppx2,ppy1,ppy2,ppz1,ppz2
    integer :: i,j,k, iHomo,zk1,zk2
!******************************************

    !print*,'nlev', nlev,nx_mg,ny_mg,zb1g,zb2g

    ppx1=pppx1      !! pppx1 is read from input file
    ppx2=pppx2
    ppy1=pppy1
    ppy2=pppy2
    ppz1=pppz1
    ppz2=pppz2

    iHomo = 0

    ! For coarse grids in multigrid iterations, set homogeneous B.C., iHomo=1
    ! Note that use of original grid variables such as dxc, dyc, etc, and
    ! the pgradx, pgrady, etc, are not compatible in array size with coarse grid variables.
    ! However, this is still fine since they only take effect when iHomo=0.

    if(nlev /= 1)  iHomo = 1 

    IF(iHomo == 1) THEN
       ppx1 = 0.0_CGReAL
       ppx2 = 0.0_CGReAL
       ppy1 = 0.0_CGReAL
       ppy2 = 0.0_CGReAL
       ppz1 = 0.0_CGReAL
       ppz2 = 0.0_CGReAL
    ENDIF

    if(it_solver_type == IT_SOLVER_TYPE_MG) then
       !iin = mgrid_I(nlev)
       !jjn = mgrid_J(nlev)
       zk1 = zc_start_mg(nlev)
       zk2 = zc_end_mg  (nlev)
    else
       zk1 = zc_start
       zk2 = zc_end
    endif

    !left boundary
    select case (pbcx1)
    case (PBC_DIRICHLET)
       do k=zk1,zk2  !1,kkn
       do j=1,ny_mg-1
          bcxp(0,j,k) = ppx1  
       enddo
       enddo
    case (PBC_NEUMANN) 
       do k=zk1,zk2  !1,kkn
       do j=1,ny_mg-1
          bcxp(0,j,k) = pres(1,j,k)   - pgradx1(j,k)*dxc(1)/2.0_CGREAL*(1-iHomo)
       enddo
       enddo
    end select

    !right boundary
    select case (pbcx2)                              
    case (PBC_DIRICHLET)
       do k=zk1,zk2  !1,kkn
       do j=1,ny_mg-1
          bcxp(1,j,k) = ppx2
       enddo
       enddo
    case (PBC_NEUMANN)
       do k=zk1,zk2  !1,kkn
       do j=1,ny_mg-1
          bcxp(1,j,k) = pres(nx_mg-1,j,k) + pgradx2(j,k)*dxc(nx_mg)/2.0_CGREAL*(1-iHomo)
       enddo
       enddo
    end select

    !front boundary
    select case (pbcy1)
    case (PBC_DIRICHLET)
       do k=zk1,zk2  !1,kkn
       do i=1,nx_mg-1
          bcyp(i,0,k) = ppy1
       enddo
       enddo
    case (PBC_NEUMANN)
       do k=zk1,zk2  !1,kkn
       do i=1,nx_mg-1
          bcyp(i,0,k) = pres(i,1,k)   - pgrady1(i,k)*dyc(1)/2.0_CGREAL*(1-iHomo)
       enddo
       enddo
    end select
       
    !back boundary
    select case (pbcy2)
    case (PBC_DIRICHLET)
       do k=zk1,zk2  !1,kkn
       do i=1,nx_mg-1
          bcyp(i,1,k) = ppy2
       enddo
       enddo
    case (PBC_NEUMANN)
       do k=zk1,zk2  !1,kkn
       do i=1,nx_mg-1
          bcyp(i,1,k) = pres(i,ny_mg-1,k) + pgrady2(i,k)*dyc(ny_mg)/2.0_CGREAL*(1-iHomo)
       enddo
       enddo
    end select

    !bottom boundary
    if(iProc==iProc_leftmost) then
       select case (pbcz1)
       case (PBC_DIRICHLET)
          do j=1,ny_mg-1
          do i=1,nx_mg-1
             bczp(i,j,0) = ppz1
          enddo
          enddo
       case (PBC_NEUMANN)
          do j=1,ny_mg-1
          do i=1,nx_mg-1
             bczp(i,j,0) = pres(i,j,1)   - pgradz1(i,j)*dzc(1)/2.0_CGREAL*(1-iHomo) 
          enddo
          enddo
       end select
    endif

    !top boundary
    if(iProc==iProc_rightmost) then       
       select case (pbcz2)
       case (PBC_DIRICHLET)
          do j=1,ny_mg-1
          do i=1,nx_mg-1
             bczp(i,j,1) = ppz2
          enddo
          enddo
       case (PBC_NEUMANN)
          do j=1,ny_mg-1
          do i=1,nx_mg-1
             bczp(i,j,1) = pres(i,j,zk2) + pgradz2(i,j)*dzc(nz)/2.0_CGREAL*(1-iHomo)
          enddo
          enddo
       end select
    endif

    ! Periodic B.C. for pressure
    do k=zk1,zk2  !1,kkn
    do j=1,ny_mg-1
       IF (bcx1  .EQ.  BC_TYPE_PERIODIC .AND. &
           bcx2  .EQ.  BC_TYPE_PERIODIC)      THEN
          bcxp(0,   j,k) = ( pres(1,j,k)+pres(nx_mg-1,j,k))/2.0_CGREAL
          bcxp(1,   j,k) =  bcxp(0,j,k)
       ENDIF
    enddo
    enddo

    do k=zk1,zk2  !1,kkn
    do i=1,nx_mg-1
       IF (bcy1  .EQ.  BC_TYPE_PERIODIC .AND. &
           bcy2  .EQ.  BC_TYPE_PERIODIC)      THEN
          bcyp(i,0,k   ) = ( pres(i,1,k)+pres(i,ny_mg-1,k))/2.0_CGREAL
          bcyp(i,1,k   ) =  bcyp(i,0,k)
       ENDIF
    enddo
    enddo

    !Periodic condition in z not implemented yet
    do i=1,nx_mg-1
    do j=1,ny_mg-1
       IF (bcz1  .EQ.  BC_TYPE_PERIODIC .AND. &
           bcz2  .EQ.  BC_TYPE_PERIODIC)      THEN
          print*,'Periodic condition in z not implemented yet.'
          stop

          bczp(i,j,0   ) = ( pres(i,j,1)+pres(i,j,nz-1))/2.0_CGREAL
          bczp(i,j,1   ) =  bczp(i,j,0)
       ENDIF
    enddo
    enddo

  END SUBROUTINE set_outer_pressure_bc
!-------------------------------------------------------------------------
  SUBROUTINE set_outer_ghost_pressure(pres,nx_mg,ny_mg,zb1g,zb2g,nlev)
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE mg_module
    USE MPI_module

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nx_mg,ny_mg,zb1g,zb2g,nlev
    REAL(KIND=CGREAL), DIMENSION(0:nx_mg+1,0:ny_mg+1,zb1g:zb2g), INTENT(INOUT) :: pres

    INTEGER             :: i,j,k, iin,jjn,kkn,zk1,zk2
    INTEGER             :: ii, jj, kk

    if(it_solver_type == IT_SOLVER_TYPE_MG) then
       !iin = mgrid_I(nlev)
       !jjn = mgrid_J(nlev)
       zk1 = zc_start_mg(nlev)
       zk2 = zc_end_mg  (nlev)
    else
       zk1 = zc_start
       zk2 = zc_end
    endif

    !Setting outer ghost boundary conditions. The grid size may not be
    ! the baseline size, as it may change in the multigrid solver.

    ! This subroutine needs work once we start to implement multigrid solver.
    !pres(0,    1:jjn,1:kkn)= 2.0_CGREAL*bcxp(0, 1:jjn,1:kkn)-pres(1,  1:jjn,1:kkn)
    !pres(iin+1,1:jjn,1:kkn)= 2.0_CGREAL*bcxp(1, 1:jjn,1:kkn)-pres(iin,1:jjn,1:kkn)

    !pres(1:iin,    0,1:kkn)= 2.0_CGREAL*bcyp(1:iin, 0,1:kkn)-pres(1:iin,   1,1:kkn)
    !pres(1:iin,jjn+1,1:kkn)= 2.0_CGREAL*bcyp(1:iin, 1,1:kkn)-pres(1:iin, jjn,1:kkn)

    !pres(1:iin,1:jjn,0    )= 2.0_CGREAL*bczp(1:iin, 1:jjn,0)-pres(1:iin,1:jjn,1   )
    !pres(1:iin,1:jjn,kkn+1)= 2.0_CGREAL*bczp(1:iin, 1:jjn,1)-pres(1:iin,1:jjn, kkn)

    DO k=zk1,zk2   !0,nz
    DO j=0,ny_mg
       ! left boundary
       ii= 0
       i = 0
       pres(i,j,k) = 2.0_CGREAL*bcxp(ii,j,k) - pres(i+1,j,k)

       ! right boundary
       ii= 1
       i = nx
       pres(i,j,k) = 2.0_CGREAL*bcxp(ii,j,k) - pres(i-1,j,k)
    ENDDO ! j
    ENDDO ! k

    DO k=zk1,zk2   !0,nz
    DO i=0,nx_mg
       !front boundary
       jj= 0
       j = 0
       pres(i,j,k) = 2.0_CGREAL*bcyp(i,jj,k) - pres(i,j+1,k)

       ! back boundary
       jj= 1
       j = ny
       pres(i,j,k) = 2.0_CGREAL*bcyp(i,jj,k) - pres(i,j-1,k)
    ENDDO ! i
    ENDDO ! k

    ! bottom boundary
    if(iProc .eq. iProc_leftmost) then  
      DO j=0,ny_mg
      DO i=0,nx_mg
        kk= 0
        k = 0
        pres(i,j,k) = 2.0_CGREAL*bczp(i,j,kk) - pres(i,j,k+1)
      ENDDO ! i
      ENDDO ! j
   endif
   
   ! top boundary
   if(iProc .eq. iProc_rightmost) then 
      DO j=0,ny_mg
      DO i=0,nx_mg
        kk= 1
        k = zk2+1
        pres(i,j,k) = 2.0_CGREAL*bczp(i,j,kk) - pres(i,j,k-1)
      ENDDO ! i
      ENDDO ! j
   endif

  END SUBROUTINE set_outer_ghost_pressure
!-------------------------------------------------------------------------
!
! Outer boundary pressure update using a second order scheme. Written by Fangbao
! Only for the homegeneous Neumann condition, dp/dn = 0
!
   SUBROUTINE set_outer_ghost_pressure_2nd(pres)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE MPI_module

    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,0:ny+1,zb1:zb2), INTENT(INOUT) :: pres
    INTEGER             :: i,j,k,k1,k2

    k1=zc_start
    k2=zc_end

    ! outer ghost boundary conditions
    pres(0, 1:ny-1,k1:k2)= ( 4.0_CGREAL*pres(1,   1:ny-1,k1:k2) - pres(2,   1:ny-1,k1:k2) )/3.0_CGREAL
    pres(nx,1:ny-1,k1:k2)= ( 4.0_CGREAL*pres(nx-1,1:ny-1,k1:k2) - pres(nx-2,1:ny-1,k1:k2) )/3.0_CGREAL

    pres(1:nx-1, 0,k1:k2)= ( 4.0_CGREAL*pres(1:nx-1,   1,k1:k2) - pres(1:nx-1,   2,k1:k2) )/3.0_CGREAL
    pres(1:nx-1,ny,k1:k2)= ( 4.0_CGREAL*pres(1:nx-1,ny-1,k1:k2) - pres(1:nx-1,ny-2,k1:k2) )/3.0_CGREAL

    if(iProc .eq. iProc_leftmost)  &
    pres(1:nx-1,1:ny-1,0 )= ( 4.0_CGREAL*pres(1:nx-1,1:ny-1,1   ) - pres(1:nx-1,1:ny-1,2   ) )/3.0_CGREAL

    if(iProc .eq. iProc_rightmost) &
    pres(1:nx-1,1:ny-1,nz)= ( 4.0_CGREAL*pres(1:nx-1,1:ny-1,nz-1) - pres(1:nx-1,1:ny-1,nz-2) )/3.0_CGREAL

   END SUBROUTINE set_outer_ghost_pressure_2nd
!-------------------------------------------------------------------------


