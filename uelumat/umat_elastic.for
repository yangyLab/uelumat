        SUBROUTINE UMAT1(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1  RPL,DDSDDT,DRPLDE,DRPLDT,
     2  STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3  NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4  CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
c
        INCLUDE 'ABA_PARAM.INC'
c
        CHARACTER*80 CMNAME
        DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1  DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2  STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3  PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4  JSTEP(4)
        double precision e,v
c
        ddsdde=0.d0
      
c     Build stiffness matrix
        E=PROPS(1)
        v=PROPS(2)
c       ddsdde(1,1)=ONE*E*(1/(1-v*v))
c       ddsdde(2,2)=ONE*E*(1/(1-v*v))
c       ddsdde(1,2)=EMU*E*(1/(1-v*v))
c       ddsdde(2,1)=EMU*E*(1/(1-v*v))
c       ddsdde(1,3)=ZERO
c       ddsdde(2,3)=ZERO
c       ddsdde(3,1)=ZERO
c       ddsdde(3,2)=ZERO
c       ddsdde(3,3)=(ONE-EMU)/TWO*E*(1/(1-EMU*EMU))
c      plane strain
c       ddsdde(1,1)=E*(1.0d0-v)/((1.0d0+v)*(1.0d0-2.0d0*v))
c       ddsdde(1,2)=E*v/((1.0d0+v)*(1.0d0-2.0d0*v))
c       ddsdde(1,3)=E*v/((1.0d0+v)*(1.0d0-2.0d0*v))
c       ddsdde(1,4)=0.0d0
c       ddsdde(2,1)=E*v/((1.0d0+v)*(1.0d0-2.0d0*v))
c       ddsdde(2,2)=E*(1.0d0-v)/((1.0d0+v)*(1.0d0-2.0d0*v))
c       ddsdde(2,3)=E*v/((1.0d0+v)*(1.0d0-2.0d0*v))
c       ddsdde(2,4)=0.0d0
c       ddsdde(3,1)=E*v/((1.0d0+v)*(1.0d0-2.0d0*v))
c       ddsdde(3,2)=E*v/((1.0d0+v)*(1.0d0-2.0d0*v))
c       ddsdde(3,3)=E*v/((1.0d0+v)*(1.0d0-2.0d0*v))
c       ddsdde(3,4)=0.0d0
c       ddsdde(4,1)=0.0d0
c       ddsdde(4,2)=0.0d0
c       ddsdde(4,3)=0.0d0
c       ddsdde(4,4)=E*(1.0d0-2.0d0*v)/(2.0d0*((1.0d0+v)*(1.0d0-2.0d0*v)))
c      plane stress
       ddsdde(1,1)=E/((1.0d0-v)*(1.0d0+v))
       ddsdde(1,2)=E*v/((1.0d0-v)*(1.0d0+v))
       ddsdde(1,3)=0.0d0
       ddsdde(1,4)=0.0d0
       ddsdde(2,1)=E*v/((1.0d0-v)*(1.0d0+v))
       ddsdde(2,2)=E/((1.0d0-v)*(1.0d0+v))
       ddsdde(2,3)=0.0d0
       ddsdde(2,4)=0.0d0
       ddsdde(3,1)=0.0d0
       ddsdde(3,2)=0.0d0
       ddsdde(3,3)=0.0d0
       ddsdde(3,4)=0.0d0
       ddsdde(4,1)=0.0d0
       ddsdde(4,2)=0.0d0
       ddsdde(4,3)=0.0d0
       ddsdde(4,4)=E/(2.0d0*(1.0d0+v))
c     
c     recover stress and strain from statev array
c
       do k1=1,ntens
         stress(k1)=statev(k1)
         stran(k1) =statev(k1+ntens)
       end do
c     
       do k1=1,ntens
          do k2=1,ntens
           stress(k2)=stress(k2)+ddsdde(k2,k1)*dstran(k1)
          end do
       end do  
      write(300,*) "stress"
       do k=1,ntens
          write(300,*) stress(k)
       end do
c       do i=1,ntens
c         stran(k1)=stran(k1)+dstran(k1)
c       end do
c
c     store stress and strain in state variable array
c
       do k1=1,ntens
         statev(k1)=stress(k1)
         statev(k1+ntens)=stran(k1) 
       end do
c      
       return 
       end 