c    Solving mean field equations
c
c Nestor Caticha.
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      common/a3/delta,beta,amu1,amu2,a,b,BB1,BB2,AA1,AA2,r1

      open(31,file='absrm.dat')
      open(33,file='specificheat.dat')
      open(34,file='absrmmaglargura.dat')
      open(35,file='absmagdelta.dat')


c    parameters
     
      delta=0.0d0      
      a=.5d0*(1.d0+delta)
      b=.5d0*(1.d0-delta)
      tol=1.d-5
      nmax=2000

      do beta=1.d0,40.d0,1.d-1
      amu1=.9 d0
      r1  =.95d0
         nit=0
c ifla=0 does not converge
 1       ifla=0
         if(nit.lt.nmax)then
c if less than max steps, calculates <h> and <|h|>
            call aintegra(denomint,anumint, abscos)
c            write(*,*)' >>>>>  ', denomint,anumint
            amu2=anumint/denomint 
            r2=abscos/denomint
c escreve estes valores intermediarios e mostra as curvas (fazer grafico sem linhas,
c so simbolos) que vao se aproximando da curva de magnetizacao  m e de r
            WRITE(33,*)beta*.25d0,amu2,r2
c            WRITE(33,*)1.d0/beta,amu2,r2
c verifica se ja convergiu
            if(abs(amu2-amu1).gt.tol)then
c nao convergiu , adota os valores calculados como os novos valores 
               amu1=amu2
               r1=r2
               nit=nit+1
            else
c convergiu , escreve os valores de beta mag, abs (h) e numero de iteracoes ate convergencia
c               write(31,*)1.d0/beta,amu2,r2,r1,nit
               write(31,*)beta*.25d0,amu2,r2,r1,nit
               ifla=1
               call flutua(cossq)
               cosdq=cossq-amu2**2
c              write(34,*)1.d0/beta,cosdq
              write(34,*)beta,cosdq
c ,cossq
               if((beta.lt.15.01d0).and.(beta.gt.14.99d0))then
                  write(35,*)delta,cosdq
                  write(*,*)delta,cosdq
               endif
            endif
         else
               write(*,*)1.d0/beta,amu2,amu1, 'nao convergiu  '
               ifla=1
         endif
         if(ifla.eq.0)goto 1
      enddo
c         write(31,*)'&'
c      enddo
 100  format(2e15.4)
      stop
      end
C ROTINA DE INTEGRAÇAO EM duas DIMESOES
      subroutine  aintegra(denomint,anumint,abscos)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      common/a3/delta,beta,amu1,amu2,a,b,BB1,BB2,AA1,AA2,r1
      pis2=dacos(-1.d0)*0.5d0
      xmax=1.d0
      xmin=0.d0
      aden=0.d0
      anum=0.d0
      absacc=0.d0
      eps=1.d0/4000.d0
      ammbr=beta*(.5*(1.+delta)*amu1-.5*(1.-delta)*r1)
      ampbr=beta*(.5*(1.+delta)*amu1+.5*(1.-delta)*r1)
      dx=(xmax-xmin)*eps
      do x=xmin,xmax,dx
c         eem=exp(amu1*beta*delta*x)-exp(amu1*beta*x)
c         eep=exp(amu1*beta*delta*x)+exp(amu1*beta*x)

         thetap=x*pis2
         thetam=thetap+pis2

         eep=dsin(thetap)**3*dexp(ammbr*dcos(thetap))
         eem=dsin(thetam)**3*dexp(ampbr*dcos(thetam))
         aden=aden+eep+eem
         anum=anum+eep*dcos(thetap)+eem*dcos(thetam)
         absacc=absacc+eep*dcos(thetap)-eem*dcos(thetam)
c     eem*dabs(dcos(thetam))
      enddo
      denomint=aden
      anumint =anum
      abscos  =absacc
c      anumint =sinh(beta*amu1)
c      denomint=cosh(beta*amu1)
      return
      end


c      SUBROUTINE parametros
c        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cc      IMPLICIT DOUBLE PRECISION (a-h,o-z)
c      common/a1/beta,beta2,deltac


c      beta=.1d0
c      delta=.2d0

c      return
c      end
c$$$      
c calculo das flutuacoes do cosseno
C ROTINA DE INTEGRAÇAO  EM duas DIMESOES
      subroutine flutua (cossq)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      common/a3/delta,beta,amu1,amu2,a,b,BB1,BB2,AA1,AA2,r1
      pis2=dacos(-1.d0)*0.5d0
      xmax=1.d0
      xmin=0.d0
      aden=0.d0
      anum=0.d0
      eps=1.d0/2000.d0
      ammbr=beta*(.5*(1.+delta)*amu1-.5*(1.-delta)*r1)
      ampbr=beta*(.5*(1.+delta)*amu1+.5*(1.-delta)*r1)
      dx=(xmax-xmin)*eps
      do x=xmin,xmax,dx
c         eem=exp(amu1*beta*delta*x)-exp(amu1*beta*x)
c         eep=exp(amu1*beta*delta*x)+exp(amu1*beta*x)

         thetap=x*pis2
         thetam=thetap+pis2

         eep=dsin(thetap)**3*dexp(ammbr*dcos(thetap))
         eem=dsin(thetam)**3*dexp(ampbr*dcos(thetam))
         aden=aden+eep+eem
         anum=anum+eep*dcos(thetap)**2+eem*dcos(thetam)**2
      enddo
c      denomint=aden*dx
c      anumint =anum*dx
      cossq=anum/aden
c      anumint =sinh(beta*amu1)
c      denomint=cosh(beta*amu1)
      return
      end
