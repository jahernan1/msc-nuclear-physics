c----------------------------------------------------------------------
c     Having computed ground-state densities in a self-consistent
c     relativistic mean-field model, we compute proton, neutron, 
c     charge, and weak-charge densities and form factors. Moreover, 
c     we also compute the above densities and form factors using 
c     2-parameter Symmetric Fermi and Helm representations.
c                                               January 30, 2012.
c----------------------------------------------------------------------
      Program NuclearFF
c 
      implicit real*8 (a-h,j-z) 
c
c     call Reader       !read p/n densities and compute FFs
c      call ChRhoFF      !compute electric charge density and FF
c      call WkRhoFF      !compute weak charge density and FF
c      call HelmSFermi   !compute equivalent SFermi and Helm
c      call Writer       !write densities and FFs
c
c      stop
      open(unit = 111, file = 'GPN.txt')
c      double precision Gp(1000), Gn(1000)
      do i = 0, 8000
            q = 0.01 * i
            print*, q
            call NucleonFF(q, GEp, GEn)
c            Gp(i) = GEp
c            Gn(i) = GEn
            write(111, *) q,GEp,GEn
111         format(A, F1.6)
      end do
      close(111)


      end
c----------------------------------------------------------------------
c     Define r-space and q-space grids. Then, read from file the
c     previously computed relavistic mean-field proton and neutron 
c     densities and compute their corresponding form factors.
c
      subroutine Reader
c
      implicit real*8 (a-h,j-z)  
      parameter(iPts=1000)
      character*10 header
      common/QGrid/qmin,qmax,dq,iQmax
      common/RGrid/rmin,rmax,dr,iRmax
      common/pDensity/Rhop(0:iPts),FFp(0:iPts)
      common/nDensity/Rhon(0:iPts),FFn(0:iPts)
c
c     Define r-space grid:
      dr=0.05000d0                    !step in r (in fm)
      rmin=0.000d0                    !minimum r (in fm)
      rmax=15.00d0                    !maximum r (in fm)
      iRmax=nint((rmax-rmin)/dr)      !number of points
c
c     Define q-space grid:
      iQmax=iRmax                      !number of points
      qmin=0.000d0                     !minimum q (in 1/fm)
      qmax=10.00d0                     !maximum q (in 1/fm)
      dq=(qmax-qmin)/iQmax             !step in q (in 1/fm)
c
c     Read the proton/neutron mean-field densities:
      open(unit=3,file='Ar40Dens.out',status='old')
      read(3,'(a)') header             !read header
      read(3,'(a)') header             !read header
      do ir=0,iRmax                    !loop over RGrid
       read(3,*) r,Rsn,Rvn,Rsp,Rvp,Rch !read from file  
       Rhop(ir)=Rvp                    !p-density 
       Rhon(ir)=Rvn                    !n-density 
      end do                           !close ir-loop  
c     Compute corresponding form factors:
      call FormFactor(Rhop,FFp)        !p form factor 
      call FormFactor(Rhon,FFn)        !n form factor 
c
      return
      end
c----------------------------------------------------------------------
c     Having read proton and neutron densities and computed 
c     their corresponding form factors, compute the charge
c     density and its corresponding form factor.
c
      subroutine ChRhoFF
c
      implicit real*8 (a-h,j-z)  
      parameter(iPts=1000)
      common/QGrid/qmin,qmax,dq,iQmax
      common/pDensity/Rhop(0:iPts),FFp(0:iPts)
      common/nDensity/Rhon(0:iPts),FFn(0:iPts)
      common/cDensity/Rhoc(0:iPts),FFc(0:iPts)
c
c     Compute nuclear charge form factor:
      do iq=0,iQmax                    !loop over QGrid
       q=qmin+iq*dq                    !select q   
       call NucleonFF(q,GEp,GEn)       !compute nucleon FFs
       FFc(iq)=Gep*FFp(iq)+Gen*FFn(iq) !charge form factor
      end do                           !close iq-loop
      call Density(FFc,Rhoc)           !compute charge density
c
      return
      end
c----------------------------------------------------------------------
c     Having read proton and neutron densities and computed 
c     their corresponding form factors, compute weak-charge
c     density and its corresponding form factor.
c
      subroutine WkRhoFF
c
      implicit real*8 (a-h,j-z)  
      parameter(iPts=1000)
      common/QGrid/qmin,qmax,dq,iQmax
      common/pDensity/Rhop(0:iPts),FFp(0:iPts)
      common/nDensity/Rhon(0:iPts),FFn(0:iPts)
      common/wDensity/Rhow(0:iPts),FFw(0:iPts)
c
c     Compute nuclear charge form factor:
      qwp=+0.0721d0                       !(1-4*sin^2)+radiative
      qwn=-0.9878d0                       !(-1)+radiative corrections
      Gs=0.000000d0                       !ignore strange quarks
      do iq=0,iQmax                       !loop over QGrid
       q=qmin+iq*dq                       !select q   
       call NucleonFF(q,GEp,GEn)          !compute nucleon FFs
       FFw1=qwp*(Gep*FFp(iq)+Gen*FFn(iq)) !"1st" contribution
       FFw2=qwn*(Gen*FFp(iq)+Gep*FFn(iq)) !"2nd" contribution
       FFw(iq)=FFw1+FFw2                  !weak form factor
      end do                              !close iq-loop
      call Density(FFw,Rhow)              !compute weak-charge density
c
      return
      end
c----------------------------------------------------------------------
c     Having computed proton, neutron, charge, and weak-charge
c     densities and form factors, compute the equivalent Helm
c     and symmetric-Fermi densities and form factors.
c
      subroutine HelmSFermi   
c
      implicit real*8 (a-h,j-z)  
      parameter(iPts=1000)
      common/QGrid/qmin,qmax,dq,iQmax
      common/RGrid/rmin,rmax,dr,iRmax
      common/pDensity/Rhop(0:iPts),FFp(0:iPts)
      common/nDensity/Rhon(0:iPts),FFn(0:iPts)
      common/cDensity/Rhoc(0:iPts),FFc(0:iPts)
      common/wDensity/Rhow(0:iPts),FFw(0:iPts)
      common/pSFDensity/RhopSF(0:iPts),FFpSF(0:iPts),cp,bp
      common/nSFDensity/RhonSF(0:iPts),FFnSF(0:iPts),cn,bn
      common/cSFDensity/RhocSF(0:iPts),FFcSF(0:iPts),cc,bc
      common/wSFDensity/RhowSF(0:iPts),FFwSF(0:iPts),cw,bw
      common/pHmDensity/RhopHm(0:iPts),FFpHm(0:iPts),R0p,sp
      common/nHmDensity/RhonHm(0:iPts),FFnHm(0:iPts),R0n,sn
      common/cHmDensity/RhocHm(0:iPts),FFcHm(0:iPts),R0c,sc
      common/wHmDensity/RhowHm(0:iPts),FFwHm(0:iPts),R0w,sw
c
c     Start by normalizing form factors to 1 at q=0:
      FFp0=FFp(0)                      !p-normalization
      FFn0=FFn(0)                      !n-normalization
      FFc0=FFc(0)                      !c-normalization
      FFw0=FFw(0)                      !w-normalization
      do iq=0,iQmax                    !loop over QGrid
       FFp(iq)=FFp(iq)/FFp0            !p-FF normalized to 1
       FFn(iq)=FFn(iq)/FFn0            !n-FF normalized to 1
       FFc(iq)=FFc(iq)/FFc0            !c-FF normalized to 1
       FFw(iq)=FFw(iq)/FFw0            !w-FF normalized to 1
      end do                           !close iq-loop      
c
c     Compute equivalent SFermi parameters:
      call FermiParams(Rhop,cp,bp)     !p density
      call FermiParams(Rhon,cn,bn)     !n density
      call FermiParams(Rhoc,cc,bc)     !c density
      call FermiParams(Rhow,cw,bw)     !w density
c
c     Compute equivalent SFermi densities and FFs:
      do in=0,iQmax                    !loop over QR-Grid
       r=rmin+in*dr                    !select r
       q=qmin+in*dq                    !select q   
       RhopSF(in)=FFp0*RhoSF(r,cp,bp)  !pSF density normalized to Z
       RhonSF(in)=FFn0*RhoSF(r,cn,bn)  !nSF density normalized to N
       RhocSF(in)=FFc0*RhoSF(r,cc,bc)  !cSF density normalized to Z
       RhowSF(in)=FFw0*RhoSF(r,cw,bw)  !wSF density normalized to Qw
       FFpSF(in)=FFSF(q,cp,bp)         !pSF Form Factor
       FFnSF(in)=FFSF(q,cn,bn)         !nSF Form Factor
       FFcSF(in)=FFSF(q,cc,bc)         !cSF Form Factor
       FFwSF(in)=FFSF(q,cw,bw)         !wSF Form Factor
      end do                           !close in-loop 
c
c     Compute equivalent Helm parameters:
      call HelmParams(FFp,R0p,sp)      !p Form Factor
      call HelmParams(FFn,R0n,sn)      !n Form Factor
      call HelmParams(FFc,R0c,sc)      !c Form Factor
      call HelmParams(FFw,R0w,sw)      !w Form Factor
c
c     Compute equivalent Helm densities and FFs:
      do in=0,iQmax                    !loop over QR-Grid
       r=rmin+in*dr                    !select r
       q=qmin+in*dq                    !select q   
       RhopHm(in)=FFp0*RhoHm(r,R0p,sp) !pHm density normalized to Z
       RhonHm(in)=FFn0*RhoHm(r,R0n,sn) !nHm density normalized to N
       RhocHm(in)=FFc0*RhoHm(r,R0c,sc) !cHm density normalized to Z
       RhowHm(in)=FFw0*RhoHm(r,R0w,sw) !wHm density normalized to Qw
       FFpHm(in)=FFHm(q,R0p,sp)        !pHm Form Factor
       FFnHm(in)=FFHm(q,R0n,sn)        !nHm Form Factor
       FFcHm(in)=FFHm(q,R0c,sc)        !cHm Form Factor
       FFwHm(in)=FFHm(q,R0w,sw)        !wHm Form Factor
      end do                           !close in-loop
c
      return
      end
c----------------------------------------------------------------------
c     Write proton, neutron, charge, and weak-charge densities and
c     form factors for the relativistic mean-field model and its
c     symmetric Fermi and Helm equivalents.
c
      subroutine Writer 
c
      implicit real*8 (a-h,j-z)  
c
      call pWriter             !write p-density and FF
      call nWriter             !write n-density and FF
      call cWriter             !write c-density and FF
      call wWriter             !write w-density and FF
c     
      return
      end
c----------------------------------------------------------------------
c     Proton RMF, Symmetric Fermi, and Helm density and Form Factor.
c
      subroutine pWriter 
c
      implicit real*8 (a-h,j-z)  
      parameter(iPts=1000)
      common/QGrid/qmin,qmax,dq,iQmax
      common/RGrid/rmin,rmax,dr,iRmax
      common/pDensity/Rhop(0:iPts),FFp(0:iPts)
      common/pSFDensity/RhopSF(0:iPts),FFpSF(0:iPts),cp,bp
      common/pHmDensity/RhopHm(0:iPts),FFpHm(0:iPts),R0p,sp
c
      open(unit=10,file="Ar40p.out")
c     RMF parameters (norm, <r^2>, and <r^4>):
      pi=4*atan(1.0d0)              !famous constant
      call RhoMoms(Rhop,n0,r2,r4)   !norm, <r^2>, and <r^4>
c
c     Write relevant bulk parameters to file:
      r2SF=(3*cp**2/5)+(7*(pi*bp)**2/5) 
      r4SF=(3*cp**4/7)+(18*(pi*bp*cp)**2/7)+(31*(pi*bp)**4/7) 
      r2Hm=(3*R0p**2/5)+(3*sp**2)
      r4Hm=(3*R0p**4/7)+(6*(R0p*sp)**2)+(15*sp**4)
      write(10,100)
      write(10,200) n0,sqrt(r2),r4**(0.250d0)
      write(10,300) cp,bp,sqrt(r2SF),sqrt(sqrt(r4SF))
      write(10,400) R0p,sp,sqrt(r2Hm),sqrt(sqrt(r4Hm))
c     
c     Write densities and FFs to file:
      write(10,'("#")')             !write to file
      write(10,500)                 !write to file
      do in=0,iRmax                 !loop over QR-grid
       r0=rmin+in*dr                !select value of r
       q0=qmin+in*dq                !select value of q
       Rho0=Rhop(in)                !RMF density
       Rho1=RhopSF(in)              !SFermi density
       Rho2=RhopHm(in)              !Helm density
       FF0=FFp(in)                  !RMF Form Factor
       FF1=FFpSF(in)                !SFermi Form Factor 
       FF2=FFpHm(in)                !Helm Form Factor
       write(10,600) r0,Rho0,Rho1,Rho2,q0,FF0,FF1,FF2 
      end do                        !close in-loop   
      close(unit=10)                !close file 
c
 100  format('#',2x,'   Proton Density and Form Factor:')
 200  format('#',2x,'   RMF:',21x,'Z=',f8.3,
     .2x,'<r^2>^(1/2)[fm]=',f5.3,2x,'<r^4>^(1/4)[fm]=',f5.3)
 300  format('#',2x,'SFermi:',2x,' c[fm]=',f5.3,2x,' beta[fm]=',f5.3,
     .2x,'<r^2>^(1/2)[fm]=',f5.3,2x,'<r^4>^(1/4)[fm]=',f5.3)
 400  format('#',2x,'  Helm:',2x,'R0[fm]=',f5.3,2x,'sigma[fm]=',f5.3,
     .2x,'<r^2>^(1/2)[fm]=',f5.3,2x,'<r^4>^(1/4)[fm]=',f5.3)
 500  format('#',8x,'r',9x,'RhoRMF',8x,'RhoSFermi',6x,'RhoHelm'
     .     10x,'q',9x,'FFRMF',9x,'FFSFermi',8x,'FFHelm')
 600  format(2x,2(' ',f10.6,2x,d13.6,2x,d13.6,2x,d13.6))
c
      return
      end
c-----------------------------------------------------------------------
c     Neutron RMF, Symmetric Fermi, and Helm density and Form Factor.
c
      subroutine nWriter 
c
      implicit real*8 (a-h,j-z)  
      parameter(iPts=1000)
      common/QGrid/qmin,qmax,dq,iQmax
      common/RGrid/rmin,rmax,dr,iRmax
      common/nDensity/Rhon(0:iPts),FFn(0:iPts)
      common/nSFDensity/RhonSF(0:iPts),FFnSF(0:iPts),cn,bn
      common/nHmDensity/RhonHm(0:iPts),FFnHm(0:iPts),R0n,sn
c
      open(unit=11,file="Ar40n.out")
c     RMF parameters (norm, <r^2>, and <r^4>):
      pi=4*atan(1.0d0)              !famous constant
      call RhoMoms(Rhon,n0,r2,r4)   !norm, <r^2>, and <r^4>
c
c     Write relevant bulk parameters to file:
      r2SF=(3*cn**2/5)+(7*(pi*bn)**2/5) 
      r4SF=(3*cn**4/7)+(18*(pi*bn*cn)**2/7)+(31*(pi*bn)**4/7) 
      r2Hm=(3*R0n**2/5)+(3*sn**2)
      r4Hm=(3*R0n**4/7)+(6*(R0n*sn)**2)+(15*sn**4)
      write(11,100)
      write(11,200) n0,sqrt(r2),r4**(0.250d0)
      write(11,300) cn,bn,sqrt(r2SF),sqrt(sqrt(r4SF))
      write(11,400) R0n,sn,sqrt(r2Hm),sqrt(sqrt(r4Hm))
c     
c     Write densities and FFs to file:
      write(11,'("#")')             !write to file
      write(11,500)                 !write to file
      do in=0,iRmax                 !loop over QR-grid
       r0=rmin+in*dr                !select value of r
       q0=qmin+in*dq                !select value of q
       Rho0=Rhon(in)                !RMF density
       Rho1=RhonSF(in)              !SFermi density
       Rho2=RhonHm(in)              !Helm density
       FF0=FFn(in)                  !RMF Form Factor
       FF1=FFnSF(in)                !SFermi Form Factor 
       FF2=FFnHm(in)                !Helm Form Factor
       write(11,600) r0,Rho0,Rho1,Rho2,q0,FF0,FF1,FF2 
      end do                        !close in-loop 
      close(unit=11)                !close file
c
 100  format('#',2x,'   Neutron Density and Form Factor:')
 200  format('#',2x,'   RMF:',21x,'N=',f8.3,
     .2x,'<r^2>^(1/2)[fm]=',f5.3,2x,'<r^4>^(1/4)[fm]=',f5.3)
 300  format('#',2x,'SFermi:',2x,' c[fm]=',f5.3,2x,' beta[fm]=',f5.3,
     .2x,'<r^2>^(1/2)[fm]=',f5.3,2x,'<r^4>^(1/4)[fm]=',f5.3)
 400  format('#',2x,'  Helm:',2x,'R0[fm]=',f5.3,2x,'sigma[fm]=',f5.3,
     .2x,'<r^2>^(1/2)[fm]=',f5.3,2x,'<r^4>^(1/4)[fm]=',f5.3)
 500  format('#',8x,'r',9x,'RhoRMF',8x,'RhoSFermi',6x,'RhoHelm'
     .     10x,'q',9x,'FFRMF',9x,'FFSFermi',8x,'FFHelm')
 600  format(2x,2(' ',f10.6,2x,d13.6,2x,d13.6,2x,d13.6))
c
      return
      end
c-----------------------------------------------------------------------
c     Charge RMF, Symmetric Fermi, and Helm density and Form Factor.
c
      subroutine cWriter 
c
      implicit real*8 (a-h,j-z)  
      parameter(iPts=1000)
      common/QGrid/qmin,qmax,dq,iQmax
      common/RGrid/rmin,rmax,dr,iRmax
      common/cDensity/Rhoc(0:iPts),FFc(0:iPts)
      common/cSFDensity/RhocSF(0:iPts),FFcSF(0:iPts),cc,bc
      common/cHmDensity/RhocHm(0:iPts),FFcHm(0:iPts),R0c,sc
c
      open(unit=12,file="Ar40c.out")
c     RMF parameters (norm, <r^2>, and <r^4>):
      pi=4*atan(1.0d0)              !famous constant
      call RhoMoms(Rhoc,n0,r2,r4)   !norm, <r^2>, and <r^4>
c
c     Write relevant bulk parameters to file:
      r2SF=(3*cc**2/5)+(7*(pi*bc)**2/5) 
      r4SF=(3*cc**4/7)+(18*(pi*bc*cc)**2/7)+(31*(pi*bc)**4/7) 
      r2Hm=(3*R0c**2/5)+(3*sc**2)
      r4Hm=(3*R0c**4/7)+(6*(R0c*sc)**2)+(15*sc**4)
      write(12,100)
      write(12,200) n0,sqrt(r2),r4**(0.250d0)
      write(12,300) cc,bc,sqrt(r2SF),sqrt(sqrt(r4SF))
      write(12,400) R0c,sc,sqrt(r2Hm),sqrt(sqrt(r4Hm))
c     
c     Write densities and FFs to file:
      write(12,'("#")')             !write to file
      write(12,500)                 !write to file
      do in=0,iRmax                 !loop over QR-grid
       r0=rmin+in*dr                !select value of r
       q0=qmin+in*dq                !select value of q
       Rho0=Rhoc(in)                !RMF density
       Rho1=RhocSF(in)              !SFermi density
       Rho2=RhocHm(in)              !Helm density
       FF0=FFc(in)                  !RMF Form Factor
       FF1=FFcSF(in)                !SFermi Form Factor 
       FF2=FFcHm(in)                !Helm Form Factor
       write(12,600) r0,Rho0,Rho1,Rho2,q0,FF0,FF1,FF2 
      end do                        !close in-loop     
      close(unit=12)                !close file
c
 100  format('#',2x,'   Charge Density and Form Factor:')
 200  format('#',2x,'   RMF:',21x,'Z=',f8.3,
     .2x,'<r^2>^(1/2)[fm]=',f5.3,2x,'<r^4>^(1/4)[fm]=',f5.3)
 300  format('#',2x,'SFermi:',2x,' c[fm]=',f5.3,2x,' beta[fm]=',f5.3,
     .2x,'<r^2>^(1/2)[fm]=',f5.3,2x,'<r^4>^(1/4)[fm]=',f5.3)
 400  format('#',2x,'  Helm:',2x,'R0[fm]=',f5.3,2x,'sigma[fm]=',f5.3,
     .2x,'<r^2>^(1/2)[fm]=',f5.3,2x,'<r^4>^(1/4)[fm]=',f5.3)
 500  format('#',8x,'r',9x,'RhoRMF',8x,'RhoSFermi',6x,'RhoHelm'
     .     10x,'q',9x,'FFRMF',9x,'FFSFermi',8x,'FFHelm')
 600  format(2x,2(' ',f10.6,2x,d13.6,2x,d13.6,2x,d13.6))
c
      return
      end
c-----------------------------------------------------------------------
c     Weak RMF, Symmetric Fermi, and Helm density and Form Factor.
c
      subroutine wWriter 
c
      implicit real*8 (a-h,j-z)  
      parameter(iPts=1000)
      common/QGrid/qmin,qmax,dq,iQmax
      common/RGrid/rmin,rmax,dr,iRmax
      common/wDensity/Rhow(0:iPts),FFw(0:iPts)
      common/wSFDensity/RhowSF(0:iPts),FFwSF(0:iPts),cw,bw
      common/wHmDensity/RhowHm(0:iPts),FFwHm(0:iPts),R0w,sw
      dimension qs(iPts),bF(iPts),cF(iPts),dF(iPts)
      dimension FFMF(iPts),FFSF(iPts),FFHm(iPts)
c
      open(unit=13,file="Ar40w.out")
c     RMF parameters (norm, <r^2>, and <r^4>):
      pi=4*atan(1.0d0)              !famous constant
      call RhoMoms(Rhow,n0,r2,r4)   !norm, <r^2>, and <r^4>
c
c     Write relevant bulk parameters to file:
      r2SF=(3*cw**2/5)+(7*(pi*bw)**2/5) 
      r4SF=(3*cw**4/7)+(18*(pi*bw*cw)**2/7)+(31*(pi*bw)**4/7) 
      r2Hm=(3*R0w**2/5)+(3*sw**2)
      r4Hm=(3*R0w**4/7)+(6*(R0w*sw)**2)+(15*sw**4)
      write(13,100)
      write(13,200) n0,sqrt(r2),r4**(0.250d0)
      write(13,300) cw,bw,sqrt(r2SF),sqrt(sqrt(r4SF))
      write(13,400) R0w,sw,sqrt(r2Hm),sqrt(sqrt(r4Hm))
c     
c     Write densities and FFs to file:
      write(13,'("#")')             !write to file
      write(13,500)                 !write to file
      do in=0,iRmax                 !loop over QR-grid
       r0=rmin+in*dr                !select value of r
       q0=qmin+in*dq                !select value of q
       Rho0=Rhow(in)                !RMF density
       Rho1=RhowSF(in)              !SFermi density
       Rho2=RhowHm(in)              !Helm density
       FF0=FFw(in)                  !RMF Form Factor
       FF1=FFwSF(in)                !SFermi Form Factor 
       FF2=FFwHm(in)                !Helm Form Factor
       write(13,600) r0,Rho0,Rho1,Rho2,q0,FF0,FF1,FF2 
      end do                        !close in-loop     
      close(unit=13)                !close file
c
 100  format('#',2x,'   Weak-Charge Density and Form Factor:')
 200  format('#',2x,'   RMF:',20x,'Qw=',f8.3,
     .2x,'<r^2>^(1/2)[fm]=',f5.3,2x,'<r^4>^(1/4)[fm]=',f5.3)
 300  format('#',2x,'SFermi:',2x,' c[fm]=',f5.3,2x,' beta[fm]=',f5.3,
     .2x,'<r^2>^(1/2)[fm]=',f5.3,2x,'<r^4>^(1/4)[fm]=',f5.3)
 400  format('#',2x,'  Helm:',2x,'R0[fm]=',f5.3,2x,'sigma[fm]=',f5.3,
     .2x,'<r^2>^(1/2)[fm]=',f5.3,2x,'<r^4>^(1/4)[fm]=',f5.3)
 500  format('#',8x,'r',9x,'RhoRMF',8x,'RhoSFermi',6x,'RhoHelm'
     .     10x,'q',9x,'FFRMF',9x,'FFSFermi',8x,'FFHelm')
 600  format(2x,2(' ',f10.6,2x,d13.6,2x,d13.6,2x,d13.6))
c
      return
      end
c-----------------------------------------------------------------------
c     Compute the form factor FF(q) on a previously defined "QGrid" 
c     as the Fourier transform of a density Rho(r) defined on a 
c     previously defined "RGrid". A simple (4-2-4...) Simpson rule 
c     is used to perform the required integral.
c
      Subroutine FormFactor(Rho,FF)

      implicit real*8 (a-h,j-z) 
      parameter(iPts=1000)
      common/QGrid/qmin,qmax,dq,iQmax
      common/RGrid/rmin,rmax,dr,iRmax
      dimension Rho(0:iPts),FF(0:iPts)
c
c     Define useful constants:
      eps=1.00d-20               !tiny number
      pi=4*atan(1.0d0)           !famous constant
c
c     Fourier transform:
      do iq=0,iQmax              !loop over QGrid
       iwgt=2                    !initialize weight
       FF(iq)=0.0d0              !initialize FF(q)
       q=qmin+iq*dq              !select q
       do ir=1,iRmax             !loop over RGrid
        iwgt=6-iwgt              !adjust Simpson weight              
        wgt=iwgt*dr/3            !define Simpson weight
        r=rmin+ir*dr             !select r
        x=q*r+eps                !argument of Bessel
        j0=sin(x)/x              !spherical Bessel0
        FF(iq)=FF(iq)+wgt*4*pi*r**2*Rho(ir)*j0 !update FF(q)
       end do                    !close ir-loop
      end do                     !close iq-loop
c
      return
      end
c-----------------------------------------------------------------------
c     Compute the density Rho(r) on a previously defined "RGrid"
c     as the Fourier transform of a form factor FF(q) defined on 
c     a previously defined "QGrid". A simple (4-2-4...) Simpson
c     rule is used to perform the required integral.
c
      Subroutine Density(FF,Rho)

      implicit real*8 (a-h,j-z) 
      parameter(iPts=1000)
      common/QGrid/qmin,qmax,dq,iQmax
      common/RGrid/rmin,rmax,dr,iRmax
      dimension Rho(0:iPts),FF(0:iPts)
c
c     Define useful constants:
      eps=1.00d-20               !tiny number
      pi=4*atan(1.0d0)           !famous constant
c
c     Fourier transform:
      do ir=0,iRmax              !loop over RGrid
       iwgt=2                    !initialize weight
       Rho(ir)=0.0d0             !initialize Rho(r)
       r=rmin+ir*dr              !select r
       do iq=1,iQmax             !loop over QGrid
        iwgt=6-iwgt              !adjust Simpson weight              
        wgt=iwgt*dq/3            !define Simpson weight
        q=qmin+iq*dq             !select q
        x=q*r+eps                !argument of Bessel
        j0=sin(x)/x              !spherical Bessel0
        Rho(ir)=Rho(ir)+wgt*4*pi*q**2*FF(iq)*j0 !update Rho(r)
       end do                    !close iq-loop
       Rho(ir)=Rho(ir)/(2*pi)**3 !final normalization
      end do                     !close ir-loop
c
      return
      end
c-----------------------------------------------------------------------
c     Compute the first three moments of the density Rho(r): 
c     normalization, <r^2>, and <r^4>. A simple (4-2-4...) 
c     Simpson rule is used to perform the required integrals.      

      subroutine RhoMoms(Rho,n0,r2,r4)
c
      implicit real*8 (a-h,j-z)  
      parameter(iPts=1000)
      dimension Rho(0:iPts)
      common/RGrid/rmin,rmax,dr,iRmax
c
c     Compute norm, <r^2>, and <r^4>:
      pi=4*atan(1.0d0)              !famous constant
      n0=0.000000000d0              !initialize norm
      r2=0.000000000d0              !initialize <r^2>        
      r4=0.000000000d0              !initialize <r^4>
      iwgt=2                        !initialize iwgt
      do ir=1,iRmax                 !loop over r-grid
       iwgt=6-iwgt                  !update Simpson weight
       wgt=iwgt*dr/3                !Simpson weight
       r0=rmin+ir*dr                !select r
       Rho0=Rho(ir)                 !extract Rho
       n0=n0+wgt*4*pi*r0**2*Rho0    !update norm
       r2=r2+wgt*4*pi*r0**4*Rho0    !update <r^2>
       r4=r4+wgt*4*pi*r0**6*Rho0    !update <r^4>
      end do                        !close ir-loop
      r2=r2/n0                      !<r^2> properly normalized
      r4=r4/n0                      !<r^4> properly normalized
c
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     Routine to compute single-nucleon charge form factors.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     Compute single-nucleon charge form factors.
c
      Subroutine NucleonFF(q0,GEp0,GEn0)
c
      implicit real*8 (a-h,j-z) 
      
c     Define nucleon mass (in Gev) and nucleon magnetic moments:
      hbarc=197.326d0           !hbarc
      q=q0*hbarc/1000           !convert q0 to GeV
      GEp0=GEp(q)               !electric FF proton
      GEn0=GEn(q)               !electric FF neutron
C      GEn0=0.00d0               !set neutron electric FF to zero
c
      return
      end
c-----------------------------------------------------------------------
c     For a given value of q (in GeV), compute the proton electric form 
c     factor GEp [Friedrich&Walcher Eur. Phys. J. A17, 607 (2003)].
c
      Function GEp(q)
c
      implicit real*8 (a-h,j-z) 
c
c     Qsquare in GeV^2:
      Q2=q**2                         
c
c     Parameters for the smooth component:
      a10=+1.041d0
      a11=+0.765d0
      a20=-0.041d0
      a21=+6.200d0
c     Parameters for the "bumpy" components:
      ab=-0.0090d0      
      Qb=+0.0700d0     
      sb=+0.2700d0     
c
c     Smooth (dipole) component:
      Dipole1=a10/(1+Q2/a11)**2        
      Dipole2=a20/(1+Q2/a21)**2        
      Gsmooth=Dipole1+Dipole2         
c
c     Bumpy (gaussian) component:
      Qm=Q-Qb                         
      Qp=Q+Qb                         
      Gaussm=exp(-0.50d0*(Qm/sb)**2)  
      Gaussp=exp(-0.50d0*(Qp/sb)**2)  
      Gbump=ab*Q2*(Gaussm+Gaussp)        
c
c     Proton electric form factor:
      GEp=Gsmooth+Gbump                
c
      return
      end
c-----------------------------------------------------------------------
c     For a given value of q (in GeV), compute the neutron electric form 
c     factor GEn [Friedrich&Walcher Eur. Phys. J. A17, 607 (2003)].
c
      Function GEn(q)
c
      implicit real*8 (a-h,j-z) 
c
c     Qsquare in GeV^2:
      Q2=q**2                         
c
c     Parameters for the smooth component:
      a10=+1.040d0
      a11=+1.730d0
      a20=-1.040d0
      a21=+1.540d0
c     Parameters for the "bumpy" components:
      ab=+0.0090d0      
      Qb=+0.2900d0     
      sb=+0.2000d0     
c
c     Smooth (dipole) component:
      Dipole1=a10/(1+Q2/a11)**2        
      Dipole2=a20/(1+Q2/a21)**2        
      Gsmooth=Dipole1+Dipole2         
c
c     Bumpy (gaussian) component:
      Qm=Q-Qb                         
      Qp=Q+Qb                         
      Gaussm=exp(-0.50d0*(Qm/sb)**2)  
      Gaussp=exp(-0.50d0*(Qp/sb)**2)  
      Gbump=ab*Q2*(Gaussm+Gaussp)        
c
c     Neutron electric form factor:
      GEn=Gsmooth+Gbump                
c
      return
      end
c-----------------------------------------------------------------------
c     For a given value of q (in GeV), compute the neutron electric form 
c     factor GEn [Geis et al., PRL101, 042501 (2008) - with NO BUMP].
c
      Function GEnPRL(q)
c
      implicit real*8 (a-h,j-z) 
c
c     Qsquare in GeV^2:
      Q2=q**2                         
c
c     Parameters for the smooth component:
      a10=+0.095d0
      a11=+2.770d0
      a20=-0.095d0
      a21=+0.339d0
c
c     Smooth (dipole) component:
      Dipole1=a10/(1+Q2/a11)**2        
      Dipole2=a20/(1+Q2/a21)**2        
      Gsmooth=Dipole1+Dipole2         
c
c     Neutron electric form factor:
      GEnPRL=Gsmooth                
c
      return
      end
c-----------------------------------------------------------------------      
c-----------------------------------------------------------------------
c     Density and Form Factor of a symmetrized Fermi Function
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     Given a density Rho(r) with mean-square radius <r^2> and 
c     mean-fourth radius <r^4>, extract the Fermi parameters c 
c     and beta of the corresponding symmetrized Fermi function.
c
      subroutine FermiParams(Rho,c,beta)
c
      implicit real*8 (a-h,j-z)  
      parameter(iPts=1000)
      dimension Rho(0:iPts)
      common/RGrid/rmin,rmax,dr,iRmax
c
c     Compute <r^2> and <r^4>:
      pi=4*atan(1.0d0)             !famous constant
      call RhoMoms(Rho,n0,r2,r4)   !norm, <r^2>, and <r^4>
c
c     Fermi parameters c and beta from r2 and r4 (see notes):
      A1=5*r2                      !definition of A 
      B2=7*r4                      !definition of B*B 
      Sq=sqrt(12*B2-3*A1**2)       !square-root factor
      c=sqrt((15*A1-7*Sq)/24)      !Fermi parameter "c"
      beta=sqrt((Sq-A1)/8)/pi      !Fermi parameter "beta"
c
      return
      end
c-----------------------------------------------------------------------
c     Compute the density of a symmetrized Fermi function at point r
c     (density has been normalized to 1). 
c
      Function RhoSF(r,c,beta) 
c
      implicit real*8 (a-h,j-z)
c
      pi=4*atan(1.0d0)                                    !famous const.  
      Rho0=(3.0d0/(4*pi*c**3))/(1+(pi*beta/c)**2)         !normalization
      RhoSF=Rho0*sinh(c/beta)/(cosh(r/beta)+cosh(c/beta)) !symmetrized 
c
      return
      end
c----------------------------------------------------------------------
c     Compute the form factor of a symmetrized Fermi function
c     at point q.
c
      Function FFSF(q,c,beta) 
c
      implicit real*8 (a-h,j-z)
c
c     Analytic Form Factor for symmetrized Fermi function:
      eps=1.00d-05                      !tiny number
      pi=4*atan(1.0d0)                  !famous constant
      k=q*beta                          !useful definition 
      xi=c/beta                         !useful definition
      if(k.lt.eps) k=eps                !avoid dividing by zero
      FF0=(3*pi)/(xi*(xi**2+pi**2))     !normalization factor     
      FF1=(pi*sin(xi*k)/tanh(pi*k)-xi*cos(xi*k))/(k*sinh(pi*k))
      FFSF=FF0*FF1                      !analytic result
c
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     Density and Form Factor of the Helm function
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     Given a Form Factor FF(q), extract the Helm parameters R0 
c     and sigma of the corresponding symmetrized Helm function.
c
      subroutine HelmParams(FF,R0,sigma)
c
      implicit real*8 (a-h,j-z)  
      parameter(iQs=1000)
      external FFSp
      dimension FF(0:iQs)
      common/QGrid/qmin,qmax,dq,iQmax
      common/FFSpline/qs(iQs),FFs(iQs),bF(iQs),cF(iQs),dF(iQs)
c
c     Spline the form factor FF(q):
      do iq=0,iQmax                          !loop over QGrid
       q0=qmin+iq*dq                         !select q
       qs(iq+1)=q0                           !store q in array
       FFs(iq+1)=FF(iq)                      !store FF in array 
      end do                                 !close iq-loop
      call csplin(iQmax,qs,FFs,bF,cF,dF)     !spline FF(q)
c
c     Find the first zero of the form factor:
      q1=0.70000d0                           !q_minimun
      q2=1.4000d0                            !q-maximum
      tol=1.00d-10                           !tolerance 
      q0=ZBRENT(FFSp,q1,q2,tol)              !find root
c
c     Compute diffraction radius R0:
      x0=4.493409d0                          !zero of FHelm(q)
      R0=x0/q0                               !diffraction radius
c
c     Compute difuseness parameter sigma:
      xm=5.600d0                             !qm near 1st maximum
      qm=xm/R0                               !select qm
      FFm=cseval(iQmax,qm,0,qs,FFs,bF,cF,dF) !FF(qm)
      j1m=(sin(xm)/xm**2)-(cos(xm)/xm)       !spherical Bessel
      Fboxm=3*j1m/xm                         !Fbox(qm)
      sigma2=(2.0d0/qm**2)*log(Fboxm/FFm)    !sigma^2
      sigma=sqrt(sigma2)                     !sigma
c
      return
      end
c----------------------------------------------------------------------
c     Use cubic spline to define form factor at an arbitrary point q.
c
      Function FFSp(q)
c
      implicit real*8 (a-h,j-z) 
      parameter(iQs=1000)
      common/QGrid/qmin,qmax,dq,iQmax
      common/FFSpline/qs(iQs),FFs(iQs),bF(iQs),cF(iQs),dF(iQs)
c
      FFSp=cseval(iQmax,q,0,qs,FFs,bF,cF,dF)  !FF(q)      
c
      return
      end
c----------------------------------------------------------------------
c     Compute the density of a Helm function at point r
c     (density has been normalized to 1). 
c
      Function RhoHm(r,R0,sigma) 
c
      implicit real*8 (a-h,j-z)
c
c     Analytic Helm Density:
      eps=1.00d-05                            !tiny number
      pi=4*atan(1.0d0)                        !famous const.
      x=r/sigma                               !useful quantity
      x0=R0/sigma                             !useful quantity
      if(x.lt.eps) x=eps                      !avoid dividing by zero   
      Rho0=(3.0d0/(4*pi*R0**3))               !normalization
      norm1=Rho0/2.000000d0                   !normalization 1
      norm2=Rho0/sqrt(2*pi)                   !normalization 2
      RhoHm1=+norm1*erf((x+x0)/sqrt(2.0d0))   !"1st" term
      RhoHm2=-norm1*erf((x-x0)/sqrt(2.0d0))   !"2nd" term 
      RhoHm3=+norm2*exp(-0.50d0*(x+x0)**2)/x  !"3rd" term
      RhoHm4=-norm2*exp(-0.50d0*(x-x0)**2)/x  !"4th" term
      RhoHm=RhoHm1+RhoHm2+RhoHm3+RhoHm4       !Helm density
c
      return
      end
c-----------------------------------------------------------------------
c     Compute the Helm form factor at point q.
c
      Function FFHm(q,R0,sigma)  
c
      implicit real*8 (a-h,j-z)
c
c     Analytic Helm Form Factor:
      eps=1.00d-05                      !tiny number
      pi=4*atan(1.0d0)                  !famous constant
      x=q*R0                            !useful quantity
      y=q*sigma                         !useful quantity
      if(x.lt.eps) x=eps                !avoid dividing by zero
      j1=(sin(x)/x**2)-(cos(x)/x)       !spherical Bessel
      FFHm=(3*j1/x)*exp(-y**2/2)        !Helm form factor
c
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     THIS SECTION INCLUDES SPLINES ROUTINES
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
C     LAST MODIFIED 4-10-76 BY W. M. COUGHRAN, JR.
      SUBROUTINE CSPLIN (N, X, Y, B, C, D)
      INTEGER N
      DOUBLE PRECISION X(N), Y(N), B(N), C(N), D(N)
C
C  THIS IS A MODIFIED VERSION OF SPLINE, DESCRIBED IN THE NOTES
C  BY FORSYTHE MALCOLM AND MOLER
C
C  THE COEFFICIENTS B(I), C(I), AND D(I), I=1,2,...,N ARE COMPUTED
C  FOR A CUBIC INTERPOLATING SPLINE
C
C    S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
C
C    FOR  X(I) .LE. X .LE. X(I+1)
C
C  INPUT..
C
C    N = THE NUMBER OF DATA POINTS OR KNOTS (N.GE.2)
C    X = THE ABSCISSAS OF THE KNOTS IN STRICTLY INCREASING ORDER
C    Y = THE ORDINATES OF THE KNOTS
C
C  OUTPUT..
C
C    B, C, D  = ARRAYS OF SPLINE COEFFICIENTS AS DEFINED ABOVE.
C
C  USING  P  TO DENOTE DIFFERENTIATION,
C
C    Y(I) = S(X(I))
C    B(I) = SP(X(I))
C    C(I) = SPP(X(I))/2
C    D(I) = SPPP(X(I))/6  (DERIVATIVE FROM THE RIGHT)
C
C  THE ACCOMPANYING FUNCTION SUBPROGRAM  CSEVAL  CAN BE USED
C  TO EVALUATE THE SPLINE, ITS DERIVATIVE OR EVEN ITS 2ND DERIVATIVE.
C
C       SEE COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS BY
C           FORSYTHE, MALCOLM, AND MOLER FOR DETAILS
C
      INTEGER LOUT, NM1, IB, I
      DOUBLE PRECISION T
      DATA LOUT/6/
C
C         CHECK INPUT FOR CONSISTENCY
C
      IF (N .GE. 2) GO TO 1
      WRITE(LOUT, 1000)N
      RETURN
C
    1 NM1 = N-1
      DO 3 I = 1, NM1
         IF (X(I) .LT. X(I+1)) GO TO 3
         WRITE(LOUT, 1001)
         RETURN
    3 CONTINUE
C
    5 IF ( N .EQ. 2 ) GO TO 50
C
C  SET UP TRIDIAGONAL SYSTEM
C
C  B = DIAGONAL, D = OFFDIAGONAL, C = RIGHT HAND SIDE.
C
      D(1) = X(2) - X(1)
      C(2) = (Y(2) - Y(1))/D(1)
      DO 10 I = 2, NM1
         D(I) = X(I+1) - X(I)
         B(I) = 2.*(D(I-1) + D(I))
         C(I+1) = (Y(I+1) - Y(I))/D(I)
         C(I) = C(I+1) - C(I)
   10 CONTINUE
C
C  END CONDITIONS.  THIRD DERIVATIVES AT  X(1)  AND  X(N)
C  OBTAINED FROM DIVIDED DIFFERENCES
C
      B(1) = -D(1)
      B(N) = -D(N-1)
      C(1) = 0.
      C(N) = 0.
      IF ( N .EQ. 3 ) GO TO 15
      C(1) = C(3)/(X(4)-X(2)) - C(2)/(X(3)-X(1))
      C(N) = C(N-1)/(X(N)-X(N-2)) - C(N-2)/(X(N-1)-X(N-3))
      C(1) = C(1)*D(1)**2/(X(4)-X(1))
      C(N) = -C(N)*D(N-1)**2/(X(N)-X(N-3))
C
C  FORWARD ELIMINATION
C
   15 DO 20 I = 2, N
         T = D(I-1)/B(I-1)
         B(I) = B(I) - T*D(I-1)
         C(I) = C(I) - T*C(I-1)
   20 CONTINUE
C
C  BACK SUBSTITUTION
C
      C(N) = C(N)/B(N)
      DO 30 IB = 1, NM1
         I = N-IB
         C(I) = (C(I) - D(I)*C(I+1))/B(I)
   30 CONTINUE
C
C  C(I) IS NOW THE SIGMA(I) OF THE TEXT
C
C  COMPUTE POLYNOMIAL COEFFICIENTS
C
      B(N) = (Y(N) - Y(NM1))/D(NM1) + D(NM1)*(C(NM1) + 2.*C(N))
      DO 40 I = 1, NM1
         B(I) = (Y(I+1) - Y(I))/D(I) - D(I)*(C(I+1) + 2.*C(I))
         D(I) = (C(I+1) - C(I))/D(I)
         C(I) = 3.*C(I)
   40 CONTINUE
      C(N) = 3.*C(N)
      D(N) = D(N-1)
      RETURN
   50 B(1) = (Y(2)-Y(1))/(X(2)-X(1))
      C(1) = 0.
      D(1) = 0.
      RETURN
 1000 FORMAT('-N < 2 IN CSPLIN CALL--',I10)
 1001 FORMAT('-X IS NOT IN ASCENDING ORDER IN CSPLIN CALL')
      END
      DOUBLE PRECISION FUNCTION CSEVAL (N, U, IDERIV, X, Y, B, C, D)
      INTEGER N, IDERIV
      DOUBLE PRECISION U, X(N), Y(N), B(N), C(N), D(N)
C
C  THIS IS A MODIFIED VERSION OF SEVAL, DESCRIBED IN THE NOTES
C  BY FORSYTHE MALCOLM AND MOLER
C
C  THIS SUBROUTINE EVALUATES THE CUBIC SPLINE FUNCTION, ITS FIRST
C  DERIVATIVE OR SECOND DERIVATIVE, NAMELY:
C
C   IDERIV=0:    Y(I)+B(I)*(U-X(I))+C(I)*(U-X(I))**2+D(I)*(U-X(I))**3
C   IDERIV=1:    B(I)+2*C(I)*(U-X(I))+3*D(I)*(U-X(I))**2
C   IDERIV=2:    2*C(I)+6*D(I)*(U-X(I))
C
C    WHERE  X(I) .LT. U .LT. X(I+1), USING HORNER'S RULE
C
C  IF  U .LT. X(1) THEN  I = 1  IS USED.
C  IF  U .GE. X(N) THEN  I = N  IS USED.
C
C  INPUT..
C
C    N = THE NUMBER OF DATA POINTS
C    U = THE ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED
C    IDERIV = 0 TO EVALUATE S(U)
C           = 1 TO EVALUATE SP(U)
C           = 2 TO EVALUATE SPP(U)
C    X,Y = THE ARRAYS OF DATA ABSCISSAS AND ORDINATES
C    B,C,D = ARRAYS OF SPLINES COEFFICIENTS COMPUTED BY SUBROUTINE CSPLI
C
C  IF  U  IS NOT IN THE SAME INTERVAL AS THE PREVIOUS CALL, THEN A
C    BINARY SEARCH IS PERFORMED TO DETERMINE THE PROPER INTERVAL.
C
      INTEGER I, J, K, LOUT
      DOUBLE PRECISION DX
      DATA I, LOUT/1, 6/
C
      IF (IDERIV .GE. 0 .AND. IDERIV .LE. 2) GO TO 5
      WRITE(LOUT, 1000)IDERIV
      RETURN
C
    5 IF ( I .GE. N ) I = 1
      IF ( U .LT. X(I) ) GO TO 10
      IF ( U .LE. X(I+1) ) GO TO 30
C
C  BINARY SEARCH
C
   10 I = 1
      J = N + 1
   20 K = (I+J)/2
      IF ( U .LT. X(K) ) J = K
      IF ( U .GE. X(K) ) I = K
      IF ( J .GT. I+1 ) GO TO 20
C
C  EVALUATE SPLINE
C
   30 DX = U - X(I)
      J = IDERIV+1
      GO TO (40,50,60), J
C  COMPUTE S(X)
   40 CSEVAL = Y(I) + DX*(B(I) + DX*(C(I) + DX*D(I)))
      RETURN
C  COMPUTE SP(X)
   50 CSEVAL = B(I) + DX*(2*C(I) + 3*DX*D(I))
      RETURN
C  COMPUTE SPP(X)
   60 CSEVAL = 2*C(I) + 6*DX*D(I)
      RETURN
 1000 FORMAT('-IDERIV IS INVALID IN CSEVAL CALL--', I10)
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     THIS SECTION INCLUDES A ROOT-FINDER ROUTINE
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      REAL*8 FUNCTION ZBRENT(FUNC,X1,X2,TOL)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ITMAX=100)
      PARAMETER (EPS=3.D-8)
C
      EXTERNAL FUNC
C
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF(FB*FA.GT.0.)THEN
        WRITE(*,*)'ZBRENT - Root is not bracketed'
        RETURN
        ENDIF
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          ZBRENT=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
11    CONTINUE
      WRITE(*,*)'ZBRENT - Exceeding maximum iterations.'
      ZBRENT=B
      RETURN
      END
c-----------------------------------------------------------------------
