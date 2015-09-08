      parameter (nfd=600,ned=600,nmd=6,nbd=50,nxy=350,ndi=72,

     &            nncn=16,nnbr=30)

c      parameter (nfd=200,ned=200,nmd=6,nbd=50,nxy=200,ndi=72,

c    &            nncn=20,nnbr=40)

c

       common  aaa,bbb,nonc,nhar,mmu,nhh,nsum,itr,ntotal(2,nmd),

     &         npo1,npo2,npo3,npor,nom1,nom2,nom,ncn,nbr,

     &         nelem,nb,nbb,nre,ncc,

     &         df,dmu0,domeg,dpi,

     &         dh(nfd,nbd),dhc(nfd,nncn+nnbr),dhl(nncn+nnbr,nfd),

     &         dhcl(nncn+nnbr,nncn+nnbr),dk(nfd+nncn+nnbr),

     &         dkk(nmd,nfd+nncn+nnbr),da(nfd+nncn+nnbr),

     &         daa(nmd,nfd+nncn+nnbr),daaa(nmd,nfd+nncn+nnbr),

     &         dnds(nnbr,3),ntabl(nnbr,5),dtabl(nnbr,5),

     &         dc(ned),dg(ned),dv(nnbr,nmd,2),

     &         dzz(nnbr,nnbr),dyy(nnbr,nnbr),das(nncn,nnbr),

     &         ds(ned,7),dcs(ned),dbx(ned,nmd,2),dby(ned,nmd,2),

     &         nod(ned,4),xy(nxy,2),dcur(100,nmd,2),dre(100),

     &         nodr(100,2),dd(ned,nmd,2,2),dn(nmd,2,2),

     &         dreluc(0:22,0:1),dsine(0:22,0:ndi),dcosi(0:22,0:ndi)
c **************************************************************        

c ***********  store the temporary solution  dstore  ***********        

c                                                                       

           subroutine dstore                                            

c                                                                       

       include 'common'

c                                                                       

         do 100 i=1,nom                                                 

           daaa(nhh,i)=da(i)                                            

  100    continue                                                       

         return                                                         

         end                                                            

c **************************************************************        

c ****************  check convergence  conv  *******************        

c                                                                       

           subroutine conv                                              

c                                                                       

       include 'common'

c                                                                       

       dimension nlist(1000)                                            

c                                                                       

         nonc=0                                                         

           k=1                                                          

c                                                                       

         damax=0.0                                                      

c                                                                       

         do 100 i=1,nom1                                                

           if(abs(daaa(1,i)).gt.damax)  damax=abs(daaa(1,i))            

100    continue                                                         

c debug                                                                 

          print *, damax                                                

c                                                                       

         damax5=damax*.01                                               

         do 210 i=1,nom1                                                

           derr=(daaa(nhh,i)-daa(nhh,i))/damax                          

           if(abs(daaa(nhh,i)).lt.damax5)         goto 220              

           if(abs(derr).gt.aaa) then                                    

                nonc=nonc+1                                             

                nlist(k)=(i+1)/2                                        

                k=k+1                                                   

           end if                                                       

           goto 210                                                     

  220    continue                                                       

           if (abs(derr).gt.aaa*5.) then                                

               nonc=nonc+1                                              

                nlist(k)=(i+1)/2                                        

                k=k+1                                                   

           end if                                                       

  210    continue                                                       

c                                                                       

*        if (nhh.gt.1)                       go to 400                  

c                                                                       

*        do 300 i=nom1+1,nom                                            

*          if (abs(daaa(nhh,i)).le.1.0e-25) then                        

*          go to 300                                                    

*           end if                                                      

c                                                                       

*            dadv=10.                                                   

*          derr=(daaa(nhh,i)-daa(nhh,i))/daaa(nhh,i)                    

*           if (abs(derr).gt.aaa*dadv) then                             

*               nonc=nonc+1                                             

*               nlist(k)=(i+1)/2                                        

*               k=k+1                                                   

*           end if                                                      

*  300    continue                                                      

*  400    continue                                                      

c                                                                       

c -----  delerating process  -------                                    

c                                                                       

         if((ntotal(1,nhh)+itr).ge.100)                                 

     &                                         bbbb=bbb*0.06            

         if(((ntotal(1,nhh)+itr).ge.90).and.                            

     &    ((ntotal(1,nhh)+itr).lt.100))        bbbb=bbb*.08             

         if(((ntotal(1,nhh)+itr).ge.80).and.                            

     &         ((ntotal(1,nhh)+itr).lt.90))    bbbb=bbb*.1              

         if(((ntotal(1,nhh)+itr).ge.70).and.                            

     &         ((ntotal(1,nhh)+itr).lt.80))    bbbb=bbb*.2              

         if(((ntotal(1,nhh)+itr).ge.60).and.                            

     &         ((ntotal(1,nhh)+itr).lt.70))    bbbb=bbb*.3              

         if(((ntotal(1,nhh)+itr).ge.40).and.                            

     &         ((ntotal(1,nhh)+itr).lt.60))    bbbb=bbb*.4              

         if(((ntotal(1,nhh)+itr).ge.20).and.                            

     &         ((ntotal(1,nhh)+itr).lt.40))    bbbb=bbb*.6              

         if(((ntotal(1,nhh)+itr).ge.10).and.                            

     &         ((ntotal(1,nhh)+itr).lt.20))    bbbb=bbb*.8              

         if((ntotal(1,nhh)+itr).lt.10)         bbbb=bbb                 

c                                                                       

         do 310 i=1,nom1                                                

             daa(nhh,i)=daa(nhh,i)+bbbb*(daaa(nhh,i)-daa(nhh,i))        

  310    continue                                                       

c                                                                       

         do 320 i=nom1,nom                                              

             daa(nhh,i)=daa(nhh,i)+2.*bbbb*(daaa(nhh,i)-daa(nhh,i))     

  320    continue                                                       

                                                                        

c                                                                       

*           print '(20i4)', (nlist(kk),kk=1,k-1)                        

c                                                                       

         return                                                         

         end                                                            

c **************************************************************        

c ****************  check maximum harmonic   *******************        

c                                                                       

           subroutine check                                             

c                                                                       

       include 'common'

c                                                                       

         nonc=0                                                         

         damax1=0.0                                                     

c                                                                       

         do 100 i=1,nom1                                                

           if(abs(daa(1,i)).gt.damax1)    damax1=abs(daa(1,i))          

  100    continue                                                       

c                                                                       

         damax2=0.0                                                     

c                                                                       

         do 200 i=1,nom1                                                

           if(abs(daa(mmu,i)).gt.damax2)  damax2=abs(daa(mmu,i))        

  200    continue                                                       

c                                                                       

         if((damax2/damax1).le.aaa*2.) then                             

             print '('' damax2/damax1='',2e14.5)', damax2/damax1,aaa    

             nonc=0                                                     

         else                                                           

             print '('' damax2/damax1='',2e14.5)', damax2/damax1,aaa    

             nonc=1                                                     

         end if                                                         

c                                                                       

          return                                                        

          end                                                           
c *******************************************************               

c ***************** read data   rdata *******************               

c                                                                       

         subroutine rdata                                               

c                                                                       

       include 'common'

c                                                                       

c ********* if double precision ******************                      

c          implicit real*8(d)                                           

       character*72 char,char1                                          

c                                                                       

c -------- read constant ----------                                     

c                                                                       

c check 2                         ********* open read file *****        

c        open (10,file='ab0169.hbfem03.data(fem400)')                   

c                                                                       

         read (10,1000) char1                                           

 1000    format(a72)                                                    

         read (10,1000) char                                            

c                                                                       

c ------ npo1, npo2, npo3, nelem, nb ----------                         

        read(10,*) npo1, npo2, npo3, nelem, nb                          

c                                                                       

        nb=2*nb                                                         

        nbb=2*nb-1                                                      

        npor=npo3-npo2                                                  

        nom1=2*npo1                                                     

c                                                                       

c -------  input coordinates  -------                                   

c                                                                       

        read(10,1000) char                                              

*        print *, char                                                  

        do 10 i=1, npo3, 2                                              

c       do 10 i=1, npo3, 5                                              

c         read(10,*) (xy(i+j,1), xy(i+j,2), j=0, 5)                     

          read(10,*) (xy(i+j,1), xy(i+j,2), j=0, 1)                     

   10   continue                                                        

c                                                                       

c ------  change unit from mm to m  ------                              

c                                                                       

c check3                                                                

          do 20 i=1, npo3                                               

            xy(i,1)=xy(i,1)*.001                                        

            xy(i,2)=xy(i,2)*.001                                        

   20     continue                                                      

c                                                                       

c ------  input node of elements  nod( , )  -----                       

c                                                                       

        read(10,1000) char                                              

*        print *, char                                                  

        do 30 i=1, nelem, 4                                             

c       do 30 i=1, nelem, 2                                             

c         read(10,*) (nod(i+j,1), nod(i+j,2),                           

c    &                nod(i+j,3), nod(i+j,4),j=0,1)                     

          read(10,*) (nod(i+j,1), nod(i+j,2),                           

     &                nod(i+j,3), nod(i+j,4),j=0,3)                     

   30   continue                                                        

c                                                                       

c ------  input magnetizing voltages  ------                            

c                                                                       

        read(10,1000) char                                              

        read(10,*) ncv                                                  

        do 40 kk=1,ncv                                                  

          read(10,*) ncoil                                              

          ncoil=ncoil-299                                               

          read(10,*) (dv(ncoil,j,1),dv(ncoil,j,2),j=1,nmd)              

   40   continue                                                        

c                                                                       

c  ---------- parameters of windings --------------                     

c                                                                       

        read (10,1000) char                                             

*        print *, char                                                  

        read(10,*) ncc                                                  

        if(ncc.eq.0)          go to 140                                 

        do 80 kk=1,ncc                                                  

          read(10,*) ncoil                                              

          ncoil=ncoil-299                                               

          read (10,*) dnds(ncoil,1),dnds(ncoil,2),dnds(ncoil,3)         

   80  continue                                                         

  140  continue                                                         

c                                                                       

c ------  input magnetizing frequency  --------------                   

c                                                                       

        read(10,1000) char                                              

*        print *, char                                                  

        read(10,*) df                                                   

        domeg=2.d0*dpi*df                                               

c                                                                       

c ------  input conductivity at node( ,4)=400  -----                    

c                                                                       

        read(10,1000) char                                              

*        print *, char                                                  

        read(10,*) nre                                                  

        if(nre.eq.0)            go to 130                               

        do 50 kk=1,nre                                                  

          read(10,*) dre(kk)                                            

   50   continue                                                        

c                                                                       

  130 continue                                                          

c ------  input periodic boundary condition  -------                    

c                                                                       

        read(10,1000) char                                              

        if(npor.eq.0)           go to 110                               

        do 60 kk=1,npor                                                 

          read(10,*) (nodr(kk,j),j=1,2)                                 

   60   continue                                                        

c                                                                       

  110     continue                                                      

c                                                                       

c ------  input data for tableau approach ------                        

c                                                                       

        read(10,1000) char                                              

        read(10,1000) char                                              

        read(10,*) ncn, nbr                                             

        read(10,1000) char                                              

        read(10,1000) char                                              

      do 150 i=1,nbr                                                    

        read(10,*) ndummy,(ntabl(i,j),j=1,3),dtabl(i,1)                 

  150 continue                                                          

c                                                                       

      do 250 i=1, nbr                                                   

        if ((ntabl(i,3).eq.2).or.(ntabl(i,3).eq.3)) then                

          dtabl(i,1)=dtabl(i,1)*domeg                                   

      end if                                                            

  250 continue                                                          

c                                                                       

        ncn=2*ncn                                                       

        nbr=2*nbr                                                       

        nom2=ncn+nbr                                                    

        nom=nom1+nom2                                                   

c                                                                       

c --------  initializated vector potentials daa( , )  ---------         

c                                                                       

c check 31                                                              

c  if nflag=1, then the initial values are set, else are read from files

       nflag=2                                                          

c                                                                       

       if(nflag.eq.2)                  go to 120                        

c                                                                       

       do 70 k=1,nom,2                                                  

         daa(1,k)  =0.0001*(1.0+0.01*float(k-int(k/10))*10.)            

         daa(1,k+1)=0.0                                                 

         do 71 l=1,nmd                                                  

           daa(l,k)  =0.0                                               

           daa(l,k+1)=0.0                                               

   71    continue                                                       

   70  continue                                                         

       go to 300                                                        

c                                                                       

  120  continue                                                         

c                                                                       

c -------  read daa( , )  ------------                                  

c                                                                       

c check 22                                                              

c                            ********* open read file ********          

c     open (20,file='ab0169.hbfem03.data(femref)')                      

c                                                                       

      read(20,1000) char                                                

      read(20,1000) char                                                

      read(20,*)  nhar1                                                 

      nn=(nhar1+1)/2                                                    

      read(20,*)  (ntotal(1,n),n=1,nn)                                  

      read(20,1000) char                                                

      read(20,1000) char                                                

c check char1                                                           

c                                                                       

      do 200 k=1,nom,2                                                  

       nflag=1                                                          

       if (nflag.eq.1)          go to 400                               

c                                                                       

       read(20,*) ndamy, (daa(i,k), i=1,nn)                             

       read(20,*)        (daa(i,k+1), i=1,nn)                           

       go to 410                                                        

c                                                                       

  400 continue                                                          

        if(nhar1.ne.3)                  go to 210                       

c *******  harmonic up to 3-order  ********                             

        read(20,*) ndamy, daa(1,k),daa(1,k+1),daa(2,k),daa(2,k+1)       

        go to 200                                                       

c                                                                       

c *******  harmonic up to 5-order  ********                             

  210 continue                                                          

        if(nhar1.ne.5)  go to 220                                       

c                                                                       

       read(20,*) ndamy, daa(1,k), daa(1,k+1), daa(2,k), daa(2,k+1)     

       read(20,*)        daa(3,k), daa(3,k+1)                           

        go to 200                                                       

c                                                                       

c *******  harmonic up to 7-order  ********                             

  220 continue                                                          

        if(nhar1.ne.7)                  go to 230                       

c                                                                       

       read(20,*) ndamy, daa(1,k), daa(1,k+1), daa(2,k), daa(2,k+1)     

       read(20,*)        daa(3,k), daa(3,k+1), daa(4,k), daa(4,k+1)     

      go to 200                                                         

c                                                                       

c *******  harmonic up to 9-order  ********                             

  230 continue                                                          

        if(nhar1.ne.9)                  go to 240                       

       read(20,*) ndamy, daa(1,k), daa(1,k+1), daa(2,k), daa(2,k+1)     

       read(20,*)        daa(3,k), daa(3,k+1), daa(4,k), daa(4,k+1)     

       read(20,*)        daa(5,k), daa(5,k+1)                           

      go to 200                                                         

c                                                                       

c *******  harmonic up to 11-order  ********                            

  240 continue                                                          

        if(nhar1.ne.11)                 go to 999                       

       read(20,*) ndamy, daa(1,k), daa(1,k+1), daa(2,k), daa(2,k+1)     

       read(20,*)        daa(3,k), daa(3,k+1), daa(4,k), daa(4,k+1)     

       read(20,*)        daa(5,k), daa(5,k+1), daa(6,k), daa(6,k+1)     

c                                                                       

  410 continue                                                          

  200 continue                                                          

*        read(20,1000) char                                             

*        print '(a72)',char                                             

c                                                                       

  300 continue                                                          

c                                                                       

      return                                                            

c                                                                       

  999 print '('' error exits in reading initial values (sub.rdata) '')' 

      return                                                            

      end                                                               

c **************************************************************        

c ***************  save and output data  odata  ****************        

c                                                                       

           subroutine odata                                             

c                                                                       

       include 'common'

c                                                                       

c ********* if double precision ******************                      

c          implicit real*8(d)                                           

c                                                                       

c check 9               **** open write file ****                       

c                                                                       

c          open (30,file='ab0169.hbfem03.data(oda4000)')                

c                                                                       

           write (30, '(a)') '# *** analysis by hbfem //femh0922// ***' 

           write (30, '(a)') '# ** iteration number **'                 

           write (30, '(i6)')  nhar                                     

           write (30, '(6i6)') (ntotal(1,n),n=1,mmu)                    

c                                                                       

c -----  save vector potential  ------                                  

c                                                                       

       write (30, '(a)') '# ** vector potential sin1 cos1 sin3 cos3     

     & sin5  cos5 **'                                                   

       write (30, '(a)') '#                     sin7 cos7 sin9 cos9     

     & sin11 cos11  '                                                   

c                                                                       

       do 100 n=1,nom,2                                                 

         if(nhar.eq.3) then                                             

           write (30,'(i4,4e14.5)') (n+1)/2, daa(1,n), daa(1,n+1),      

     &                                       daa(2,n), daa(2,n+1)       

         else if (nhar.eq.5) then                                       

           write (30,'(i4,4e14.5)') (n+1)/2, daa(1,n), daa(1,n+1),      

     &                                       daa(2,n), daa(2,n+1)       

           write (30,'(4x,2e14.5)')          daa(3,n), daa(3,n+1)       

         else if (nhar.eq.7) then                                       

           write (30,'(i4,4e14.5)') (n+1)/2, daa(1,n), daa(1,n+1),      

     &                                       daa(2,n), daa(2,n+1)       

           write (30,'(4x,4e14.5)')          daa(3,n), daa(3,n+1),      

     &                                       daa(4,n), daa(4,n+1)       

         else if (nhar.eq.9) then                                       

           write (30,'(i4,4e14.5)') (n+1)/2, daa(1,n), daa(1,n+1),      

     &                                       daa(2,n), daa(2,n+1)       

           write (30,'(4x,4e14.5)')          daa(3,n), daa(3,n+1),      

     &                                       daa(4,n), daa(4,n+1)       

           write (30,'(4x,2e14.5)')          daa(5,n), daa(5,n+1)       

         else                                                           

           write (30,'(i4,4e14.5)') (n+1)/2, daa(1,n), daa(1,n+1),      

     &                                       daa(2,n), daa(2,n+1)       

           write (30,'(4x,4e14.5)')          daa(3,n), daa(3,n+1),      

     &                                       daa(4,n), daa(4,n+1)       

           write (30,'(4x,4e14.5)')          daa(5,n), daa(5,n+1),      

     &                                       daa(6,n), daa(6,n+1)       

         end if                                                         

c                                                                       

  100    continue                                                       

c                                                                       

c -------  save flux density dbx( , ) & dby( , )  -------               

c                                                                       

c       nflag=0                                                         

        nflag=1                                                         

c                                                                       

        if (nflag.eq.0) goto 310                                        

c                                                                       

       write (30, '(a)') '# ** flux density, sin1-x cos1-x sin3-x... **'

       write (30, '(a)') '#                  sin1-y cos1-y sin3-y... '  

c                                                                       

       do 200 i=1,nelem                                                 

         write (30, '(i4,6e12.4)') i,(dbx(i,k,1), k=1,mmu)              

         write (30, '(4x,6e12.4)')   (dbx(i,k,2), k=1,mmu)              

         write (30, '(4x,6e12.4)')   (dby(i,k,1), k=1,mmu)              

         write (30, '(4x,6e12.4)')   (dby(i,k,2), k=1,mmu)              

  200  continue                                                         

c                                                                       

  310  continue                                                         

      write(30, '(a)')   '#  *****  comment  *****  '                   

      write(30, '(3e14.4)')  aaa, bbb, df                               

      write(30, '(7i6)') nelem,nb,nbb,nre,ncc,ncn,nbr                   

      write(30, '(6e9.2)') (dv(1,j,1),dv(1,j,2),j=1,nmd)                

      write(30, '(10e14.4)') (dre(i),i=1,nre)                           

c                                                                       

       return                                                           

       end                                                              


c ***************************************************************       

c ************  arrange the system matrix  fem  *****************       

c                                                                       

           subroutine fem                                               

c                                                                       

       include 'common'

c                                                                       

         dimension ddd(nmd,2,2), dss(3,3)                               

c                                                                       

c **********  loop for element  ****************                        

c                                                                       

      do 400 i=1, nom2                                                  

        dk(nom1+i)=dkk(nhh,nom1+i)                                      

  400 continue                                                          

c                                                                       

      do 100 ne=1, nelem                                                

        if(int(nod(ne,4)/100).ne.3)   go to 190                         

c  *** coil region ****                                                 

        if(nod(ne,4).ge.350) then                                       

          nee=nod(ne,4)-349                                             

        else                                                            

          nee=nod(ne,4)-299                                             

        end if                                                          

c                                                                       

        jj1=ncn+2*(nee-1)                                               

        jj2=2*(nee-1)                                                   

        do 110 i=1,3                                                    

          ii=nod(ne,i)                                                  

          if ((ii.gt.npo1).and.(ii.le.npo2))   go to 110                

c                                                                       

          if (ii.le.npo1) then                                          

            nlow=2*(ii-1)                                               

            dsig1= 1.0                                                  

          else                                                          

            ii=nodr(ii-npo2,2)                                          

            nlow=2*(ii-1)                                               

            dsig1= -1.0                                                 

          end if                                                        

c                                                                       

          do 120 j=1,2                                                  

            dhc(nlow+j,jj1+j)=dhc(nlow+j,jj1+j)-dg(ne)*dsig1            

  120     continue                                                      

c                                                                       

          do 140 j=1,2                                                  

            do 150 k=1,2                                                

              if(dn(nhh,j,k).eq.0.0)                    go to 150       

              dhl(jj2+j,nlow+k)=dhl(jj2+j,nlow+k)                       

     &                         -dn(nhh,j,k)*dc(ne)*dsig1*2.             

  150       continue                                                    

  140     continue                                                      

  110   continue                                                        

c                                                                       

  190   continue                                                        

c                                                                       

c ---------  set ddd(mmu,2,2)  ----------                               

c                                                                       

       do 210 i=1, 2                                                    

         do 220 j=1, 2                                                  

           do 230 k=1, mmu                                              

             ddd(k,i,j)=dd(ne,k,i,j)                                    

  230      continue                                                     

  220    continue                                                       

  210  continue                                                         

c                                                                       

c --------  dss( , )  -------------                                     

c                                                                       

        nss=1                                                           

        do 240 j=1,3                                                    

          do 250 k=j,3                                                  

            dss(j,k)=ds(ne,nss)                                         

            dss(k,j)=ds(ne,nss)                                         

            nss=nss+1                                                   

  250     continue                                                      

  240   continue                                                        

c                                                                       

c -------  insert dss(,) & ddd(,) into dh#(,)  ------------             

c                                                                       

        do 300 j=1,3                                                    

          ii=nod(ne,j)                                                  

          if ((ii.gt.npo1).and.(ii.le.npo2))   go to 300                

c                                                                       

          if (ii.le.npo1) then                                          

            nlow=2*(ii-1)+1                                             

            dsig1= 1.0                                                  

          else                                                          

            ii=nodr(ii-npo2,2)                                          

            nlow=2*(ii-1)+1                                             

            dsig1= -1.0                                                 

          end if                                                        

c                                                                       

c --------  scalor potential  ----------                                

c                                                                       

c check                                                                 

c             if (int(nod(ne,4)/100).ne.4)          go to 340           

c           ie=nod(ne,4)-399                                            

c           dcons=dre(ie)*dcs(ne)*(2.0*float(j1)-1.0)*domeg*dsig1/3.0d0 

c           dk(nlow)  =dk(nlow)  -dcc(ie,j1,1)*dcons                    

c           dk(nlow+1)=dk(nlow+1)-dcc(ie,j1,2)*dcons                    

c  340  continue                                                        

c                                                                       

            do 310 k=1,3                                                

              if(k.eq.j) then                                           

                dnn=2.0d0                                               

              else                                                      

                dnn=1.0d0                                               

              end if                                                    

c                                                                       

            jj=nod(ne,k)                                                

            if ((jj.gt.npo1).and.(jj.le.npo2))     go to 310            

            if (jj.le.npo1) then                                        

              ncol=2*(jj-1)+1                                           

              dsig2 = 1.0                                               

            else                                                        

              jj=nodr(jj-npo2,2)                                        

              ncol=2*(jj-1)+1                                           

              dsig2 = -1.0                                              

            end if                                                      

c                                                                       

            dnc1=dss(j,k)*dsig1*dsig2                                   

            dnc2=dnn*dsig1*dsig2                                        

            do 320 m=1, 2                                               

              nl=nlow+m-1                                               

              do 330 n=1, 2                                             

                nc =ncol+n-1                                            

                nc1=nc-nl+nb                                            

c                                                                       

                dh(nl,nc1) =dh(nl,nc1)  + ddd(1,m,n)*dnc1               

     &                                  + dn(nhh,m,n)*ds(ne,7)*dnc2     

                do 330 l=2,mmu                                          

                  dk(nl)=dk(nl)-ddd(l,m,n)*dnc1*                        

     &                          daa(mod(l+nhh-2,mmu)+1,nc)              

  330           continue                                                

  320       continue                                                    

  310     continue                                                      

  300   continue                                                        

  100 continue                                                          

c                                                                       

      return                                                            

      end                                                               
                  program fh0922                                       

c *************************************************************         

c *************************************************************         

c ***((( fort77 )))                                         ***         

c ***  time-periodic nonlinear magnetic field analysis      ***         

c ***    by harmonic balance finite element method          ***         

c ***          femh0922              ver.3.3                ***         

c ***                                                       ***         

c ***  1. generalized harmonic balance fem                  ***         

c ***  2. harmonic up to 11th can be considered             ***         

c ***  3. the maximu harmonics is automatically defined     ***         

c ***  4. method of block-structure iteration is employed   ***         

c ***  5. the external circuit is considered                ***         

c ***                                                       ***         

c ***                                                       ***         

c ***                                       91.9.10         ***         

c *************************************************************         

c *************************************************************         

c   fh0922.f    (main program)                                          

c   pripro     pripro  init                                             

c   data.f     rdata odata                                              

c   fem.f      smat  dmat  cmat  fem                                    

c   gauss3.f   gauss3                                                   

c   flux.f     flux                                                     

c   rfour.f    rfour                                                    

c   postpr.f   dstore  conv  postpr  ckeck                              

c                                                                       

c---------------- list of variable ----------------------------         

c                                                                       

c    aaa            criterion of conversion                             

c    bbb            decelerating constant                               

c    nonc           status on return                                    

c    nhar           numbers of harmonic to consider (3,5,7,9,11)        

c    mmu            number of harmonic to be considered (mmu=(nhar+1)/2)

c    npo1           number of unknown potentials (degree of freedom)    

c    nom            size of system matrix                               

c    nom1           nom1=mmu*npo1                                       

c    nom2           nom2=mmu*(ncn+nbr)                                  

c    npo2           last number of node on dirichlet's condition        

c    npo3           number of node                                      

c    npor           count of node on periodic boundary                  

c                      (npor=npo3-npo2)                                 

c    ncn            number of circuit node                              

c    nbr            number of branch                                    

c    ncc            number of coil                                      

c                                                                       

c    df             frequency                                           

c ***** definition of node number *****                                 

c 1+++(unknown npo1)+++(dirichlet npo2)+++(periodic npo3)               

c                                                                       

c    nelem          number of element                                   

c    nb             bound width of system matrix                        

c    nod( , )       node numner of triangle element                     

c                   and material                                        

c                                                                       

c -----------------------------------------------------------------     

c    dh( , )        system matrix  h                                    

c    dhc( , )       system matrix  hc                                   

c    dhl( , )       system matrix  hl                                   

c    dhcl( , )      system matrix  hcl                                  

c    dk( )          forced column vector for calculation                

c    dkk( )         ibid for memory                                     

c    da( )          calculated values of vector potential               

c    daa( )         initial values of vector potential                  

c    daaa( )        temporary vector potential                          

c    ds( ,1-6)      matrix for on eelement                              

c    ds( ,7)        coefficient for eddy current                        

c    dcs( )         cross section of element                            

c    dbx( , , )     flux density in x-direction                         

c    dby( , , )     flux density in y-direction                         

c    dd( , )        reluctivity matrix at each elemnent                 

c    xy( , )        coordinates of node                                 

c    dr( , )        conductivity                                        

c    dcur( , , )    current density (sin1 cos1 sin3 cos3                

c                                    sin5 cos5 sin7 cos7 ...)           

c    nodr( , )      corresponding node on periodic boundary             

c    dn( , )        harmonic matrix                                     

c    dcc( , , )     scalor potential                                    

c    dssum( )       sum of conductor in which eddy currents flows       

c    neddy( , )     node # and element # in the conductor               

c    deddy( , )     eddy-current density at each node                   

c    dnds(nnbr,3)   parameters of coils                                 

c                                                                       

c    dzz( , )       impedance matrix                                    

c    dyy( , )       admittance matrix                                   

c    das( , )       matrix                                              

c    ntabl( , )     data for tableau approach                           

c    dtabl( , )     data for branch element                             

c *******************************************************************   

c ******************** main program *********************************   

c                                                                       

c    nfd            freedom degree of matrix                            

c    ned            number of element                                   

c    nmd            number of harmonic to be considered                 

c                    (nmd(=4) means 1,3,5 and 7-order)                  

c    nbd            band width                                          

c    nxy            number of node                                      

c    ndi            number of divided area for simpson's integration    

c    nncn           number of circuit node for matrix size              

c    nnbr           number of branch for matrix size                    

c                                                                       

      include 'common'

c                                                                       

c ************ if double precision *********************                

c        implicit  real*8(d)                                            

c                                                                       

c check1                                                                

c ************ set condition **********************                     

c                                                                       

c   1. maximum ineration            max                                 

                                    max=150                             

c   2. criterion of convergence     aaa                                 

                                    aaa=0.005                           

c   3. deceleration constant        bbb                                 

                                    bbb=0.2                             

c   4. number of harmonic to be considered (3,5,7,9,11)                 

c                     maximum       mmu < mnd                           

                                    mmu=4                               

                                    nhar=2*mmu-1                        

c                                                                       

c *********** set coeficients ************                              

c   1. pi                           dpi                                 

                                    dpi=3.14159265358979                

c   2. permeability of air          dmu0                                

           dmu0=4.0d-7*dpi                                              

c   3. angular frequency            domeg                               

c          domeg=2.0d0*dpi*df                                           

c      which is calculated in subroutine ' rdata '                      

c                                                                       

c *********** 1. pri program *************                              

c                                                                       

c 1-1 read data                                                         

*       print '('' rdata '')'                                           

              call rdata                                                

c                                                                       

c 1-2 set coefficients and initial values                               

*       print '('' pripro '')'                                          

              call pripro                                               

c                                                                       

c 1-3 calculate matrix for each element                                 

*       print '('' smat '')'                                            

              call smat                                                 

c                                                                       

c set miscellaneous coefficients                                        

c     initial harmonic set number of iterations                         

      nsum=0                                                            

c                                                                       

c ##########  2. loop of calculation  ##########                        

c                                                                       

  100       continue                                                    

c                                                                       

         nhar=2*mmu-1                                                   

c                                                                       

      do 200 nloop=1, mmu                                               

c                                                                       

c ******* loop for harmonic *******                                     

c                                                                       

         itr=0                                                          

         nhh=nloop                                                      

c                                                                       

  210 continue                                                          

         itr=itr+1                                                      

         nsum=nsum+1                                                    

c                                                                       

         if ((itr.gt.max).or.(nsum.gt.max)) then                        

               print '(''stop because iteration is over'',i4)',max      

               goto 9999                                                

         end if                                                         

c                                                                       

c 2-1 clear variables                                                   

*      print '('' init '')'                                             

          call init                                                     

c                                                                       

c 2-2 calculate flux density                                            

*      print '('' flux '')'                                             

          call flux                                                     

c                                                                       

c 2-3 calculate reluctivity matrix                                      

*      print '('' dmat '')'                                             

          call dmat                                                     

c                                                                       

c 2-4-1 set system matrix                                               

*      print '('' fem '',i2)',nloop                                     

          call fem                                                      

c                                                                       

c 2-4-2 set system matrix                                               

*      print '('' cmat '',i2)',nloop                                    

          call cmat                                                     

c                                                                       

c 2-5 calculate system equation                                         

*      print '('' gauss3 '',i2)',nloop                                  

          call gauss3                                                   

c                                                                       

c 2-6 store each harmonic data                                          

*      print '('' dstore '',i2)',nloop                                  

          call dstore                                                   

c                                                                       

c 2-7 check convergence                                                 

*      print '('' conv'',i2)', nloop                                    

            call conv                                                   

c                                                                       

       print '(''number of non-conv, point '',i4,3x,3i4)',              

     &              nonc,   mmu,nloop,itr                               

        if (nonc.ne.0)          goto 210                                

c                                                                       

      ntotal(2,nloop)=itr                                               

      ntotal(1,nloop)=ntotal(1,nloop)+ntotal(2,nloop)                   

c                                                                       

  200 continue                                                          

c                                                                       

c ########## check the convergence of all harmonics ##########          

c                                                                       

      nonc=0                                                            

      do 120 i=1, mmu                                                   

        if (ntotal(2,i).ne.1)     nonc=nonc+1                           

  120 continue                                                          

c                                                                       

      do 130 i=1, mmu                                                   

        ntotal(2,i)=0                                                   

  130 continue                                                          

c                                                                       

           if (nonc.ne.0)   goto 100                                    

c                                                                       

c ##### number of non-convergend potential #####                        

c                                                                       

c 3-1 check  convergence                                                

*      print '('' check  '')'                                           

          call check                                                    

c                                                                       

*      if (nonc.eq.0)        goto 300                                   

*        if (mmu.ge.nmd) then                                           

*           print '(''stop because over number of harmonic'')'          

*         goto 9999                                                     

*      end if                                                           

c                                                                       

*          mmu=mmu+1                                                    

*          goto 100                                                     

c                                                                       

c ################# loop end ##################                         

c                                                                       

* 300   continue                                                        

c ********** when converged **********                                  

c                                                                       

        print '(''*******************'')'                               

        print '(''*** convergence ***'')'                               

        print '(''*******************'')'                               

c                                                                       

c ******* 3. output calculation results *****                           

c                                                                       

 9999        continue                                                   

c                                                                       

c 3-1 calculate data                                                    

*      print '('' postpr '')'                                           

              call postpr                                               

c                                                                       

c 3-2 print out data                                                    

*      print '('' odata '')'                                            

              call odata                                                

      stop                                                              

      end                                                               


c ******  calculate flux density for each element  flux  ******         

c                                                                       

           subroutine flux                                              

c                                                                       

       include 'common'

c                                                                       

       dimension dx(3),dy(3),dq(3),dr(3),dda(nmd,3,2)                   

c                                                                       

      do 100 ne=1,nelem                                                 

        do 110 j=1,3                                                    

          n=nod(ne,j)                                                   

          dx(j)=xy(n,1)                                                 

          dy(j)=xy(n,2)                                                 

          do 120 k=1, mmu                                               

            if (n.le.npo1)            then                              

              dda(k,j,1)=daa(k,2*n-1)                                   

              dda(k,j,2)=daa(k,2*n)                                     

            else if (n.le.npo2)       then                              

              dda(k,j,1)=0.0                                            

              dda(k,j,2)=0.0                                            

            else if (n.le.npo3)       then                              

              n2=nodr(n-npo2,2)                                         

              dda(k,j,1)=daa(k,2*n2-1)                                  

              dda(k,j,2)=daa(k,2*n2)                                    

            else                                                        

              print '('' error in data of element, prog(flux)'',2i3)',  

     &       ne,n                                                       

            end if                                                      

  120     continue                                                      

  110   continue                                                        

c                                                                       

        dq(1)=dy(2) -dy(3)                                              

        dr(1)=dx(3) -dx(2)                                              

        dq(2)=dy(3) -dy(1)                                              

        dr(2)=dx(1) -dx(3)                                              

        dq(3)=dy(1) -dy(2)                                              

        dr(3)=dx(2) -dx(1)                                              

c                                                                       

        do 140 j=1,mmu                                                  

          dbx(ne,j,1)=0.0                                               

          dbx(ne,j,2)=0.0                                               

          dby(ne,j,1)=0.0                                               

          dby(ne,j,2)=0.0                                               

c                                                                       

          do 150 k=1,3                                                  

            dbx(ne,j,1)=dbx(ne,j,1)+dr(k)*dda(j,k,1)                    

            dbx(ne,j,2)=dbx(ne,j,2)+dr(k)*dda(j,k,2)                    

            dby(ne,j,1)=dby(ne,j,1)-dq(k)*dda(j,k,1)                    

            dby(ne,j,2)=dby(ne,j,2)-dq(k)*dda(j,k,2)                    

  150     continue                                                      

          dbx(ne,j,1)=dbx(ne,j,1)*.5/dcs(ne)                            

          dbx(ne,j,2)=dbx(ne,j,2)*.5/dcs(ne)                            

          dby(ne,j,1)=dby(ne,j,1)*.5/dcs(ne)                            

          dby(ne,j,2)=dby(ne,j,2)*.5/dcs(ne)                            

  140   continue                                                        

c                                                                       

  100 continue                                                          

      return                                                            

      end                                                               

c **************************************************************        

c *********  fourier expansion of reluctivity, rfour ***********        

c                                                                       

c            ne : number of element                                     

c                                                                       

           subroutine rfour(ne)                                         

c                                                                       

       include 'common'

c                                                                       

           dimension dab(0:144)                                         

c                                                                       

c ----  set parameters -----                                            

c                                                                       

           dcon=1.0/dpi                                                 

c                                                                       

c ****  ndi points simpson's integral  ****                             

c                                                                       

         ddh=dpi/float(ndi)                                             

         do 100 nt=0,ndi                                                

c                                                                       

c check 6                                                               

c                                                                       

         dbbx=0.0                                                       

         dbby=0.0                                                       

         do 110 i=1, mmu                                                

           dbbx = dbbx+dbx(ne,i,1)*dsine(2*i-1,nt)                      

     &                +dbx(ne,i,2)*dcosi(2*i-1,nt)                      

           dbby = dbby+dby(ne,i,1)*dsine(2*i-1,nt)                      

     &                +dby(ne,i,2)*dcosi(2*i-1,nt)                      

  110    continue                                                       

           dbt=sqrt(dbbx*dbbx+dbby*dbby)                                

c                                                                       

c check 7                                                               

c -------  dbx/dt, dby/dt  -------                                      

c                                                                       

         nflug=1                                                        

         if (nflug.eq.1)                         go to 130              

c                                                                       

          dbbxt=0.0                                                     

          dbbyt=0.0                                                     

          do 120 i=1,mmu                                                

          dbbxt = dbbxt+(dbx(ne,i,1)*dcosi(2*i-1,nt)                    

     &                 - dbx(ne,i,2)*dsine(2*i-1,nt))*(2*i-1)*domeg     

          dbbyt = dbbyt+(dby(ne,i,1)*dcosi(2*i-1,nt)                    

     &                 - dby(ne,i,2)*dsine(2*i-1,nt))*(2*i-1)*domeg     

  120      continue                                                     

c                                                                       

  130      continue                                                     

c                                                                       

c **************************************************************        

c *******  magnetizing characteristic of magnetic core  ********        

c **************************************************************        

c                                                                       

c  dbt          flux density                    b                       

c  dbxt         flux density                    bx                      

c  dbyt         flux density                    by                      

c  dbbxt        derivative of flux density      dbx/dt                  

c  dbbyt        derivative of flux density      dby/dt                  

c                                                                       

c -----  select magnetizing characteristics  -------                    

c    if characteristic i,then go to 9#0                                 

         go to 940                                                      

c                                                                       

  910  continue                                                         

c                                                                       

c---------------------(1)-----------------------------------            

c                                                                       

c 1. hysteresis model                                                   

c                                                                       

       dab(nt)=128.+100.*dbt**8                                         

c     &        +(0.1+0.01*dbt**2)*(dbbx*dbbxt+dbby*dbbyt)/(dbt*dbt)     

c                                                                       

       go to 999                                                        

c                                                                       

c---------------------(2)-----------------------------------            

c                                                                       

c 2.silicon steel characteristic                                        

c                                                                       

  920       continue                                                    

c                                                                       

       dab(nt)=100. +40.4*dbt**8                                        

c                                                                       

       go to 999                                                        

c                                                                       

c---------------------(3)------------------------------------           

c                                                                       

c 3. ferrite core     h7c1 tdk                                          

c                                                                       

  930       continue                                                    

c                                                                       

       if (abs(dbt).gt.0.4)                      go to 931              

         dab(nt)=84.55+1582.4*dbt**4                                    

         go to 938                                                      

  931  if (abs(dbt).gt.0.435)                    go to 932              

         dab(nt)=50.+1428.*(abs(dbt)-0.4)                               

         go to 939                                                      

  932  if (abs(dbt).gt.0.460)                    go to 933              

         dab(nt)=100.+4.0e3*(abs(dbt)-0.435)                            

         go to 939                                                      

  933  if(abs(dbt).gt.0.490)                     go to 934              

         dab(nt)=200.+1.0e4*(abs(dbt)-0.460)                            

         go to 939                                                      

  934  if (abs(dbt).gt.0.510)                    go to 935              

         dab(nt)=500.+2.0e4*(abs(dbt)-0.490)                            

         go to 939                                                      

  935  if (abs(dbt).gt.0.520)                    go to 936              

         dab(nt)=900.+6.0e4*(abs(dbt)-0.510)                            

         go to 939                                                      

  936    dab(nt)=1500.+1.0e5*(abs(dbt)-0.520)                           

c                                                                       

  939  dab(nt)=dab(nt)/abs(dbt)                                         

  938  continue                                                         

c                                                                       

       go to 999                                                        

c                                                                       

c---------------------(4)-----------------------------------            

c                                                                       

c 4. silicon steel characteristic                                       

c                                                                       

  940       continue                                                    

c                                                                       

       dab(nt)=100. +40.*dbt**8                                         

c      dab(nt)=100.                                                     

c                                                                       

       go to 999                                                        

c                                                                       

c*********************                                                  

  999       continue                                                    

c ****** end of magnetization curve representation *****                

c                                                                       

  100  continue                                                         

c                                                                       

c -----  fourier expansion by simpson integral method  ---------        

c                                                                       

         n1=ndi-1                                                       

         dsc=dab(0)                                                     

         do 200 n=1,n1,2                                                

           dsc=dsc+4.*dab(n)+2.*dab(n+1)                                

  200    continue                                                       

c   **** 2*dreluc(0,1) ****                                             

         dreluc(0,1)=2.*(dsc-dab(ndi))*ddh*dcon/3.                      

         dreluc(0,0)=0.0                                                

c                                                                       

c check 8                                                               

c                                                                       

        nend=2*(mmu+nhh-1)                                              

         do 210 n=2, nend, 2                                            

           dsc=dab(0)                                                   

           dss=0.0                                                      

c                                                                       

           do 220 nk=1, n1, 2                                           

             dys=dab(nk)*dsine(n,nk)                                    

               dss=dss+4.*dys                                           

             dys=dab(nk+1)*dsine(n,nk+1)                                

               dss=dss+2.*dys                                           

             dyc=dab(nk)*dcosi(n,nk)                                    

               dsc=dsc+4.*dyc                                           

             dyc=dab(nk+1)*dcosi(n,nk+1)                                

               dsc=dsc+2.*dyc                                           

  220      continue                                                     

           dreluc(n,0)=(dss-dys)*ddh*dcon*2.0/3.0                       

           dreluc(n,1)=(dsc-dyc)*ddh*dcon*2.0/3.0                       

  210    continue                                                       

c                                                                       

       return                                                           

       end         
c ***************************************************************       

c ********  gaussian elemination method, gauss  ver.3.2 *********       

c                                                                       

      subroutine gauss3                                                 

c                                                                       

       include 'common'

c                                                                       

c  number of matrix width                                               

      ndia=nb                                                           

c                                                                       

c *****************************************************************     

c     dk(nom)        : forced vector                                    

c     dhc(nom1,nom2) : side matrix of dh matrix                         

c     dhl(nom2,nom1) : under matrix of dh matrix                        

c     dhcl(nom2,nom2): coner matrix of dh matrix                        

c     daa(nom)       : potential      (nob)                             

c     nom            : number of unknown potential                      

c              nom=nom1+nom2                                            

c     nom1           :                                                  

c     nom2           :                                                  

c     ndia           :                                                  

c              ndia=nb                                                  

c     nbb            : band width                                       

c              nbb=2*nb-1                                               

c                                                                       

c ************** eliminate l of matrix dh(i,j) ********************     

c                                                                       

              nbb=2*nb-1                                                

c                                                                       

      do 100 nn=1, nom-1                                                

c                                                                       

        if (nn.gt.nom1)                      go to 110                  

        dt1=1.0/dh(nn,ndia)                                             

c ***** dh ******                                                       

c (nm,nn)                                                               

        do 200 nm=nn+1, nom                                             

          if(nm.gt.nom1)                     go to 210                  

c ***** dh-dh ******                                                    

          j1=ndia-nm+nn                                                 

          i1=nm                                                         

          if(j1.lt.1)                        go to 200                  

c                                                                       

          if(dh(i1,j1).eq.0.0)               go to 200                  

          dt=dt1*dh(i1,j1)                                              

c                                                                       

          dk(i1)=dk(i1)-dt*dk(nn)                                       

c (nm,nk)                                                               

          ii1=nm                                                        

          ii2=nn                                                        

          do 300 nk=nn, nom                                             

c                                                                       

            if(nk.gt.nom1)                   go to 310                  

c  ****** dh-dh *****                                                   

            jj1=ndia-nm+nk                                              

            jj2=ndia-nn+nk                                              

            if((jj1.gt.nbb).or.(jj1.lt.1))    go to 300                 

            if((jj2.gt.nbb).or.(jj2.lt.1))    go to 300                 

c                                                                       

            dh(ii1,jj1)=dh(ii1,jj1)-dt*dh(ii2,jj2)                      

            go to 300                                                   

c   ***** hdc-hdc *****                                                 

  310       continue                                                    

            jj1=nk-nom1                                                 

            jj2=jj1                                                     

            dhc(ii1,jj1)=dhc(ii1,jj1)-dt*dhc(ii2,jj2)                   

c      -- nk --                                                         

  300     continue                                                      

c                                                                       

        go to 200                                                       

c                                                                       

c   ***** dh-dhl, dhc-dhlc *****                                        

  210     continue                                                      

            ii1=nm-nom1                                                 

            ii2=nn                                                      

            jj1=nn                                                      

c                                                                       

            if(dhl(ii1,jj1).eq.0.0)           go to 200                 

            dt=dt1*dhl(ii1,jj1)                                         

c                                                                       

            i1=ii1+nom1                                                 

            dk(i1)=dk(i1)-dt*dk(ii2)                                    

c (nm,nk)                                                               

          do 400 nk=nn, nom                                             

            if(nk.gt.nom1)                    go to 410                 

c  ****** dh-dhl *****                                                  

            jj1=nk                                                      

            jj2=ndia-nn+nk                                              

            if(jj2.gt.nbb)                    go to 400                 

c                                                                       

            if(dh(ii2,jj2).eq.0.0)            go to 400                 

            dhl(ii1,jj1)=dhl(ii1,jj1)-dt*dh(ii2,jj2)                    

            go to 400                                                   

c   ***** hdc-hdcl *****                                                

  410     continue                                                      

            jj1=nk-nom1                                                 

            jj2=jj1                                                     

            if(dhc(ii2,jj2).eq.0.0)           go to 400                 

            dhcl(ii1,jj1)=dhcl(ii1,jj1)-dt*dhc(ii2,jj2)                 

c      -- nk --                                                         

  400     continue                                                      

c      -- nm --                                                         

  200   continue                                                        

c                                                                       

        go to 100                                                       

c                                                                       

c ****** dhlc *******                                                   

  110   continue                                                        

c (nm,nn)                                                               

        nnn=nn-nom1                                                     

c ***** pipotting *****                                                 

        l=nnn                                                           

        al=abs(dhcl(nnn,nnn))                                           

        do 20 j=nnn+1, nom2                                             

          if(abs(dhcl(j,nnn)).gt.al) then                               

            l=j                                                         

            al=abs(dhcl(j,nnn))                                         

          end if                                                        

   20   continue                                                        

c                                                                       

        if (nnn.ne.l) then                                              

          do 30 k=nnn,nom2                                              

            dtmp=dhcl(nnn,k)                                            

            dhcl(nnn,k)=dhcl(l,k)                                       

            dhcl(l,k)=dtmp                                              

   30     continue                                                      

          dtmp=dk(nom1+nnn)                                             

          dk(nom1+nnn)=dk(nom1+l)                                       

          dk(nom1+l)=dtmp                                               

        end if                                                          

c       print '( i4, e14.4)',nnn,dhcl(nnn,nnn)                          

c       if(abs(dhcl(nnn,nnn)).gt.1.0e-20)  go to 1110                   

c       print '(''bunnbo zero'',i4)',nnn                                

 1110   continue                                                        

c                                                                       

        dt1=1.0/dhcl(nnn,nnn)                                           

        do 500 nm=nn+1,nom                                              

c                                                                       

          jj2=nn-nom1                                                   

          ii2=nm-nom1                                                   

          if(dhcl(ii2,jj2).eq.0.0)                go to 500             

          dt=dt1*dhcl(ii2,jj2)                                          

c                                                                       

          ii1=nm                                                        

          dk(ii1)=dk(ii1)-dt*dk(nn)                                     

c (nm,nk)                                                               

            ii1=nm-nom1                                                 

            ii2=nn-nom1                                                 

          do 600 nk=nn, nom                                             

            jj1=nk-nom1                                                 

            jj2=jj1                                                     

c                                                                       

            dhcl(ii1,jj1)=dhcl(ii1,jj1)-dt*dhcl(ii2,jj2)                

c      -- nk --                                                         

  600     continue                                                      

c      -- nm --                                                         

  500   continue                                                        

c      -- nn --                                                         

  100 continue                                                          

c                                                                       

c ******************* calculate da(i) *******************               

c                                                                       

      da(nom)=dk(nom)/dhcl(nom2,nom2)                                   

c                                                                       

      do 800 nn=1, nom-1                                                

        ii1=nom-nn                                                      

        if(ii1.le.nom1)                          go to 810              

c                                                                       

        ii2=ii1-nom1                                                    

        do 700 nk=ii2+1, nom2                                           

          jj2=nk+nom1                                                   

          dk(ii1)=dk(ii1)-dhcl(ii2,nk)*da(jj2)                          

  700 continue                                                          

        da(ii1)=dk(ii1)/dhcl(ii2,ii2)                                   

        go to 800                                                       

c                                                                       

  810   continue                                                        

        do 900 nk=ii1+1, nom                                            

          if(nk.gt.nom1)                         go to 910              

c  ***** dh *****                                                       

          jj1=ndia-ii1+nk                                               

          if(jj1.gt.nbb)                          go to 900             

          dk(ii1)=dk(ii1)-dh(ii1,jj1)*da(nk)                            

          go to 900                                                     

c                                                                       

  910   continue                                                        

          jj1=nk-nom1                                                   

          dk(ii1)=dk(ii1)-dhc(ii1,jj1)*da(nk)                           

c   -- nk --                                                            

  900 continue                                                          

         if(dh(ii1,ndia).eq.0.0)       then                             

            print '('' matrix error (3)'',i5)',nn                       

         end if                                                         

         da(ii1)=dk(ii1)/dh(ii1,ndia)                                   

c   -- nn --                                                            

  800 continue                                                          

       return                                                           

       end     
                                                         
                                                     
c *******************************************************************   

c ********** calculate matrix for each element  smat ****************   

c                                                                       

           subroutine smat                                              

c                                                                       

       include 'common'

c                                                                       

       dimension dx(3),dy(3),dq(3),dr(3)                                

c                                                                       

c ********* if double precision ******************                      

c          implicit real*8(d)                                           

c                                                                       

c ----- calculate the matrix ds( , ) ------                             

      do 100 ne=1,nelem                                                 

        do 110 j=1,3                                                    

          dx(j)=xy(nod(ne,j),1)                                         

          dy(j)=xy(nod(ne,j),2)                                         

  110   continue                                                        

c ---- calculate the cross section at at each element dcs( ) ----       

        dcs(ne)=(dx(2)*dy(3)+dx(1)*dy(2)+dx(3)*dy(1)                    

     &          -dx(2)*dy(1)-dx(1)*dy(3)-dx(3)*dy(2))*.5d0              

c ---------------------------------------------------------------       

         dq(1)=dy(2) -dy(3)                                             

         dr(1)=dx(3) -dx(2)                                             

         dq(2)=dy(3) -dy(1)                                             

         dr(2)=dx(1) -dx(3)                                             

         dq(3)=dy(1) -dy(2)                                             

         dr(3)=dx(2) -dx(1)                                             

c ---------------------------------------------------------------       

        n=1                                                             

        do 120 j=1,3                                                    

          do 130 k=j,3                                                  

            ds(ne,n)=0.25d0/dcs(ne)*(dr(j)*dr(k)+dq(j)*dq(k))           

            n=n+1                                                       

  130     continue                                                      

  120   continue                                                        

c ------------  calculate dg( ),dc( )  -----------------                

c                                                                       

c       do 150 j=1,3                                                    

          if(int(nod(ne,4)/100).eq.3)      then                         

            if (nod(ne,4).ge.350) then                                  

              nee=nod(ne,4)-349                                         

              dg(ne)=-dcs(ne)*.3333333/dnds(nee,3)                      

              dc(ne)=domeg*dg(ne)*dnds(nee,2)                           

            else                                                        

              nee=nod(ne,4)-299                                         

              dg(ne)=dcs(ne)*.3333333/dnds(nee,3)                       

              dc(ne)=domeg*dg(ne)*dnds(nee,2)                           

            end if                                                      

          else                                                          

            dg(ne)=0.0                                                  

            dc(ne)=0.0                                                  

          end if                                                        

c 150   continue                                                        

c                                                                       

c -------  calculate coefficients for eddy currents  -----              

c                                                                       

        if (int(nod(ne,4)/100).eq.4) then                               

          ds(ne,7)=dre(nod(ne,4)-399)*domeg*dcs(ne)/12.d0               

        else                                                            

          ds(ne,7)=0.0                                                  

        end if                                                          

c                                                                       

  100   continue                                                        

c                                                                       

        do 200 i=1,nbr/2                                                

          do 210 j=1,nmd                                                

            jj1=nom1+2*i-1                                              

            dkk(j,jj1)=dv(i,j,1)                                        

            dkk(j,jj1+1)=dv(i,j,2)                                      

  210     continue                                                      

  200   continue                                                        

        return                                                          

          end                                                           

c *******************************************************************   

c ********** calculate matrix for external circuit  cmat ************   

c                                                                       

           subroutine cmat                                              

c                                                                       

       include 'common'

c                                                                       

c ********* if double precision ******************                      

c          implicit real*8(d)                                           

c                                                                       

      do 100 i=1,nbr                                                    

        do 110 j=1,nbr                                                  

          dzz(i,j)=0.0                                                  

          dyy(i,j)=0.0                                                  

  110   continue                                                        

  100 continue                                                          

c                                                                       

      do 200 i=1, ncn                                                   

        do 210 j=1, nbr                                                 

          das(i,j)=0.0                                                  

  210   continue                                                        

  200 continue                                                          

c                                                                       

      do 220 i=1,nbr/2                                                  

        if (ntabl(i,3).eq.1) then                                       

          ii=2*i-1                                                      

          dzz(ii,ii)=-dtabl(i,1)                                        

          dzz(ii+1,ii+1)=-dtabl(i,1)                                    

        else if (ntabl(i,3).eq.2) then                                  

          ii=2*i-1                                                      

          dzz(ii,ii+1)=-dn(nhh,1,2)*dtabl(i,1)                          

          dzz(ii+1,ii)=-dn(nhh,2,1)*dtabl(i,1)                          

        else if (ntabl(i,3).eq.3) then                                  

          ii=2*i-1                                                      

          dzz(ii,ii)=-1.0                                               

          dzz(ii+1,ii+1)=-1.0                                           

        end if                                                          

  220 continue                                                          

c                                                                       

      do 230 i=1,nbr/2                                                  

        if (ntabl(i,3).ne.3) then                                       

          ii=2*i-1                                                      

          dyy(ii,ii)=1.0                                                

          dyy(ii+1,ii+1)=1.0                                            

        else                                                            

          ii=2*i-1                                                      

          dyy(ii,ii+1)=dn(nhh,1,2)*dtabl(i,1)                           

          dyy(ii+1,ii)=dn(nhh,2,1)*dtabl(i,1)                           

        end if                                                          

  230 continue                                                          

c                                                                       

      do 240 i=1,nbr/2                                                  

        if (ntabl(i,1).ne.0) then                                       

          ii=2*i-1                                                      

          jj=2*ntabl(i,1)-1                                             

          das(jj,ii)=1.0                                                

          das(jj+1,ii+1)=1.0                                            

        end if                                                          

        if (ntabl(i,2).ne.0) then                                       

          ii=2*i-1                                                      

          jj=2*ntabl(i,2)-1                                             

          das(jj,ii)=-1.0                                               

          das(jj+1,ii+1)=-1.0                                           

        end if                                                          

  240 continue                                                          

c                                                                       

c ***** make dhcl matrix ****                                           

c                                                                       

      do 300 i=1,nbr                                                    

        do 310 j=1,ncn                                                  

          dhcl(i,j)=0.0                                                 

          do 320 k=1,nbr                                                

            dhcl(i,j)=dhcl(i,j)+dyy(i,k)*das(j,k)                       

  320     continue                                                      

  310   continue                                                        

        do 330 j=1,nbr                                                  

          dhcl(i,ncn+j)=dzz(i,j)                                        

  330   continue                                                        

  300 continue                                                          

c                                                                       

      do 340 i=1,ncn                                                    

        do 350 j=1,nbr                                                  

          dhcl(nbr+i,ncn+j)=das(i,j)                                    

  350   continue                                                        

  340 continue                                                          

c                                                                       

       return                                                           

       end                                                              

c ****************************************************************      

c ***********  calculate reluctivity matrix  dmat  ***************      

c                                                                       

           subroutine dmat                                              

c                                                                       

       include 'common'

c                                                                       

      do 100 ne=1, nelem                                                

        if (int(nod(ne,4)/100).eq.2)         goto 200                   

c                                                                       

c -------  if material is air or coil  -----------                      

c                                                                       

        do 110 k=1,mmu                                                  

          if (k.eq.1) then                                              

            dd(ne,k,1,1)=1.0d0/dmu0                                     

            dd(ne,k,1,2)=0.0                                            

            dd(ne,k,2,1)=0.0                                            

            dd(ne,k,2,2)=1.0d0/dmu0                                     

          else                                                          

            dd(ne,k,1,1)=0.0                                            

            dd(ne,k,1,2)=0.0                                            

            dd(ne,k,2,1)=0.0                                            

            dd(ne,k,2,2)=0.0                                            

          end if                                                        

  110   continue                                                        

c                                                                       

        goto 100                                                        

c                                                                       

c -------  magnetic core  200<=nod(i,4)<=299  -------                   

c                                                                       

  200     continue                                                      

c                                                                       

c -------  fourier expansion of reluctivity  --------                   

c                                                                       

          call rfour(ne)                                                

c                                                                       

c -------  make the reluctivity matrix  -------------                   

c                                                                       

        kk=nhh                                                          

          kk1=nhh                                                       

        do 210 k=1,mmu                                                  

          if (kk1.le.mmu) then                                          

            dd(ne,k,1,1)=(dreluc(2*(kk-nhh),1)                          

     1                  - dreluc(2*(kk+nhh-1),1))*.5                    

            dd(ne,k,1,2)=(-dreluc(2*(kk-nhh),0)                         

     1                  + dreluc(2*(kk+nhh-1),0))*.5                    

            dd(ne,k,2,1)=(dreluc(2*(kk-nhh),0)                          

     1                  + dreluc(2*(kk+nhh-1),0))*.5                    

            dd(ne,k,2,2)=(dreluc(2*(kk-nhh),1)                          

     1                  + dreluc(2*(kk+nhh-1),1))*.5                    

          else                                                          

            dd(ne,k,1,1)=(dreluc(2*(nhh-kk),1)                          

     1                  - dreluc(2*(nhh+kk-1),1))*.5                    

            dd(ne,k,2,1)=(-dreluc(2*(nhh-kk),0)                         

     1                  + dreluc(2*(nhh+kk-1),0))*.5                    

            dd(ne,k,1,2)=(dreluc(2*(nhh-kk),0)                          

     1                  + dreluc(2*(nhh+kk-1),0))*.5                    

            dd(ne,k,2,2)=(dreluc(2*(nhh-kk),1)                          

     1                  + dreluc(2*(nhh+kk-1),1))*.5                    

         end if                                                         

           kk=kk+1                                                      

           kk1=kk1+1                                                    

           if (kk.gt.mmu)    then                                       

                kk=kk-mmu                                               

           end if                                                       

  210   continue                                                        

c                                                                       

  100  continue                                                         

       return                                                           

       end                                                              
c *******************************************************               

c **************** pri-process pripro ******************                

c                                                                       

         subroutine pripro                                              

c                                                                       

       include 'common'

c                                                                       

c ********* if double precision ******************                      

c          implicit real*8(d)                                           

c                                                                       

c ---------- set harmonic matrix dn( , , ) ----------                   

       do 10 j=1,nmd                                                    

         dn(j,1,1)=0.0                                                  

         dn(j,1,2)=-2.0d0*float(j)+1.0                                  

         dn(j,2,1)= 2.0d0*float(j)-1.0                                  

         dn(j,2,2)=0.0                                                  

   10  continue                                                         

c ---------- clear matrix ----------                                    

      do 20 j=1, nmd                                                    

        do 21 k=1,nom                                                   

          dkk(j,k)=0.0                                                  

   21   continue                                                        

   20 continue                                                          

c                                                                       

      ddh=dpi/float(ndi)                                                

      do 30 n=0,2*(2*nmd-1)                                             

        do 31 nt=0,ndi                                                  

          dsine(n,nt)=sin(float(n*nt)*ddh)                              

          dcosi(n,nt)=cos(float(n*nt)*ddh)                              

   31   continue                                                        

   30 continue                                                          

c                                                                       

      return                                                            

      end                                                               

c ******************************************                            

c ********** clear matrix  init ************                            

c                                                                       

           subroutine init                                              

c                                                                       

       include 'common'

c                                                                       

c ********* if double precision ******************                      

c          implicit real*8(d)                                           

c                                                                       

       do 100 i=1, nom1                                                 

         do 110 j=1, nbb                                                

           dh(i,j) =0.0                                                 

  110    continue                                                       

         do 120 j=1, nom2                                               

           dhc(i,j)=0.0                                                 

           dhl(j,i)=0.0                                                 

  120    continue                                                       

         dk(i)=0.0                                                      

  100  continue                                                         

c                                                                       

       do 200 i=1,nom2                                                  

         do 210 j=1,nom2                                                

           dhcl(i,j)=0.0                                                

  210    continue                                                       

         dk(nom1+i)=0.0                                                 

  200  continue                                                         

c                                                                       

       return                                                           

       end                                                              
c **************************************************************        

c ***********  calculate other data  postpro  ******************        

c                                                                       

           subroutine postpr                                            

c                                                                       

       include 'common'

c                                                                       

c 1. calculate flux density                                             

c                                                                       

           call flux                                                    

c                                                                       

c 3. other                                                              

c                                                                       

           return                                                       

           end                                                          
