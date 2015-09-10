c ***************************************************************       
c ************  assemble the system matrix  fem  *****************       
c                                                                       
           subroutine fem                                               
c                                                                       
       include 'common.inc'
c                                                                       
       dimension ddd(nmd,2,2), dss(3,3)
		 
c
        do i=1, nom2                                                  
          dk(nom1+i)=dkk(nhh,nom1+i)                                      
        enddo    		 
c                                                                       
c **********  loop for element  ****************                        
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

          do  j=1,2                                                  
            dhc(nlow+j,jj1+j)=dhc(nlow+j,jj1+j)-dg(ne)*dsig1            
          enddo                                                     
c                                                                       
          do  j=1,2                                                  
            do 150 k=1,2                                                
c              if(dn(nhh,j,k).eq.0.0)                    go to 150    
            if(dn(nhh,j,k).ne.0.0)   
              dhl(jj2+j,nlow+k)=dhl(jj2+j,nlow+k)                       
     &                         -dn(nhh,j,k)*dc(ne)*dsig1*2.             
            enddo                                                  
           enddo                                                     

  110   continue                                                        
c                                                                       
  190   continue                                                        
c                                                                       
c ---------  set ddd(mmu,2,2)  ----------                               
c                                                                       
       do i=1, 2                                                    
         do j=1, 2                                                  
           do k=1, mmu                                              

             ddd(k,i,j)=dd(ne,k,i,j)                                    

           enddo                                                     
         enddo                                                      
       enddo                                                          
c                                                                       
c --------  dss( , )  -------------                                     
c                                                                       
        nss=1                                                           

        do 240 j=1,3                                                    
          do 250 k=j,3                                                  

            dss(j,k)=ds(ne,nss)                                         
            dss(k,j)=ds(ne,nss)                                         
            nss=nss+1         
			
          enddo                                                     
        enddo                                                     
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

            do  m=1, 2                                               
              nl=nlow+m-1                                               
              do  n=1, 2                                             
                nc =ncol+n-1                                            
                nc1=nc-nl+nb                                            
c                                                                       
                dh(nl,nc1) =dh(nl,nc1)  + ddd(1,m,n)*dnc1               
     &                                  + dn(nhh,m,n)*ds(ne,7)*dnc2     

                do l=2,mmu                                          
                  dk(nl)=dk(nl)-ddd(l,m,n)*dnc1*                        
     &                          daa(mod(l+nhh-2,mmu)+1,nc)              

                enddo
              enddo				
            enddo                                                    

  310     continue                                                      

  300   continue                                                        
C.....100 loop over elements
  100 continue                                                          
c                                                                       
      return                                                            

      end                                
	  
c *******************************************************************   
c ********** calculate matrix for each element  smat ****************   
c                                                                       
           subroutine smat                                              
c                                                                       
       include 'common.inc'
c                                                                       
       dimension dx(3),dy(3),dq(3),dr(3)                                
c                                                                       
c ********* if double precision ******************                      
c          implicit real*8(d)                                           
c                                                                       
c ----- calculate the matrix ds( , ) ------                             
      do 100 ne=1,nelem                                                 
        do j=1,3                                                    
          dx(j)=xy(nod(ne,j),1)                                         
          dy(j)=xy(nod(ne,j),2)                                         
        enddo                                                        

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
		
        do  j=1,3                                                    
          do  k=j,3                                                  
            ds(ne,n)=0.25d0/dcs(ne)*(dr(j)*dr(k)+dq(j)*dq(k))           
            n=n+1                                                       
          enddo                                                       
        enddo                                                         
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
        do i=1,nbr/2                                                
          do j=1,nmd                                                
            jj1=nom1+2*i-1                                              
            dkk(j,jj1)=dv(i,j,1)                                        
            dkk(j,jj1+1)=dv(i,j,2)                                      
          enddo                                                  
        enddo                                                         

        return                                                          

        end                                                           

c *******************************************************************   
c ********** calculate matrix for external circuit  cmat ************   
c                                                                       
           subroutine cmat                                              
c                                                                       
       include 'common.inc'
c                                                                       
c ********* if double precision ******************                      
c          implicit real*8(d)                                           
c                                                                       
      do i=1,nbr                                                    
        do j=1,nbr                                                  
          dzz(i,j)=0.0                                                  
          dyy(i,j)=0.0                                                  
        enddo                                                        
      enddo                                                         
c                                                                       
      do i=1, ncn                                                   
        do j=1, nbr                                                 
          das(i,j)=0.0                                                  
        enddo                                                        
      enddo                                                         
c                                                                       
      do i=1,nbr/2                                                  
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
      enddo                                                          
c                                                                       
      do i=1,nbr/2                                                  
        if (ntabl(i,3).ne.3) then                                       
          ii=2*i-1                                                      
          dyy(ii,ii)=1.0                                                
          dyy(ii+1,ii+1)=1.0                                            
        else                                                            
          ii=2*i-1                                                      
          dyy(ii,ii+1)=dn(nhh,1,2)*dtabl(i,1)                           
          dyy(ii+1,ii)=dn(nhh,2,1)*dtabl(i,1)                           
        end if                                                          
      enddo                                                          
c                                                                       
      do  i=1,nbr/2                                                  
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
      enddo                                                          
c                                                                       
c ***** make dhcl matrix ****                                           
c                                                                       
      do 300 i=1,nbr     
	  
        do j=1,ncn                                                  
          dhcl(i,j)=0.0                                                 
          do k=1,nbr                                                
            dhcl(i,j)=dhcl(i,j)+dyy(i,k)*das(j,k)                       
          enddo                                                     
        enddo

        do j=1,nbr                                                  
          dhcl(i,ncn+j)=dzz(i,j)                                        
        enddo                                                        

  300 continue                                                          

c                                                                       
      do  i=1,ncn                                                    
        do  j=1,nbr                                                  
          dhcl(nbr+i,ncn+j)=das(i,j)                                    
        enddo                                                         
      enddo                                                         
c                                                                       
       return                                                           
       end                                                              

c ****************************************************************      
c ***********  calculate reluctivity matrix  dmat  ***************      
c                                                                       
           subroutine dmat                                              
c                                                                       
       include 'common.inc'
c                                                                       
      do 100 ne=1, nelem                                                
        if (int(nod(ne,4)/100).eq.2)         goto 200                   
c                                                                       
c -------  if material is air or coil  -----------                      
c                                                                       
        do k=1,mmu                                                  
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
        enddo                                                         
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
c ....      nhh = nloop    (number of harmonics)                                                                
        kk  = nhh                                                          
        kk1 = nhh                                                       
		  
        do 210 k=1,mmu  
		
          if (kk1.le.mmu) then                                          
            dd(ne,k,1,1) = (dreluc(2*(kk-nhh),1)                          
     1                  - dreluc(2*(kk+nhh-1),1))*.5                    
            dd(ne,k,1,2) = (-dreluc(2*(kk-nhh),0)                         
     1                  + dreluc(2*(kk+nhh-1),0))*.5                    
            dd(ne,k,2,1) = (dreluc(2*(kk-nhh),0)                          
     1                  + dreluc(2*(kk+nhh-1),0))*.5                    
            dd(ne,k,2,2) = (dreluc(2*(kk-nhh),1)                          
     1                  + dreluc(2*(kk+nhh-1),1))*.5                    
          else                                                          
            dd(ne,k,1,1) = (dreluc(2*(nhh-kk),1)                          
     1                  - dreluc(2*(nhh+kk-1),1))*.5                    
            dd(ne,k,2,1) = (-dreluc(2*(nhh-kk),0)                         
     1                  + dreluc(2*(nhh+kk-1),0))*.5                    
            dd(ne,k,1,2) = (dreluc(2*(nhh-kk),0)                          
     1                  + dreluc(2*(nhh+kk-1),0))*.5                    
            dd(ne,k,2,2) = (dreluc(2*(nhh-kk),1)                          
     1                  + dreluc(2*(nhh+kk-1),1))*.5                    
         end if                                                         

           kk=kk+1                                                      
           kk1=kk1+1                                                    

           if (kk.gt.mmu)    then                                       
                kk=kk-mmu                                               
           end if                                                       

  210   continue                                                        
c......100 loop over elements                                                                       
  100  continue                                                         

       return                                                           

       end                                                              
c ******************************************************* 