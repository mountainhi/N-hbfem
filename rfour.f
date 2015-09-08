c **************************************************************        
c *********  fourier expansion of reluctivity, rfour ***********        
c                                                                       
c            ne : number of element                                     
c                                                                       
           subroutine rfour(ne)                                         
c                                                                       
           include 'common.inc'
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

         do i=1, mmu                                                
           dbbx = dbbx+dbx(ne,i,1)*dsine(2*i-1,nt)                      
     &                +dbx(ne,i,2)*dcosi(2*i-1,nt)                      
           dbby = dbby+dby(ne,i,1)*dsine(2*i-1,nt)                      
     &                +dby(ne,i,2)*dcosi(2*i-1,nt)                      
         enddo                                                      

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

          do i=1,mmu                                                
          dbbxt = dbbxt+(dbx(ne,i,1)*dcosi(2*i-1,nt)                    
     &                 - dbx(ne,i,2)*dsine(2*i-1,nt))*(2*i-1)*domeg     
          dbbyt = dbbyt+(dby(ne,i,1)*dcosi(2*i-1,nt)                    
     &                 - dby(ne,i,2)*dsine(2*i-1,nt))*(2*i-1)*domeg     
          enddo                                                    
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

         do n=1,n1,2                                                
           dsc=dsc+4.*dab(n)+2.*dab(n+1)                                
         enddo                                                     
c   **** 2*dreluc(0,1) ****                                             
         dreluc(0,1)=2.*(dsc-dab(ndi))*ddh*dcon/3.                      
         dreluc(0,0)=0.0                                                
c                                                                       
c check 8                                                               
c                                                                       
        nend=2*(mmu+nhh-1)                                              

         do n=2, nend, 2                                            
           dsc=dab(0)                                                   
           dss=0.0                                                      
c                                                                       
            do nk=1, n1, 2                                           
              dys=dab(nk)*dsine(n,nk)                                    
                dss=dss+4.*dys                                           
              dys=dab(nk+1)*dsine(n,nk+1)                                
                dss=dss+2.*dys                                           
              dyc=dab(nk)*dcosi(n,nk)                                    
                dsc=dsc+4.*dyc                                           
              dyc=dab(nk+1)*dcosi(n,nk+1)                                
                dsc=dsc+2.*dyc                                           
            enddo                                                   

           dreluc(n,0)=(dss-dys)*ddh*dcon*2.0/3.0                       
           dreluc(n,1)=(dsc-dyc)*ddh*dcon*2.0/3.0                       

         enddo                                                      
c                                                                       
       return                                                           

       end         