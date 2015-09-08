c *******************************************************               
c **************** pre-process prepro ******************                
c                                                                       
         subroutine prepro                                              
c                                                                       
       include 'common.inc'
c                                                                       
c ********* if double precision ******************                      
c          implicit real*8(d)                                           
c                                                                       
c ---------- set harmonic matrix dn( , , ) ----------                   
       do  j=1,nmd                                                    
         dn(j,1,1)=0.0                                                  
         dn(j,1,2)=-2.0d0*float(j)+1.0                                  
         dn(j,2,1)= 2.0d0*float(j)-1.0                                  
         dn(j,2,2)=0.0                                                  
       enddo	   
c ---------- clear matrix ----------                                    
      do  j=1, nmd                                                    
        do  k=1,nom                                                   
          dkk(j,k)=0.0                                                  
        enddo                                                        
      enddo                                                         
c                                                                       
      ddh=dpi/float(ndi)       
	  
      do n=0,2*(2*nmd-1)                                             
        do nt=0,ndi                                                  
          dsine(n,nt)=sin(float(n*nt)*ddh)                              
          dcosi(n,nt)=cos(float(n*nt)*ddh)                              
        enddo                                                        
      enddo                                                           
c                                                                       
      return                                                            

      end                                               
	  
c ******************************************                            
c ********** initialize matrix  init ************                            
c                                                                       
           subroutine init                                              
c                                                                       
       include 'common.inc'
c                                                                       
c ********* if double precision ******************                      
c          implicit real*8(d)                                           
c                                                                       
       do i=1, nom1                                                 
         do  j=1, nbb                                                
           dh(i,j) =0.0                                                 
         enddo                                                       

         do  j=1, nom2                                               
           dhc(i,j)=0.0                                                 
           dhl(j,i)=0.0                                                 
         enddo                                                       

         dk(i)=0.0                                                      
       enddo                                                         
c                                                                       
       do i=1,nom2                                                  
         do j=1,nom2                                                
           dhcl(i,j)=0.0                                                
         enddo                                                         

         dk(nom1+i)=0.0                                                 

       enddo                                                           
c                                                                       
       return                                                           

       end                                                