c **************************************************************        
c ****************  check maximum harmonic   *******************        
c                                                                       
           subroutine check                                             
c                                                                       
         include 'common.inc'
c                                                                       
         nonc=0                                                         
         damax1=0.0                                                     
c                                                                       
         do i=1,nom1                                                
           if(abs(daa(1,i)).gt.damax1)    damax1=abs(daa(1,i))          
         enddo                                                       
c                                                                       
         damax2=0.0                                                     
c                                                                       
         do i=1,nom1                                                
           if(abs(daa(mmu,i)).gt.damax2)  damax2=abs(daa(mmu,i))        
         enddo                                                      
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