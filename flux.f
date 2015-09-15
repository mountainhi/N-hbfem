c ******  calculate flux density for each element  flux  ******         
c                                                                       
           subroutine flux                                              
c                                                                       
       include 'common.inc'
c   
c ***** definition of node number *****                                 
c 1---(unknown npo1)---(dirichlet npo2)---(periodic npo3)  
                         
C      daa() : initial values of vector potential   
                                       
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
              print '('' error in data of element, sub(flux)'',2i3)',  
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
        do  j=1,mmu                                                  
          dbx(ne,j,1)=0.0                                               
          dbx(ne,j,2)=0.0                                               
          dby(ne,j,1)=0.0                                               
          dby(ne,j,2)=0.0                                               
c                                                                       
             do k=1,3                                                  
               dbx(ne,j,1)=dbx(ne,j,1)+dr(k)*dda(j,k,1)                    
               dbx(ne,j,2)=dbx(ne,j,2)+dr(k)*dda(j,k,2)                    
               dby(ne,j,1)=dby(ne,j,1)-dq(k)*dda(j,k,1)                    
               dby(ne,j,2)=dby(ne,j,2)-dq(k)*dda(j,k,2)                    
             enddo                                                     

          dbx(ne,j,1)=dbx(ne,j,1)*.5/dcs(ne)                            
          dbx(ne,j,2)=dbx(ne,j,2)*.5/dcs(ne)                            
          dby(ne,j,1)=dby(ne,j,1)*.5/dcs(ne)                            
          dby(ne,j,2)=dby(ne,j,2)*.5/dcs(ne)                            

        enddo         
		
c.....loop over elements                                                                       
  100 continue                                                          

      return                                                            

      end                                                               
