c **************************************************************        
c ****************  check convergence  conv  *******************        
c                                                                       
           subroutine conv                                              
c                                                                       
       include 'common.inc'
c                                                                       
       dimension nlist(1000)                                            
c                                                                       
         nonc=0                                                         
         k=1                                                          
c                                                                       
         damax=0.0                                                      
c                                                                       
         do i=1,nom1                                                
           if(abs(daaa(1,i)).gt.damax)  damax=abs(daaa(1,i))            
         enddo                                                        
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
C        if (nhh.gt.1)                       go to 400                  
c                                                                       
c        do 300 i=nom1+1,nom                                            
c*          if (abs(daaa(nhh,i)).le.1.0e-25) then                        
c*          go to 300                                                    
c*           end if                                                      
c                                                                       
c*            dadv=10.                                                   
c*          derr=(daaa(nhh,i)-daa(nhh,i))/daaa(nhh,i)                    

c*           if (abs(derr).gt.aaa*dadv) then                             
c*               nonc=nonc+1                                             
c*               nlist(k)=(i+1)/2                                        
c*               k=k+1                                                   
c*           end if                                                      
c*  300    continue                                                      
c*  400    continue                                                      
c                                                                       

c -----  delerating process  -------                                    
c                                                                       
         if((ntotal(1,nhh)+itr).ge.100)                                 
     &                      bbbb=bbb*0.06            

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

         do i=1,nom1                                                
             daa(nhh,i)=daa(nhh,i)+bbbb*(daaa(nhh,i)-daa(nhh,i))        
         enddo                                                       
c                                                                       

         do i=nom1,nom                                              
             daa(nhh,i)=daa(nhh,i)+2.*bbbb*(daaa(nhh,i)-daa(nhh,i))     
         enddo                                                       
                                                                       
c                                                                       
c*           print '(20i4)', (nlist(kk),kk=1,k-1)                        
c                                                                       
         return                                                         

         end                                                            