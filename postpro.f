c **************************************************************        
c ***********  calculate other data  postpro  ******************        
c                                                                       
           subroutine postpr                                            
c                                                                       
       include 'common.inc'
c                                                                       
c   calculate flux density                                             
c                                                                       
           call flux                                                    
c                                                                                                                                    
c                                                                       
           return                                                       

           end                                                          

c **************************************************************        
c ***************  save and output data  odata  ****************        
c                                                                       
           subroutine odata                                             
c                                                                       
       include 'common.inc'
c                                                                       
c ********* if double precision ******************                      
c          implicit real*8(d)                                           

c                                                                       
c check 9               **** open write file ****                       
c                                                                       
           open (30,file='ab0169.hbfem03.data(oda4000)')                
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
       do  i=1,nelem                                                 
         write (30, '(i4,6e12.4)') i,(dbx(i,k,1), k=1,mmu)              
         write (30, '(4x,6e12.4)')   (dbx(i,k,2), k=1,mmu)              
         write (30, '(4x,6e12.4)')   (dby(i,k,1), k=1,mmu)              
         write (30, '(4x,6e12.4)')   (dby(i,k,2), k=1,mmu)              
       enddo                                                         
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