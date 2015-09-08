c *******************************************************               
c ***************** read data   rdata *******************               
c                                                                       
         subroutine rdata                                               
c                                                                       
       include 'common.inc'
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

c*        print *, char                                                  
        do i=1, npo3, 2                                              

c       do 10 i=1, npo3, 5                                              

c         read(10,*) (xy(i+j,1), xy(i+j,2), j=0, 5)                     

          read(10,*) (xy(i+j,1), xy(i+j,2), j=0, 1)                     

        enddo                                                     

c                                                                       
c ------  change unit from mm to m  ------                              
c                                                                       
c check3                                                                

          do i=1, npo3                                               
            xy(i,1)=xy(i,1)*.001                                        
            xy(i,2)=xy(i,2)*.001                                        
          enddo                                                      
c                                                                       
c ------  input node of elements  nod( , )  -----                       
c                                                                       
        read(10,1000) char                                              

c*        print *, char                                                 


        do 30 i=1, nelem, 4                                             
c       do 30 i=1, nelem, 2                                             
c         read(10,*) (nod(i+j,1), nod(i+j,2),                           
c    &                nod(i+j,3), nod(i+j,4),j=0,1)                     

          read(10,*) (nod(i+j,1), nod(i+j,2),                           
     &                nod(i+j,3), nod(i+j,4),j=0,3)                     
        enddo                                                     

c                                                                       
c ------  input magnetizing voltages  ------                            
c                                                                       
        read(10,1000) char                                              

        read(10,*) ncv                                                  

        do 40 kk=1,ncv                                                  
          read(10,*) ncoil                                              
          ncoil=ncoil-299                                               
          read(10,*) (dv(ncoil,j,1),dv(ncoil,j,2),j=1,nmd)              
        enddo                                                       
c                                                                       
c  ---------- parameters of windings --------------                     
c                                                                       
        read (10,1000) char                                             

c*        print *, char                                                 
        read(10,*) ncc                                                  

        if(ncc.eq.0)          go to 140                                 
		
        do kk=1,ncc                                                  
          read(10,*) ncoil                                              
          ncoil=ncoil-299                                               
          read (10,*) dnds(ncoil,1),dnds(ncoil,2),dnds(ncoil,3)         
        enddo                                                      

  140  continue                                                         

c                                                                       
c ------  input magnetizing frequency  --------------                   
c                                                                       
        read(10,1000) char                                              

c*        print *, char                                                  

        read(10,*) df                                                   
        domeg=2.d0*dpi*df                                               
c                                                                       
c ------  input conductivity at node( ,4)=400  -----                    
c                                                                       
        read(10,1000) char                                              

c*        print *, char                                                  

        read(10,*) nre                                                  
		
        if(nre.eq.0)            go to 130                               

        do kk=1,nre                                                  
          read(10,*) dre(kk)                                            
        enddo                                                        
c                                                                       
  130 continue                                                          

c ------  input periodic boundary condition  -------                    
c                                                                       

        read(10,1000) char                                              
        if(npor.eq.0)           go to 110                               

        do kk=1,npor                                                 
          read(10,*) (nodr(kk,j),j=1,2)                                 
        enddo                                                        
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

       do i=1,nbr                                                    
         read(10,*) ndummy,(ntabl(i,j),j=1,3),dtabl(i,1)                 
       enddo                                                          

c                                                                       
       do i=1, nbr                                                   
         if ((ntabl(i,3).eq.2).or.(ntabl(i,3).eq.3)) then                
           dtabl(i,1)=dtabl(i,1)*domeg                                   
         end if                                                            
       enddo                                                         
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
         do l=1,nmd                                                  
           daa(l,k)  =0.0                                               
           daa(l,k+1)=0.0                                               
        enddo                                                  
       enddo                                                        

       go to 300                                                        
c                                                                       
  120  continue                                                         

c                                                                       
c -------  read daa( , )  ------------                                  
c                                                                       
c check 22                                                              
c                            ********* open read file ********          
      open (20,file='ab0169.hbfem03.data(femref)')                      
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

c*        read(20,1000) char                                             

c*        print '(a72)',char                                             

c                                                                       
  300 continue                                                          
c                                                                       
      return                                                            
c                                                                       

  999 print '('' error exits in reading initial values (sub.rdata) '')' 
      return                                                            

      end                                                