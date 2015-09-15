                 program fh0922                                       
c *************************************************************         
c *************************************************************         
c ***((( fort77 )))                                         ***         
c ***  time-periodic nonlinear magnetic field analysis      ***         
c ***    by harmonic balance finite element method          ***         

c ***          femh0922              ver.3.4                ***         
c ***                                                       ***         
c ***  1. generalized harmonic balance fem                  ***         
c ***  2. harmonic up to 11th can be considered             ***         
c ***  3. the maximum harmonics is automatically defined    ***         
c ***  4. method of block-structure iteration is employed   ***         
c ***  5. the external circuit is considered                ***         
c ***                                                       ***         
c ***                                                       ***         
c ***                                       91.9.10         ***   
c ***                                       15.7.31         ***        
c *************************************************************         
c *************************************************************         
c   main.f    (main program)                                          

c   prepro     prepro  init                                             

c   init.f     rdata                                               

c   fem.f      smat  dmat  cmat  fem                                    

c   gauss3.f   gauss3                                                   

c   flux.f     flux                                                     

c   rfour.f    rfour   

C   conv.f     conv check                                          

c   postpro.f   odata  postpr  ckeck                              

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

c    nod( , )       node number of triangle element                     

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
      include 'common.inc'
c                                                                       
c ************ if double precision *********************                
c        implicit  real*8(d)                                            
c                                                                       

c check1                                                                
c ************ set condition **********************                     
c                                                                       
c   1. maximum iteration            max                                 
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

c *********** set coefficients ************                              
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
              print '('' read input data--subroutine rdata '')'                                           

              call rdata                                                
c                                                                       
c 1-2 set coefficients and initial values                               
              print '('' set coefficients --prepro '')'                                          

              call prepro                                               
c                                                                       
c 1-3 calculate matrix for each element                                 
              print '('' smat '')'                                            

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
               print '(''stop--- Maximum iterations reached'',i4)',max      
               goto 9999                                                
         end if                                                         
c                                                                       
c 2-1 clear variables                                                   
          print '('' init '')'                                             

          call init                                                     
c                                                                       
c 2-2 calculate flux density                                            

          print '('' flux '')'                                             
           call flux                                                     
c                                                                       
c 2-3 calculate reluctivity matrix                                      
          print '('' dmat '')'                                             
          call dmat                                                     
c                                                                       
c 2-4-1 set system matrix                                               
          print '('' fem '',i2)',nloop                                     
          call fem                                                      
c                                                                       
c 2-4-2 set system matrix                                               
          print '('' cmat '',i2)',nloop                                    
          call cmat                                                     
c                                                                       
c 2-5 calculate system equation                                         
          print '('' gauss3 '',i2)',nloop                                  

          call gauss3                                                   

c                                                                       
c 2-6 store each harmonic data                                          
          print '('' dstore '',i2)',nloop                                  
          call dstore                                                   
c                                                                       
c 2-7 check convergence                                                 
          print '('' conv'',i2)', nloop                                    
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

      do  i=1, mmu                                                   
        if (ntotal(2,i).ne.1)     nonc=nonc+1                           
      enddo                                                       

c                                                                       
      do  i=1, mmu                                                   
        ntotal(2,i)=0                                                   
      enddo                                                          
c                                                                       
           if (nonc.ne.0)   goto 100                                    
c                                                                       
c ##### number of non-converged potential #####                        
c                                                                       
c 3-1 check  convergence                                                
          print '('' check  '')'                                           

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
             print '('' postpr '')'                                             
              call postpr                                               
c                                                                       
c 3-2 print out data                                                    
             print '('' odata '')'                                            
              call odata                                                

      stop                                                              

      end                          